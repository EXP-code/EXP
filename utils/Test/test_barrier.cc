/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Test mpi barrier wrapper
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 02/05/04
 *
 ***************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <memory>
#include <random>

#include <cstdlib>
#include <cmath>

#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <BarrierWrapper.H>

namespace po = boost::program_options;

int myid=0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;
boost::mt19937 random_gen;

using namespace std;

int main(int argc, char **argv)
{
  //--------------------------------------------------
  // MPI preliminaries
  //--------------------------------------------------

  int numprocs, slaves, myid, proc_namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

  //--------------------------------------------------
  // Declare the supported options.
  //--------------------------------------------------

  bool barrier_check, barrier_light, barrier_quiet;
  bool barrier_extra, barrier_debug, barrier_label;
  int  times, vsize, ndice;

  po::variables_map vm;

  po::options_description desc("Available options");
  desc.add_options()
    ("help,h",
     "Produce help message")
    ("mistake,m",
     "Make an intentional mistake to break synchronization")
    ("broken,b",
     "Intentionally change label to test synchronization error")
    ("times,n", po::value<int >(&times)->default_value(10),            "number of iterations")
    ("dice,d",  po::value<int >(&ndice)->default_value(1),             "number of dice to roll")
    ("size,s",  po::value<int >(&vsize)->default_value(1000),          "data size")
    ("check,C", po::value<bool>(&barrier_check)->default_value(true),  "check the barrier")
    ("label,l", po::value<bool>(&barrier_label)->default_value(true),  "use labeling")
    ("light,L", po::value<bool>(&barrier_light)->default_value(false), "light-weight barrier")
    ("extra,V", po::value<bool>(&barrier_extra)->default_value(false),  "extra verbose")
    ("debug,D", po::value<bool>(&barrier_debug)->default_value(false),  "turn on debugging")
    ;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch(boost::program_options::error& e){
    std::cerr << "Invalid_option_value exception thrown parsing config file:"
	      << std::endl << e.what() << std::endl;
    return 2;
  } catch(std::exception e){
    std::cerr <<"Exception thrown parsing config file:" 
	      << std::endl << e.what() << std::endl;
    return 2;
  }

  if (vm.count("help")) {
    if (myid==0) {
      std::cout << std::setfill('-') << setw(76) << '-' 
		<< std::endl << std::setfill(' ')
		<< "Barrier wrapper test routine" << std::endl
		<< std::endl
		<< std::setfill('-') << setw(76) << '-' 
		<< std::endl << std::setfill(' ')
		<< "Usage: " << argv[0] << " [options] file" << std::endl
		<< std::setfill('-') << setw(76) << '-' 
		<< std::endl << std::setfill(' ')
		<< desc << std::endl;
    }
    
    MPI_Finalize();
    return 0;
  }

  bool mistake = false;
  if (vm.count("mistake")) {
    if (myid==1) {
      mistake = true;
      std::cout << "I will intentionally make a mistake at Iteration "
		<< times/2 << " which will cause an endless error loop"
		<< std::endl;
    }
  }

  bool badlabel = false;
  if (vm.count("broken")) {
    if (myid==0) {
      badlabel = true;
      std::cout << "I will intentionally change the label at Iteration "
		<< times-1 << " which will cause an endless error loop"
		<< std::endl;
    }
  }

  //--------------------------------------------------
  // Do stuff
  //--------------------------------------------------

  typedef std::shared_ptr<BarrierWrapper> BWptr;
  BWptr barrier(new BarrierWrapper(MPI_COMM_WORLD, barrier_label));

  if (barrier_check) barrier->on();
  else               barrier->off();
  if (barrier_light) barrier->setLightWeight();
  else               barrier->setHeavyWeight();
  if (barrier_quiet) BarrierWrapper::verbose       = false;
  else               BarrierWrapper::verbose       = true;
  if (barrier_extra) BarrierWrapper::extra_verbose = true;
  else               BarrierWrapper::extra_verbose = false;
  if (barrier_debug) BarrierWrapper::debugging     = true;
  else               BarrierWrapper::debugging     = false;

  std::default_random_engine generator(1+myid);
  std::uniform_int_distribution<int> distribution(1,6);
  auto dice = std::bind ( distribution, generator );

  for (int i=0; i<times; i++) {

    // Call the barrier wrapper
    //
    {
      std::ostringstream sout;
      sout << "In loop BEG, count = " << i+1;
      (*barrier)(sout.str(), __FILE__, __LINE__);
    }

    std::vector<double> vec1(vsize), vec0(vsize);
    double x = 0.0, dx = 1.1;
    for (auto & v : vec1) { v =  x; x += dx; }

    // Do an MPI thing
    //
    if (!mistake or i!=times/2)
      MPI_Reduce(&vec1[0], &vec0[0], vsize, MPI_DOUBLE, MPI_SUM, 
		 0, MPI_COMM_WORLD);
    else {
      std::cout << "Node " << myid << " is skipping MPI_Reduce" << std::endl;
    }

    // Wait a random amount of time determined by throwing ndice
    //
    int secs = 0;
    for (int j=0; j<ndice; j++) secs += dice();
    sleep(secs);

    // Call the barrier wrapper
    //
    {
      std::ostringstream sout;
      sout << "In loop END, count = " << i+1;
      if (badlabel and i+1==times and myid==0) sout << " ===>oops<===";
      (*barrier)(sout.str(), __FILE__, __LINE__);
    }
  }
    
  //--------------------------------------------------
  // MPI cleanup
  //--------------------------------------------------

  MPI_Finalize();

  return 0;
}

