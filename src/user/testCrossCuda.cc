/*
  Test cuda cross section implementation
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <tuple>

#include <boost/program_options.hpp>

#include "atomic_constants.H"
#include "Ion.H"
#include "Elastic.H"

namespace po = boost::program_options;

#include <mpi.h>

extern "C" void phfit2_(int* nz, int* ne, int* is, float* e, float* s);

// For sorting tuples
//
template<int M, template<typename> class F = std::less>
struct TupleCompare
{
  template<typename T>
  bool operator()(T const &t1, T const &t2)
  {
    return F<typename std::tuple_element<M, T>::type>()(std::get<M>(t1), std::get<M>(t2));
  }
};

int numprocs, myid;
std::string outdir(".");
std::string runtag("run");
char threading_on = 0;
pthread_mutex_t mem_lock;

int main (int ac, char **av)
{
  //===================
  // MPI preliminaries 
  //===================

  MPI_Init(&ac, &av);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  std::string cmd_line;
  for (int i=0; i<ac; i++) {
    cmd_line += av[i];
    cmd_line += " ";
  }

  unsigned short Z, C;
  double emin, emax, logL = 10.0, kdel;
  std::string RRtype;
  bool eVout = false;
  int num;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("eV",		"print results in eV")
    ("Z,Z",		po::value<unsigned short>(&Z)->default_value(2),
     "atomic number")
    ("C,C",		po::value<unsigned short>(&C)->default_value(3),
     "ionic charge")
    ("Emin,e",		po::value<double>(&emin)->default_value(0.001),
     "minimum energy (Rydbergs)")
    ("Emax,E",		po::value<double>(&emax)->default_value(100.0),
     "maximum energy (Rydbergs)")
    ("Num,N",		po::value<int>(&num)->default_value(200),
     "number of evaluations")
    ("kdel,k",          po::value<double>(&kdel)->default_value(0.25),
     "default logarithmic spacing for k grid")
    ("RRtype,R",	po::value<std::string>(&RRtype)->default_value("Verner"),
     "cross-section type")
    ;



  po::variables_map vm;

  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    return -1;
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    std::cout << "Example: Helium II recombination" << std::endl;
    std::cout << "\t" << av[0]
	      << " -Z 2 -C 2" << std::endl;
    MPI_Finalize();
    return 1;
  }

  if (vm.count("eV")) {
    eVout = true;
  }

  if (vm.count("kdel")) {
    Ion::kdel = kdel;
  }

  std::string prefix("crossSection");
  std::string cmdFile = prefix + ".cmd_line";
  std::ofstream out(cmdFile.c_str());
  if (!out) {
    std::cerr << "testCrossSection: error opening <" << cmdFile
	      << "> for writing" << std::endl;
  } else {
    out << cmd_line << std::endl;
  }

  // Initialize CHIANTI
  //

  std::set<unsigned short> ZList = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16};

  if (ZList.find(Z) == ZList.end()) {
    if (myid==0) {
      std::cout << "Z=" << Z 
		<< " is not in element list.  Current list contains:";
      for (auto z : ZList) std::cout << " " << z;
      std::cout << std::endl;
    }
    MPI_Finalize();
    exit(1);
  }

  chdata ch;

  ch.createIonList(ZList);

  return 0;
}
