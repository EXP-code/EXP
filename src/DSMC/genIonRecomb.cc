/*
  Test radiative cross section 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <tuple>

#include <boost/random/mersenne_twister.hpp>
#include <boost/program_options.hpp>

#include "Ion.H"

namespace po = boost::program_options;

#include <mpi.h>

// Global variables (not all used but needed for EXP libraries)

int numprocs, myid;
std::string outdir(".");
std::string runtag("run");
char threading_on = 0;
pthread_mutex_t mem_lock;
boost::mt19937 random_gen;

int main (int ac, char **av)
{
  //===================
  // MPI preliminaries 
  //===================

  MPI_Init(&ac, &av);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  unsigned short Z;
  double T;
  std::string RRtype;
  int norder;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("Z,Z",		po::value<unsigned short>(&Z)->default_value(2),
     "atomic number")
    ("Temp,T",		po::value<double>(&T),
     "temperature")
    ("norder,n",	po::value<int>(&norder)->default_value(20),
     "Laguerre order")
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
    std::cout << "Example: Helium ionization fractions at T=40,000 K" << std::endl;
    std::cout << "\t" << av[0]
	      << " -Z 2 -T 40000" << std::endl;
    MPI_Finalize();
    return 1;
  }

  std::string prefix("IonRecombFrac");
  std::ofstream out(prefix + ".data");
  if (out.bad()) {
    std::cerr << "testCrossSection: error opening <" << prefix + ".data"
	      << "> for writing" << std::endl;
    exit(-1);
  }

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

  PeriodicTable pt;

  Ion::setRRtype(RRtype);

  atomicData ad;

  ad.createIonList(ZList);

  if (myid) {
    MPI_Finalize();
    return 0;
  }

  std::map<unsigned short, double> values = ad.fraction(Z, T, norder);
    
  for (auto v : values) out << std::setw(20) << std::setprecision(12) << v.second;
  out << std::endl;
  
  MPI_Finalize();

  return 0;
}
