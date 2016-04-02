/*
  Test radiative cross section 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <tuple>

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

  unsigned short Z;
  double Tmin, Tmax;
  std::string RRtype;
  int numT, norder;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("Z,Z",		po::value<unsigned short>(&Z)->default_value(2),
     "atomic number")
    ("Tmin,t",		po::value<double>(&Tmin)->default_value(1000.0),
     "minimum temperature")
    ("Tmax,T",		po::value<double>(&Tmax)->default_value(1000000.0),
     "maximum temperature")
    ("NumT,N",		po::value<int>(&numT)->default_value(200),
     "number of temperature points")
    ("norder,n",		po::value<int>(&norder)->default_value(20),
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
  std::string cmdFile = prefix + ".cmd_line";
  std::ofstream out(cmdFile.c_str());
  if (!out) {
    std::cerr << "testCrossSection: error opening <" << cmdFile
	      << "> for writing" << std::endl;
  } else {
    out << cmd_line << std::endl;
  }

  out.close();
  out.open(prefix + ".data", ios::out);
  if (out.bad()) {
    std::cerr << "testCrossSection: error opening <" << prefix + ".data"
	      << "> for writing" << std::endl;
    exit(-1);
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

  PeriodicTable pt;

  chdata ch;

  ch.createIonList(ZList);

  if (myid) {
    MPI_Finalize();
    return 0;
  }

  Ion::setRRtype(RRtype);

  double tmin = log(Tmin);
  double tmax = log(Tmax);
  double dT   = (tmax - tmin)/numT;

  for (int nt=0; nt<=numT; nt++) {

    double T = exp(tmin + dT*nt);

    std::map<unsigned short, double> values = ch.fraction(Z, T, norder);
    
    out << std::setw(16) << T;
    for (auto v : values) out << std::setw(16) << v.second;
    out << std::endl;
  }

  MPI_Finalize();

  return 0;
}
