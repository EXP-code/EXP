/*
  Test radiative cross section 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <tuple>

#include "Ion.H"
#include <cxxopts.H>
#include <libvars.H>

#include <mpi.h>

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

  cxxopts::Options options(av[0], "Test radiative cross section");

  options.add_options()
    ("h,help", "produce help message")
    ("Z,Z", "atomic number",
     cxxopts::value<unsigned short>(Z)->default_value("2"))
    ("T,Temp", "temperature",
     cxxopts::value<double>(T))
    ("n,norder", "Laguerre order",
     cxxopts::value<int>(norder)->default_value("20"))
    ("R,RRtype", "cross-section type",
     cxxopts::value<std::string>(RRtype)->default_value("Verner"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
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
