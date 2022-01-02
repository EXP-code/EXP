/*
  Test cuda cross section implementation
*/

/* Manual compile string:

mpiCC -g -o testCrossCuda testCrossCuda.o Ion.o cudaIon.o TopBase.o spline.o phfit2.o -lexputil -lexpgnu -lvtkCommonCore-7.1 -lvtkCommonDataModel-7.1 -lvtkIOXML-7.1 -lmpi -lcuda -lcudart

*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <tuple>

#include "atomic_constants.H"
#include "Ion.H"
#include "Elastic.H"
#include <cxxopts.H>		// Option parsing
#include <libvars.H>		// EXP global variables

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

  int num;
  double emin, emax;
  bool eVout = false;
  std::string scaling;

  cxxopts::Options options(av[0], "Cross-section test routine for cuda");

  options.add_options()
    ("h,help", "produce help message")
    ("vanHoof", "Use Gaunt factor from van Hoof et al.")
    ("KandM", "Use original Koch & Motz cross section")
    ("eV", "print results in eV")
    ("c,compare", "for comparison with CPU version")
    ("N,Num", "number of evaluations",
     cxxopts::value<int>(num)->default_value("200"))
    ("e,Emin", "minimum energy (Rydbergs)",
     cxxopts::value<double>(emin)->default_value("0.001"))
    ("E,Emax", "maximum energy (Rydbergs)",
     cxxopts::value<double>(emax)->default_value("100.0"))
    ("S,scaling", "cross-section scaling (born, mbarn, (null))",
     cxxopts::value<std::string>(scaling))
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
    MPI_Finalize();
    return 1;
  }

  if (vm.count("eV")) {
    eVout = true;
  }

  if (vm.count("vanHoof")) {
    Ion::use_VAN_HOOF = true;
  }

  if (vm.count("KandM")) {
    Ion::use_VAN_HOOF = false;
  }

  std::string prefix("crossCuda");
  std::string cmdFile = prefix + ".cmd_line";
  std::ofstream out(cmdFile.c_str());
  if (!out) {
    std::cerr << "testCrossCuda: error opening <" << cmdFile
	      << "> for writing" << std::endl;
  } else {
    out << cmd_line << std::endl;
  }

  // Initialize CHIANTI
  //
  // std::set<unsigned short> ZList = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16};

  std::set<unsigned short> ZList = {1, 2};

  atomicData ad;

  // Using CUDA-----------+
  //                      |
  //                      v
  ad.createIonList(ZList, true);

  std::cout << "# Ions = " << ad.IonList.size() << std::endl;
  if (vm.count("compare"))
    ad.testCrossCompare(num, emin, emax, eVout, scaling);
  else
    ad.testCross(num);
  
  MPI_Finalize();

  return 0;
}
