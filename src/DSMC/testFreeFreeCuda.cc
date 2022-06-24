/*
  Test free-free photon generation
*/

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <tuple>

#include "atomic_constants.H"
#include "Ion.H"
#include "interactSelect.H"
#include <AsciiHisto.H>
#include <cxxopts.H>		// Option parsing
#include <libvars.H>		// EXP global variables

#include <mpi.h>

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

  int num, nhisto;
  double energy;
  unsigned short Z, C;

  cxxopts::Options options(av[0], "Free-free photon test routine for cuda");

  options.add_options()
    ("h,help", "produce help message")
    ("vanHoof", "Use Gaunt factor from van Hoof et al.")
    ("KandM", "Use original Koch & Motz cross section")
    ("N,Num", "number of evaluations",
     cxxopts::value<int>(num)->default_value("200"))
    ("n,Nhisto", "number of histogram bins",
     cxxopts::value<int>(nhisto)->default_value("20"))
    ("Z,ZZ", "atomic number",
     cxxopts::value<unsigned short>(Z)->default_value("1"))
    ("C,CC", "ionic charge",
     cxxopts::value<unsigned short>(C)->default_value("2"))
    ("e,energy", "minimum energy (Rydbergs)",
     cxxopts::value<double>(energy)->default_value("1.0"))
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

  if (vm.count("vanHoof")) {
    Ion::use_VAN_HOOF = true;
  }

  if (vm.count("KandM")) {
    Ion::use_VAN_HOOF = false;
  }

  std::string prefix("free_free_cuda");
  std::string cmdFile = prefix + ".cmd_line";
  std::ofstream out(cmdFile.c_str());
  if (!out) {
    std::cerr << "testFreeFreeCuda: error opening <" << cmdFile
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


				// Convert from Rydberg to eV
  const double ryd = 27.2113845/2.0;
  double eV = energy * ryd;

  std::cout << "# Q=(" << Z << ", " << C << ")" << std::endl
	    << "#" << std::endl;


  std::vector<double> xc, ph;
  ad.testFreeFreeGen(Z, C, eV, num, xc, ph);

  for (int i=0; i<num; i++)
    std::cout << std::setw(18) << xc[i]
	      << std::setw(18) << ph[i]
	      << std::endl;

  if (myid==0) {
    std::ofstream out(prefix + ".histo");
    if (out) {
      AsciiHisto spect(ph, nhisto, 0.01, true);
      spect(out);
    } else {
      std::cerr << "testFreeFree: error opening <" << prefix + ".histo" << ">"
		<< std::endl;
    }

    auto minmax = std::minmax_element(ph.begin(), ph.end());
    std::cout << "Min ph=" << *minmax.first  << std::endl;
    std::cout << "Max ph=" << *minmax.second << std::endl;
  }

  MPI_Finalize();

  return 0;
}
