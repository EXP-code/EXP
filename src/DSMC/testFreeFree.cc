/*
  Test free-free photon generation
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <memory>
#include <numeric>
#include <tuple>

#include "atomic_constants.H"
#include "Ion.H"
#include "interactSelect.H"
#include <AsciiHisto.H>
#include <cxxopts.H>		// Option parsing
#include <libvars.H>		// EXP library globals

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

  unsigned short Z, C;
  double energy, logL = 10.0, kdel;
  std::string RRtype;
  int num, nhisto;

  cxxopts::Options options(av[0], "Free-free photo test routine");

  options.add_options()
    ("h,help", "produce help message")
    ("eV", "print results in eV")
    ("vanHoof", "Use Gaunt factor from van Hoof et al.")
    ("KandM", "Use original Koch & Motz cross section")
    ("Z,ZZ", "atomic number",
     cxxopts::value<unsigned short>(Z)->default_value("1"))
    ("C,CC", "ionic charge",
     cxxopts::value<unsigned short>(C)->default_value("2"))
    ("e,energy", "minimum energy (Rydbergs)",
     cxxopts::value<double>(energy)->default_value("1.0"))
    ("N,Num", "number of evaluations",
     cxxopts::value<int>(num)->default_value("200"))
    ("n,Nhisto", "number of histogram bins",
     cxxopts::value<int>(nhisto)->default_value("20"))
    ("k,kdel", "default logarithmic spacing for k grid",
     cxxopts::value<double>(kdel)->default_value("0.25"))
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
    std::cout << "Example: Helium II recombination" << std::endl;
    std::cout << "\t" << av[0]
	      << " -Z 2 -C 2" << std::endl;
    MPI_Finalize();
    return 1;
  }

  if (vm.count("kdel")) {
    Ion::kdel = kdel;
  }

  if (vm.count("vanHoof")) {
    Ion::use_VAN_HOOF = true;
  }

  if (vm.count("KandM")) {
    Ion::use_VAN_HOOF = false;
  }

  std::string prefix("free_free");
  std::string cmdFile = prefix + ".cmd_line";
  std::ofstream out(cmdFile.c_str());
  if (!out) {
    std::cerr << "testFreeFree: error opening <" << cmdFile
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

  const double mrat = me/mp;

  PeriodicTable pt;

  atomicData ad;

  ad.createIonList(ZList);

  if (myid) {
    MPI_Finalize();
    return 0;
  }

  // Print cross section for requested ion
  //
				// Bohr cross section (pi*a_0^2) in nm
  const double b_cross = 0.00879735542978;

				// Convert from Rydberg to eV
  const double ryd      = 27.2113845/2.0;

  double scale = 1.0;

  if (vm.count("Born")) {
    scale = 1.0/b_cross;
  }

  if (vm.count("Megabarn")) {
    scale = 1.0e+04;
  }

  lQ Q(Z, C), QL(Z, C-1);

  int nf = 11;

  std::cout << "# Q=(" << Z << ", " << C << ")" << std::endl
	    << "#" << std::endl;

  double EeV = energy * ryd;

  
  InteractSelect IS;

  std::vector<double> ph;

  for (int i; i<num; i++) {
    auto ffre = ad.IonList[Q]->freeFreeCross(EeV, 0);
    ph.push_back(IS.selectFFInteract(ffre));
    std::cout << std::setw(18) << ffre.first * scale
	      << std::setw(18) << ph.back() << std::endl;
  }

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
