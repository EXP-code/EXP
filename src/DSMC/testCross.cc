/*
  Test radiative cross section 
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
#include "Elastic.H"
#include <cxxopts.H>		// Option parsing
#include <libvars.H>		// EXP library globals

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

  unsigned short Z, C;
  double emin, emax, logL = 10.0, kdel;
  std::string RRtype;
  bool eVout = false;
  int num;

  cxxopts::Options options(av[0], "Test radiative cross section computation");

  options.add_options()
    ("h,help",   "produce help message")
    ("eV",       "print results in eV")
    ("vanHoof",  "Use Gaunt factor from van Hoof et al.")
    ("KandM",    "Use original Koch & Motz cross section")
    ("Z,ZZ",     "atomic number",
     cxxopts::value<unsigned short>(Z)->default_value("2"))
    ("C,CC",     "ionic charge",
     cxxopts::value<unsigned short>(C)->default_value("3"))
    ("e,Emin",   "minimum energy (Rydbergs)",
     cxxopts::value<double>(emin)->default_value("0.001"))
    ("E,Emax",   "maximum energy (Rydbergs)",
     cxxopts::value<double>(emax)->default_value("100.0"))
    ("N,Num",    "number of evaluations",
     cxxopts::value<int>(num)->default_value("200"))
    ("k,kdel",   "default logarithmic spacing for k grid",
     cxxopts::value<double>(kdel)->default_value("0.25"))
    ("R,RRtype", "cross-section type",
     cxxopts::value<std::string>(RRtype)->default_value("Badnell"))
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

  if (vm.count("eV")) {
    eVout = true;
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

  const double mrat = me/mp;

  PeriodicTable pt;

  atomicData ad;

  ad.createIonList(ZList);

  if (myid) {
    MPI_Finalize();
    return 0;
  }

  Ion::setRRtype(RRtype);

  Geometric geometric;
  Elastic   elastic;


  // Print cross section for requested ion
  //
				// Bohr cross section (pi*a_0^2) in nm
  const double b_cross = 0.00879735542978;

				// Convert from Rydberg to eV
  const double ryd      = 27.2113845/2.0;
  const double eVerg    = 1.60217733e-12;

  emin = log(emin);
  emax = log(emax);

  double dE = (emax - emin)/(num - 1);

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

  std::cout << std::left
	    << std::setw(16) << "# Energy"
	    << std::setw(16) << "Wave number"
	    << std::setw(16) << "Geometric"
	    << std::setw(16) << "Elastic"
	    << std::setw(16) << "Rutherford"
	    << std::setw(16) << "Free-free"
	    << std::setw(16) << "Collisional"
	    << std::setw(16) << "Ionization"
	    << std::setw(16) << "Recombination"
	    << std::setw(16) << "Photoionization"
	    << std::setw(16) << "Photo states"
	    << std::endl;

  std::cout << std::setw(16) << "# [1]";
  for (int i=1; i<nf; i++) {
    std::ostringstream sout; sout << " [" << i+1  << "]";
    std::cout << std::setw(16) << sout.str();
  }
  std::cout << std::endl << std::setw(16) << "# --------";
  for (int i=1; i<nf; i++) {
    std::cout << std::setw(16) << "---------";
  }
  std::cout << std::endl;

  for (int i=0; i<num; i++) {
    double E   = exp(emin + dE*i);
    double EeV = E * ryd;

    double geom = geometric(Z);
    double elas = elastic(Z, EeV);
    double coul = 0.0;

    {
      const double ips = 1000.0;
      double b90 = 0.5*esu*esu*(C-1) / (EeV * eVerg) * 1.0e7; // nm
      b90 = std::min<double>(b90, ips);
      coul = M_PI*b90*b90 * 4.0*mrat/pt[Z]->weight() * logL;
    }

    std::pair<double, double> ffre = ad.IonList[Q]->freeFreeCross(EeV, 0);

    double ionz = ad.IonList[Q]->directIonCross(EeV, 0);

    Ion::collType       CE1 = ad.IonList[Q]->collExciteCross(EeV, 0);
    
    std::vector<double> RE1 = ad.IonList[Q]->radRecombCross(EeV, 0);

    std::vector<double> PI1 = ad.IonList[Q]->photoIonizationCross(EeV, 0);

    double sum = ad.VernerXC.cross(Q, EeV);

    int ZZ  = Z;
    int Nel = Z - C + 1;
    float ee = EeV, cs, csum=0.0;

    for (int S=0; S<=7; S++) {
      phfit2_(&ZZ, &Nel, &S, &ee, &cs);
      if (S>0) csum += cs;
    }

    double coll = CE1.back().first;

    std::cout << std::setw(16) << (eVout ? EeV : E)
	      << std::setw(16) << 0.0001239841842144513*1.0e8/EeV
	      << std::setw(16) << geom       * scale
	      << std::setw(16) << elas       * scale
	      << std::setw(16) << coul       * scale
	      << std::setw(16) << ffre.first * scale
	      << std::setw(16) << coll       * scale
	      << std::setw(16) << ionz       * scale
	      << std::setw(16) << RE1.back() * scale
	      << std::setw(16) << PI1.back() * scale
	      << std::setw(16) << csum
	      << std::setw(16) << sum
	      << std::endl;
  }

  MPI_Finalize();

  return 0;
}
