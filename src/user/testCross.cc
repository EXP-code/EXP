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
  double emin, emax, logL = 10.0;
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

  constexpr double mrat = me/mp;

  PeriodicTable pt;

  chdata ch;

  ch.createIonList(ZList);

  if (myid) {
    MPI_Finalize();
    return 0;
  }

  Ion::setRRtype(RRtype);

  Geometric geometric;
  Elastic   elastic;


  // Print cross section for requested ion
  //
				// Convert from Rydberg to eV
  const double ryd      = 27.2113845/2.0;
  const double eVerg    = 1.60217733e-12;

  emin = log(emin);
  emax = log(emax);

  double dE = (emax - emin)/(num - 1);

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

    double ffre = ch.IonList[Q]->freeFreeCross(EeV, 0);

    double ionz = ch.IonList[Q]->directIonCross(EeV, 0);

    Ion::collType       CE1 = ch.IonList[Q]->collExciteCross(EeV, 0);
    
    std::vector<double> RE1 = ch.IonList[Q]->radRecombCross(EeV, 0);

    std::vector<double> PI1 = ch.IonList[Q]->photoIonizationCross(EeV, 0);

    std::vector< std::tuple<int, double> >
      REv = ch.IonList[Q]->recombCrossV(EeV, 0);

    std::sort(begin(REv), end(REv), TupleCompare<0>());
    double sum = 0.0;
    std::vector<double> cum;
    for (auto t : REv) {
      double val = std::get<1>(t);
      sum += val;
      cum.push_back(val);
    }
    for (auto & t : cum) t /= sum;

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
	      << std::setw(16) << geom       * 1.0e+04 // Mb
	      << std::setw(16) << elas       * 1.0e+04 // Mb
	      << std::setw(16) << coul       * 1.0e+04 // Mb
	      << std::setw(16) << ffre       * 1.0e+04 // Mb
	      << std::setw(16) << coll       * 1.0e+04 // Mb
	      << std::setw(16) << ionz       * 1.0e+04 // Mb
	      << std::setw(16) << RE1.back() * 1.0e+04 // Mb
	      << std::setw(16) << PI1.back() * 1.0e+04 // Mb
	      << std::setw(16) << csum;		       // Mb
    for (auto t : cum) std::cout << std::setw(16) << t;
    std::cout << std::endl;
  }

  MPI_Finalize();

  return 0;
}
