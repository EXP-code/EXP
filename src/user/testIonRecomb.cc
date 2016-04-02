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

#include <gaussQ.h>

#include "atomic_constants.H"
#include "Ion.H"
#include "Elastic.H"

namespace po = boost::program_options;

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

double df(double v, double T)
{
  const double b = 0.5*me/boltz;
  const double c = 4.0/sqrt(M_PI)*pow(b, 1.5);
  double v2 = v*v;
  return c * pow(T, -1.5) * v2 * exp(-b*v2);
}

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
  bool rates;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("rates,r",         "print rate coefficients")
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

  if (vm.count("rates")) rates = true;
  else                   rates = false;

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    std::cout << "Example: Helium II recombination" << std::endl;
    std::cout << "\t" << av[0]
	      << " -Z 2 -C 2" << std::endl;
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

  std::map<lQ, double> ionize, recomb;

  for (unsigned short C=1; C<=Z+0; C++) ionize[lQ(Z, C)] = 0.0;
  for (unsigned short C=2; C<=Z+1; C++) recomb[lQ(Z, C)] = 0.0;

  LaguQuad lq(norder);

  const double beta = 0.5*me/boltz;

  for (int nt=0; nt<=numT; nt++) {

    double T = exp(tmin + dT*nt);

    for (auto & v : ionize) v.second = 0.0;
    for (auto & v : recomb) v.second = 0.0;

    for (int i=1; i<=norder; i++) {
      double y = lq.knot(i);
      double w = lq.weight(i);

      // Ionize
      for (auto & v : ionize) {
	double ionE = ch.ipdata[v.first];  // ionization energy in eV
	double ab   = ionE/(boltzEv*T);	  
	double eab  = exp(-ab);	           // Boltzmann factor for ionization
	double EeV  = (y + ab)*boltzEv*T;  // eval energy
	double vel  = sqrt((y + ab)/beta); // eval velocity squared

	v.second += w * eab * vel *
	  ch.IonList[v.first]->directIonCross(EeV, 0);
      }

      // Recomb
      for (auto & v : recomb) {
	double EeV  = y*boltzEv*T;
	double vel  = sqrt(y/beta);

	v.second += w * vel *
	  ch.IonList[v.first]->radRecombCross(EeV, 0).back();
      }
    }

    out << std::setw(16) << T;
    if (rates) {
      for (auto & v : ionize) out << std::setw(16) << v.second;
      for (auto & v : recomb) out << std::setw(16) << v.second;
    }

    std::vector<double> nn(Z+1, 1);
    double norm = 1.0;
    for (unsigned short C=1; C<=Z; C++) {
      nn[C] = nn[C-1] * ionize[lQ(Z,C)]/recomb[lQ(Z,C+1)];
      norm += nn[C];
    }
    for (auto n : nn) out << std::setw(16) << n/norm;
    out << std::endl;
  }

  MPI_Finalize();

  return 0;
}
