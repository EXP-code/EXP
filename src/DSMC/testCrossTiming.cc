/*
  Test cross section timing
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <chrono>
#include <random>

#include "atomic_constants.H"
#include "Ion.H"
#include "Elastic.H"
#include <cxxopts.H>

#include <mpi.h>

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

  cxxopts::Options options(av[0], "Test runtime performance of cross section evaluation");

  options.add_options()
    ("h,help", "produce help message")
    ("eV", "print results in eV")
    ("Z,Z", "atomic number",
     cxxopts::value<unsigned short>(&Z)->default_value("2"))
    ("C,C", "ionic charge",
     cxxopts::value<unsigned short>(&C)->default_value("3"))
    ("e,Emin", "minimum energy (Rydbergs)",
     cxxopts::value<double>(&emin)->default_value("0.001"))
    ("E,Emax", "maximum energy (Rydbergs)",
     cxxopts::value<double>(&emax)->default_value("100.0"))
    ("N,Num", "number of evaluations",
     cxxopts::value<int>(&num)->default_value("1000000"))
    ("k,kdel", "default logarithmic spacing for k grid",
     cxxopts::value<double>(&kdel)->default_value("0.25"))
    ("R,RRtype", "cross-section type",
     cxxopts::value<std::string>(&RRtype)->default_value("Verner"))
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
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

  std::string prefix("crossSectionTiming");
  std::string cmdFile = prefix + ".cmd_line";
  std::ofstream out(cmdFile.c_str());
  if (!out) {
    std::cerr << "crossSectionTiming: error opening <" << cmdFile
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
				// Convert from Rydberg to eV
  const double ryd      = 27.2113845/2.0;
  const double eVerg    = 1.60217733e-12;

  emin = log(emin);
  emax = log(emax);

  double dE = emax - emin;

  lQ Q(Z, C), QL(Z, C-1);

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> u(0.0, 1.0);
  
  std::map<std::string, double> data;

  for (int i=0; i<num; i++) {
    double E   = exp(emin + dE*u(mt));
    double EeV = E * ryd;

    double geom = geometric(Z);
    double elas = elastic(Z, EeV);
    double coul = 0.0;

    {
      auto timer = std::chrono::high_resolution_clock::now();
      const double ips = 1000.0;
      double b90 = 0.5*esu*esu*(C-1) / (EeV * eVerg) * 1.0e7; // nm
      b90 = std::min<double>(b90, ips);
      coul = M_PI*b90*b90 * 4.0*mrat/pt[Z]->weight() * logL;

      std::chrono::duration<double>  delT =
      (std::chrono::high_resolution_clock::now() - timer);
      data["Coulomb"] += delT.count();
    }

    {
      auto timer = std::chrono::high_resolution_clock::now();
      std::pair<double, double> ffre = ad.IonList[Q]->freeFreeCross(EeV, 0);
      std::chrono::duration<double>  delT =
      (std::chrono::high_resolution_clock::now() - timer);
      data["Free free"] += delT.count();
    }

    {
      auto timer = std::chrono::high_resolution_clock::now();
      double ionz = ad.IonList[Q]->directIonCross(EeV, 0);
      std::chrono::duration<double>  delT =
      (std::chrono::high_resolution_clock::now() - timer);
      data["Coll ionize"] += delT.count();
    }

    {
      auto timer = std::chrono::high_resolution_clock::now();
      Ion::collType       CE1 = ad.IonList[Q]->collExciteCross(EeV, 0);
      std::chrono::duration<double>  delT =
      (std::chrono::high_resolution_clock::now() - timer);
      data["Coll excite"] += delT.count();
    }
    
    {
      auto timer = std::chrono::high_resolution_clock::now();
      std::vector<double> RE1 = ad.IonList[Q]->radRecombCross(EeV, 0);
      std::vector< std::tuple<int, double> >
	REv = ad.IonList[Q]->recombCrossV(EeV, 0);
      std::chrono::duration<double>  delT =
      (std::chrono::high_resolution_clock::now() - timer);
      data["Rad recomb"] += delT.count();
    }

    {
      auto timer = std::chrono::high_resolution_clock::now();
      std::vector<double> PI1 = ad.IonList[Q]->photoIonizationCross(EeV, 0);
      std::chrono::duration<double>  delT =
      (std::chrono::high_resolution_clock::now() - timer);
      data["Ph ionize"] += delT.count();
    }
    
  }

  MPI_Finalize();

  double sum = 0.0;
  for (auto v : data) sum += v.second;

  std::cout << std::setw(20) << "Cross section" << " | "
	    << std::setw(20) << "Time"
	    << " | " << "%%" << std::endl;
  std::cout << std::setw(20) << std::setfill('-') << '-' << "-+-"
	    << std::setw(20) << std::setfill('-') << '-'
	    << "-+-" << "-----" << std::endl << std::setfill(' ');
  for (auto v : data)
    std::cout << std::setw(20) <<  v.first
	      << " | " << std::setw(20) << std::setprecision(6) << v.second
	      << " | " << std::setprecision(3) << 100.0*v.second/sum
	      << std::endl;

  return 0;
}
