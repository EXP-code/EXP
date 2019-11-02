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
  double Tmin, Tmax, Emin, Emax;
  std::string RRtype;
  int numT, norder, numE;
  size_t nout = 1;
  bool rates = false, use_log = false;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("long,l",		"long-form output")
    ("Z,Z",		po::value<unsigned short>(&Z)->default_value(2),
     "atomic number")
    ("Tmin,t",		po::value<double>(&Tmin)->default_value(1000.0),
     "minimum temperature")
    ("Tmax,T",		po::value<double>(&Tmax)->default_value(1000000.0),
     "maximum temperature")
    ("Emin,e",		po::value<double>(&Emin)->default_value(0.01),
     "minimum energy in eV")
    ("Emax,E",		po::value<double>(&Emax)->default_value(100.0),
     "maximum energy in eV")
    ("NumT,N",		po::value<int>(&numT)->default_value(200),
     "number of temperature points")
    ("NumE,M",		po::value<int>(&numE)->default_value(1),
     "number of incremental energy points")
    ("rates",
     "plot VALUES of rate coefficients")
    ("use_log",
     "use logarithmic spacing in Legendre integration")
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

  if (vm.count("long")) nout = 3;

  if (vm.count("rates")) rates = true;

  if (vm.count("log")) use_log = true;

  std::string prefix("IonRecombFrac");
  std::string cmdFile = prefix + ".cmd_line";
  std::ofstream cmd(cmdFile.c_str());
  if (!cmd) {
    std::cerr << "testCrossSection: error opening <" << cmdFile
	      << "> for writing" << std::endl;
  } else {
    cmd << cmd_line << std::endl;
  }

  cmd.close();

  std::vector<std::ofstream> out(nout);
  out[0].open(prefix + ".data", ios::out);
  if (out[0].bad()) {
    std::cerr << "testIonRecomb: error opening <" << prefix + ".data"
	      << "> for writing" << std::endl;
    exit(-1);
  }

  if (nout==3) {
    out[1].open(prefix + ".ionize", ios::out);
    if (out[1].bad()) {
      std::cerr << "testIonRecomb: error opening <" << prefix + ".ionize"
		<< "> for writing" << std::endl;
      exit(-1);
    }

    out[2].open(prefix + ".recomb", ios::out);
    if (out[2].bad()) {
      std::cerr << "testIonRecomb: error opening <" << prefix + ".recomb"
		<< "> for writing" << std::endl;
      exit(-1);
    }
  }

  std::ofstream mat, tab;
  if (numE>1) {
    mat.open(prefix + ".matrix", ios::out);
    if (mat.bad()) {
      std::cerr << "testIonRecomb: error opening <" << prefix + ".matrix"
		<< "> for writing" << std::endl;
      exit(-1);
    }

    if (rates) {

      tab.open(prefix + ".table", ios::out);
      if (tab.bad()) {
	std::cerr << "testIonRecomb: error opening <" << prefix + ".table"
		  << "> for writing" << std::endl;
	exit(-1);
      }
    }
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

  bool logE = false;
  if (Emin>0.0) {
    Emin = log(Emin);
    Emax = log(Emax);
    logE = true;
  }
  const size_t NE = 100;
  double dE = (Emax - Emin)/NE;

  for (int nt=0; nt<=numT; nt++) {

    double T = exp(tmin + dT*nt);

    std::map<unsigned short, std::vector<double> > val1 = ch.recombEquil(Z, T, norder);

    double emin0 = Emin;
    double emax0 = Emax;
    if (logE) {
      emin0 = exp(emin0);
      emax0 = exp(emax0);
    }
    std::map<unsigned short, std::vector<double> > val2 = ch.recombEquil(Z, T, emin0, emax0, norder, use_log);

    std::map<unsigned short, std::vector<double> > val3;
    for (size_t ne=0; ne<NE; ne++) {
      double Eb = Emin + dE*ne;
      double Ef = Emin + dE*(ne+1);
      if (logE) {
	Eb = exp(Eb);
	Ef = exp(Ef);
      }
      if (val3.size()) {
	std::map<unsigned short, std::vector<double> > valT =
	  ch.recombEquil(Z, T, Eb, Ef, norder, use_log);
	for (auto v : val3) {
	  unsigned short C = v.first;
	  size_t sz = v.second.size();
	  for (size_t k=0; k<sz; k++) val3[C][k] += valT[C][k];
	}
      } else {
	val3 = ch.recombEquil(Z, T, Eb, Ef, norder, use_log);
      }
    }

    // Renormalize
    std::vector<double> nn(Z+1, 1);
    double norm = 1.0;
    for (unsigned short C=1; C<Z+1; C++) {
      if (val3[C+1][2]>0.0)
	nn[C] = nn[C-1] * val3[C][1]/val3[C+1][2];
      else
	nn[C] = 0.0;
      norm += nn[C];
    }
    for (unsigned short C=1; C<Z+2; C++) val3[C][0] = nn[C-1]/norm;
    
    for (size_t n=0; n<nout; n++) {
      out[n] << std::setw(16) << T;
      for (auto v : val1) out[n] << std::setw(16) << v.second[n];
      for (auto v : val2) out[n] << std::setw(16) << v.second[n];
      for (auto v : val3) out[n] << std::setw(16) << v.second[n];
      out[n] << std::endl;
    }

    if (numE>1) {

      typedef std::map<unsigned short, std::vector<double> > rateMap;

      rateMap valH, val0;
      std::vector<rateMap> val1(numE);
      
      if (rates) valH = ch.recombEquil(Z, T, norder);

      for (int ne=0; ne<numE; ne++) {
	
	double Eb = Emin + dE*ne;
	double Ef = Emin + dE*(ne+1);
	if (logE) {
	  Eb = exp(Eb);
	  Ef = exp(Ef);
	}

	if (val0.size()) {
	  std::map<unsigned short, std::vector<double> >
	    valT = ch.recombEquil(Z, T, Eb, Ef, norder);
	  for (auto v : val0) {
	    unsigned short C = v.first;
	    size_t        sz = v.second.size();
	    for (size_t j=0; j<sz; j++) val0[C][j] += valT[C][j];
	  }
	} else {
	  val0 = ch.recombEquil(Z, T, Ef, Eb, norder);
	}
	val1[ne] = val0;
      }
      
      for (int ne=0; ne<numE; ne++) {
	
	double Ef = Emin + dE*(ne+1);
	if (logE) Ef = exp(Ef);

	mat << std::setw(16) << T 
	    << std::setw(16) << Ef;

	for (auto v : val1[ne]) {
	  unsigned short C = v.first;
	  for (size_t j=1; j<3; j++) {
	    if (rates) {
	      mat << std::setw( 3) << C
		  << std::setw(16) << valH[C][j]
		  << std::setw(16) << val0[C][j]
		  << std::setw(16) << val1[ne][C][j];
	    } else if (val0[C][j]>0.0) {
	      mat << std::setw(16) << val1[ne][C][j]/val0[C][j];
	    } else {
	      mat << std::setw(16) << 0.0;
	    }
	  }
	}
	mat << std::endl;
      }
      mat << std::endl;

      if (rates) {

	tab << std::setw(16) << T ;
	for (auto v : val1.back()) {
	  unsigned short C = v.first;
	  for (size_t j=1; j<3; j++) {
	    tab << std::setw( 3) << C
		<< std::setw(16) << valH[C][j]
		<< std::setw(16) << val0[C][j]
		<< std::setw(16) << v.second[j];
	  }
	}
	tab << std::endl;
      }
    }
  }


  MPI_Finalize();

  return 0;
}
