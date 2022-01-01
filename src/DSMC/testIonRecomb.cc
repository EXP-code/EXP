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

// Global variables (not all used but needed for EXP libraries)

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

  cxxopts::Options options(av[0], "Test recombination cross section");

  options.add_options()
   ("h,help", "produce help message")
   ("l,long", "long-form output")
   ("Z,Z", "atomic number",
     cxxopts::value<unsigned short>(Z)->default_value("2"))
   ("t,Tmin", "minimum temperature",
     cxxopts::value<double>(Tmin)->default_value("1000.0"))
   ("T,Tmax", "maximum temperature",
     cxxopts::value<double>(Tmax)->default_value("1000000.0"))
   ("e,Emin", "minimum energy in eV",
     cxxopts::value<double>(Emin)->default_value("0.001"))
   ("E,Emax", "maximum energy in eV",
     cxxopts::value<double>(Emax)->default_value("1000.0"))
   ("N,NumT", "number of temperature points",
     cxxopts::value<int>(numT)->default_value("200"))
   ("M,NumE", "number of incremental energy points",
     cxxopts::value<int>(numE)->default_value("1"))
   ("rates", "plot VALUES of rate coefficients")
   ("use_log", "use logarithmic spacing in Legendre integration")
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

  bool first_table = true;


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

  atomicData ad;

  ad.createIonList(ZList);

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

    std::map<unsigned short, std::vector<double> > val1 = ad.recombEquil(Z, T, norder);

    double emin0 = Emin;
    double emax0 = Emax;
    if (logE) {
      emin0 = exp(emin0);
      emax0 = exp(emax0);
    }
    std::map<unsigned short, std::vector<double> > val2 = ad.recombEquil(Z, T, emin0, emax0, norder, use_log);

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
	  ad.recombEquil(Z, T, Eb, Ef, norder, use_log);
	for (auto v : val3) {
	  unsigned short C = v.first;
	  size_t sz = v.second.size();
	  for (size_t k=0; k<sz; k++) val3[C][k] += valT[C][k];
	}
      } else {
	val3 = ad.recombEquil(Z, T, Eb, Ef, norder, use_log);
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
      
      if (rates) valH = ad.recombEquil(Z, T, norder);

      for (int ne=0; ne<numE; ne++) {
	
	double Eb = Emin + dE*ne;
	double Ef = Emin + dE*(ne+1);
	if (logE) {
	  Eb = exp(Eb);
	  Ef = exp(Ef);
	}

	if (val0.size()) {
	  std::map<unsigned short, std::vector<double> >
	    valT = ad.recombEquil(Z, T, Eb, Ef, norder);
	  for (auto v : val0) {
	    unsigned short C = v.first;
	    size_t        sz = v.second.size();
	    for (size_t j=0; j<sz; j++) val0[C][j] += valT[C][j];
	  }
	} else {
	  val0 = ad.recombEquil(Z, T, Ef, Eb, norder);
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

	if (first_table) {
	  tab << "#" << std::endl
	      << "#" << std::setw(15) << "temp";
	  for (auto v : val1.back()) {
	    tab << std::setw(6) << "C";
	    std::ostringstream sout;
	    for (int i=1; i<=3; i++) {
	      sout.str(""); sout << "ionize[" << i << "]";
	      tab << std::setw(16) << sout.str();
	    }
	    for (int i=1; i<=3; i++) {
	      sout.str(""); sout << "recomb[" << i << "]";
	      tab << std::setw(16) << sout.str();
	    }
	  }
	  tab << std::endl;

	  int cnt = 1;
	  tab << "#" << std::setw(15) << cnt++;
	  for (auto v : val1.back()) {
	    std::ostringstream sout;
	    sout.str(""); sout << "[" << cnt++ << "]";
	    tab << std::setw(6) << sout.str();
	    for (size_t j=0; j<6; j++) {
	      sout.str(""); sout << "[" << cnt++ << "]";
	      tab << std::setw(16) << sout.str();
	    }
	  }
	  tab << std::endl;

	  tab << "#" << std::setw(15) << "----------";
	  for (auto v : val1.back()) {
	    tab << std::setw(6) << "----";
	    for (size_t j=0; j<6; j++)
	      tab << std::setw(16) << "----------";
	  }
	  tab << std::endl;

	  first_table = false;
	}

	tab << std::setw(16) << T;
	for (auto v : val1.back()) {
	  unsigned short C = v.first;
	  tab << std::setw(6) << C;
	  for (size_t j=1; j<3; j++) {
	    tab << std::setw(16) << valH[C][j]
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
