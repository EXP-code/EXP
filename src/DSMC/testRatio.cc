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
#include <boost/filesystem.hpp>

#include <errno.h>
#include <sys/stat.h>

#include "Ion.H"

namespace po = boost::program_options;

#include <mpi.h>

bool file_exists(const std::string& fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}


// Write ChiantiPy script, if necessary
//
void writeScript()
{
  const char *py =
#include "recomb_py.h"
    ;

  const std::string file("recomb.py");
    
  if (not file_exists(file)) {
    std::ofstream out(file);
    out << py;
    if (chmod(file.c_str(), S_IWUSR | S_IRUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
      perror("Error in chmod:");
    }
  }
}


// Global variables (not all used but needed for EXP libraries)

int numprocs, myid;
std::string outdir(".");
std::string runtag("run");
char threading_on = 0;
pthread_mutex_t mem_lock;

int main (int ac, char **av)
{
  writeScript();

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

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("long,l",          "long output: print each rate and ratio")
    ("Z,Z",		po::value<unsigned short>(&Z)->default_value(2),
     "atomic number")
    ("Tmin,t",		po::value<double>(&Tmin)->default_value(1000.0),
     "minimum temperature")
    ("Tmax,T",		po::value<double>(&Tmax)->default_value(10000000.0),
     "maximum temperature")
    ("Emin,e",		po::value<double>(&Emin)->default_value(0.001),
     "minimum energy in eV")
    ("Emax,E",		po::value<double>(&Emax)->default_value(1000.0),
     "maximum energy in eV")
    ("NumT,N",		po::value<int>(&numT)->default_value(200),
     "number of temperature points")
    ("NumE,M",		po::value<int>(&numE)->default_value(100),
     "number of incremental energy points")
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

  std::string prefix("testRatio");
  std::string cmdFile = prefix + ".cmd_line";
  std::ofstream cmd(cmdFile.c_str());
  if (!cmd) {
    std::cerr << "testCrossSection: error opening <" << cmdFile
	      << "> for writing" << std::endl;
  } else {
    cmd << cmd_line << std::endl;
  }

  cmd.close();

  // Get recombationation rates from ChiantiPy
  //
  std::string inFile("testRatio.in");
  std::ostringstream sout;
  sout << "$PYTHON ./recomb.py"
       << " -Z " << Z
       << " -t " << Tmin
       << " -T " << Tmax
       << " -n " << numT
       << " -o " << inFile << '\0';

  std::cout << "Command: " << &sout.str()[0] << std::endl;
  int sret = system(&sout.str()[0]);

  std::ifstream in(inFile);

  // For temperature data
  //
  std::vector<double> temp(numT);

  // For recombination rate data
  //
  std::vector<std::vector<double>> cdata(Z);
  for (int C=2; C<=Z+1; C++) cdata[C-2].resize(numT);

  // Get data from generated file
  //
  for (int nt=0; nt<numT; nt++) {
    std::string line;
    if (std::getline(in, line)) {
      std::istringstream ins(line);
      ins >> temp[nt];
      for (int C=2; C<=Z+1; C++) ins >> cdata[C-2][nt];
    }
  }

  // Open output file
  //
  std::ofstream rcb(prefix + ".ratio", ios::out);
  if (rcb.bad()) {
    std::cerr << "testRatio: error opening <" << prefix + ".ratio"
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

  bool logE = false;
  if (Emin>0.0) {
    Emin = log(Emin);
    Emax = log(Emax);
    logE = true;
  }

  typedef std::map<unsigned short, std::vector<double> > rateMap;

  numE = std::max<int>(numE, 10);
  
  double dE = (Emax - Emin)/numE;

  for (int nt=0; nt<numT; nt++) {

    double T = temp[nt];

    rateMap val0;
    std::vector<rateMap> val1(numE);
    
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
      
    rcb << std::setw(16) << T;
    for (auto v : val1.back()) {
      unsigned short C = v.first;
      if (C>1) {
	if (cdata[C-2][nt]>0.0) {
	  if (vm.count("long"))
	    rcb << std::setw(16) << v.second[2]*1.0e-14
		<< std::setw(16) << cdata[C-2][nt];
	  rcb 	<< std::setw(16) << v.second[2]*1.0e-14/cdata[C-2][nt];
	} else {
	  if (vm.count("long"))
	    rcb << std::setw(16) << v.second[2]*1.0e-14
		<< std::setw(16) << cdata[C-2][nt];
	  rcb 	<< std::setw(16) << 0.0;
	}
      }
    }
    rcb << std::endl;
  }

  MPI_Finalize();

  return 0;
}
