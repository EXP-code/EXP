/*
  Make recombination ratio caches
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <random>
#include <tuple>

#include <yaml-cpp/yaml.h>

#include <errno.h>
#include <sys/stat.h>

#include <cxxopts.H>
#include "Ion.H"

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
std::mt19937 random_gen;

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

  const double Emin0 = 1.0e-3;	// Minimum energy in eV
  const double Emax0 = 4000.0;	// Maximum energy in eV
  const int    numE0 = 200;	// Intervals in log energy
  const int   norder = 100;	// Number of Legendre knots
  
  unsigned short Z;
  double Tmin, Tmax;
  std::string RRtype, prefix;
  int numT;

  cxxopts::Options options(av[0], "Make recombination ratio caches");

  options.add_options()
   ("h,help", "produce help message")
   ("l,long", "long output: print each rate and ratio")
   ("nogrid", "turn off recombination cross section cache")
   ("Z,Z", "atomic number",
     cxxopts::value<unsigned short>(Z)->default_value("2"))
   ("t,Tmin", "minimum temperature",
     cxxopts::value<double>(Tmin)->default_value("1000.0"))
   ("T,Tmax", "maximum temperature",
     cxxopts::value<double>(Tmax)->default_value("10000000.0"))
   ("N,NumT", "number of temperature points",
     cxxopts::value<int>(numT)->default_value("400"))
   ("R,RRtype", "cross-section type",
     cxxopts::value<std::string>(RRtype)->default_value("Verner"))
   ("p,prefix", "cache name prefix",
     cxxopts::value<std::string>(prefix)->default_value(".recomb_ratio_test"))
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
    std::cout << "\t" << av[0]
	      << " -Z 1" << std::endl;
    MPI_Finalize();
    return 1;
  }

  std::ostringstream fout;
  fout << prefix << "_" << Z;

  if (file_exists(fout.str())) {
    std::cout << av[0] << ": output file <" << fout.str() << "> exists"
	      << std::endl;
    exit(-1);
  }

  std::ofstream out(fout.str());

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

  // Close and remove temporary file
  //
  in.close();

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

  if (vm.count("nogrid")) Ion::useRadRecombGrid = false;
  Ion::setRRtype(RRtype);

  atomicData ad;

  ad.createIonList(ZList);

  std::cout << "RRtype: " << Ion::getRRtype() << std::endl;

  // Log ranges for temperature/ratio grid
  //
  double Tmn = log(Tmin);
  double Tmx = log(Tmax);
  double dT  = (Tmx - Tmn)/(numT-1);

  // Log ranges for energy integration with a minimum interval
  //
  double Emin = log(Emin0);
  double Emax = log(Emax0);
  int    numE = std::max<int>(numE0, 10);
    
  typedef std::map<unsigned short, std::vector<double> > rateMap;

  double dE = (Emax - Emin)/numE;

  for (int nt=0; nt<numT; nt++) {

    double T = temp[nt];
    
    rateMap val0;
    std::vector<rateMap> val1(numE);
    
    for (int ne=0; ne<numE; ne++) {
      
      double Eb = Emin + dE*ne;
      double Ef = Emin + dE*(ne+1);

      Eb = exp(Eb);
      Ef = exp(Ef);
      
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
      
    for (auto v : val1.back()) {
      unsigned short C = v.first;
      if (C>1) {
	if (v.second[2]>0.0)
	  cdata[C-2][nt] = cdata[C-2][nt]/(v.second[2]*1.0e-14);
	else
	  cdata[C-2][nt] = 0.0;
      }
    }
  }

  // Write cache file
  //
  YAML::Node node;
    
  node["Z"   ] = Z;
  node["Tmin"] = Tmin;
  node["Tmax"] = Tmax;
  node["numT"] = numT;
    
  // Serialize the node
  //
  YAML::Emitter y; y << node;
  
  // Get the size of the string
  //
  unsigned int hsize = strlen(y.c_str());
  
  // Write YAML string size
  //
  out.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));
  
  // Write YAML string
  //
  out.write(reinterpret_cast<const char *>(y.c_str()), hsize);
  
  // Write data base of values
  //
  for (int nt=0; nt<numT; nt++) {
    for (int C=2; C<=Z+1; C++)
      out.write(reinterpret_cast<const char *>(&cdata[C-2][nt]), sizeof(double));
  }
  //
  // Cache written!

  return 0;
}
