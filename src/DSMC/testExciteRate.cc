/*
  Test collisional excitation rate
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <tuple>

#include <cxxopts.H>
#include "atomic_constants.H"
#include "Ion.H"
#include "Elastic.H"

#include <mpi.h>

int numprocs, myid;
std::string outdir(".");
std::string runtag("run");
char threading_on = 0;
pthread_mutex_t mem_lock;
std::mt19937 random_gen;

int main (int ac, char **av)
{
  //===================
  // MPI preliminaries 
  //===================

  MPI_Init(&ac, &av);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  //===================
  // Parse command line
  //===================

  std::string cmd_line;
  for (int i=0; i<ac; i++) {
    cmd_line += av[i];
    cmd_line += " ";
  }

  unsigned short Z, C;
  double tmin, tmax; 
  int num, zmin, zmax, nc;

  cxxopts::Options options(av[0], "Test collisional excitation rate");

  options.add_options()
   ("h,help", "produce help message")
   ("eV", "print results in eV")
   ("zmin", "minimum atomic number",
     cxxopts::value<int>(zmin)->default_value("1"))
   ("zmax", "maximum atomic number",
     cxxopts::value<int>(zmax)->default_value("2"))
   ("Z,Z", "atomic number",
     cxxopts::value<unsigned short>(Z)->default_value("2"))
   ("C,C", "ionic charge",
     cxxopts::value<unsigned short>(C)->default_value("3"))
   ("t,Tmin", "minimum temperature",
     cxxopts::value<double>(tmin)->default_value("10000.0"))
   ("T,Tmax", "maximum temperature",
     cxxopts::value<double>(tmax)->default_value("1000000.0"))
   ("N,Num", "number of temperatures",
     cxxopts::value<int>(num)->default_value("200"))
   ("n,knots", "number of Laguerre knots",
     cxxopts::value<int>(nc)->default_value("40"))
    ;


  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    std::cout << "Example: Helium II recombination" << std::endl;
    std::cout << "\t" << av[0]
	      << " -Z 2 -C 2" << std::endl;
    return 1;
  }

  // Sanity check
  if (Z<zmin or Z>zmax) {
    std::cout << "Your requested value for the base element must in your requested range of [" << zmin << ", " << zmax << "]" << std::endl;
    return -2;
  }

  std::ofstream out("testExciteRate.out");

  std::map<unsigned short, double> abund = {{1, 0.76}, {2, 0.24}};

  std::set<unsigned short> ZList = {1, 2};

  PeriodicTable pt;

  atomicData ad;

  ad.createIonList(ZList);

  double Tmin = log(tmin);
  double Tmax = log(tmax);
  double dT   = (Tmax - Tmin)/(num-1);

  for (int n=0; n<num; n++) {

    double T = exp(Tmin + dT*n);
    
    const std::string ioneq("testExciteRate.ioneq");
    std::ostringstream sout;
    sout << "./genIonization"
         << " -1 " << zmin
	 << " -2 " << zmax
	 << " -T " << T << " -o " << ioneq;
  
    int ret = system(sout.str().c_str());
      
    if (ret) {
      std::cout << "System command  = " << sout.str() << std::endl;
      std::cout << "System ret code = " << ret << std::endl;
    }
    
    typedef std::vector<std::string> vString;
      
    std::vector< std::vector<double> > frac;

    std::string   inLine;
    std::ifstream sFile(ioneq.c_str());
      
    if (sFile.is_open()) {
	
      // Read and discard the headers
      std::getline(sFile, inLine); 

      for (int z=zmin; z<=zmax; z++) {

	std::getline(sFile, inLine);
    
	for (int i=0; i<2; i++) {
	  std::getline(sFile, inLine);

	  // Differential only
	  if (i==0) {
	    vString s;

	    std::istringstream iss(inLine);
	    std::copy(std::istream_iterator<std::string>(iss), 
		      std::istream_iterator<std::string>(), 
		      std::back_inserter<vString>(s));
	    
	    std::vector<double> v;
	    for (auto i : s) v.push_back(::atof(i.c_str()));
	    frac.push_back(v);
	  }
	}
      }
      
      double n_ion = abund[Z]*frac[Z-1][C-1]/pt[Z]->weight();
      double n_elc = 0.0;
      
      double ffE   = 0.0;

      for (int z=zmin; z<=zmax; z++) {
	for (int c=1; c<z+1; c++) {
	  double nd = abund[z] * frac[z-1][c]/pt[z]->weight();
	  n_elc += nd*c;
	  ffE   += nd*ad.freeFreeEmiss(z, c+1, T);
	}
      }
	
      out << std::setw(16) << T
	  << std::setw(16) << n_ion
	  << std::setw(16) << n_elc
	  << std::setw(16) << ad.collEmiss(Z, C, T, 100.0, 40)*n_ion*n_elc
	  << std::setw(16) << ffE*n_elc
	  << std::endl;
    }
  }

  return 0;
}
