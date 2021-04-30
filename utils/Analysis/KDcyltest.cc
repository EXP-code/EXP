/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  k-d density test for cylinder
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/27/20
 *
 ***************************************************************************/

				// C++/STL headers
#include <numeric>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <queue>
#include <map>

				// Boost stuff

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/make_unique.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <Progress.H>

namespace po = boost::program_options;

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <Vector.h>
#include <numerical.h>
#include "Particle.h"
#include <PSP2.H>
#include <foarray.H>
#include <KDtree.H>

// Variables not used but needed for linking
//
int VERBOSE = 4;
int nthrds = 1;
int myid = 0;
int this_step = 0;
unsigned multistep = 0;
unsigned maxlev = 100;
int mstep = 1;
int Mstep = 1;
vector<int> stepL(1, 0), stepN(1, 1);
char threading_on = 0;
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
std::string outdir, runtag;
double tpos = 0.0;
double tnow = 0.0;
  
int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double ROUT, ZOUT;
  int NICE, indx, Ndens, NOUT;
  std::string dir("./"), cname, prefix;
  bool ignore;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "K-D density estimate test" << std::endl
       << std::string(60, '-') << std::endl << std::endl
       << "Allowed options";
  
  po::options_description desc(sout.str());
  desc.add_options()
    ("help,h",
     "Print this help message")
    ("OUT",
     "assume original, single binary PSP files as input")
    ("SPL",
     "assume new split binary PSP files as input")
    ("Ndens,K",             po::value<int>(&Ndens)->default_value(32),
     "KD density estimate count (use 0 for expansion estimate)")
    ("NICE",                po::value<int>(&NICE)->default_value(0),
     "system priority")
    ("NOUT",                po::value<int>(&NOUT)->default_value(40),
     "Number of grid points for surface output")
    ("ROUT",                po::value<double>(&ROUT)->default_value(0.1),
     "Maximum radius for output")
    ("ZOUT",                po::value<double>(&ZOUT)->default_value(0.02),
     "Maximum height for output")
    ("prefix",              po::value<string>(&prefix)->default_value("crossval"),
     "Filename prefix")
    ("runtag",              po::value<string>(&runtag)->default_value("run1"),
     "Phase space file")
    ("outdir",              po::value<string>(&outdir)->default_value("."),
     "Output directory path")
    ("indx",                po::value<int>(&indx)->default_value(0),
     "PSP index")
    ("dir,d",               po::value<std::string>(&dir),
     "directory for SPL files")
    ("cname",
     po::value<std::string>(&cname)->default_value("star disk"),
     "component name")
    ;
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // ==================================================
  // Print help message and exit
  // ==================================================

  if (vm.count("help")) {
    if (myid==0) std::cout << std::endl << desc << std::endl;
    MPI_Finalize();
    return 0;
  }

  bool SPL = false;
  if (vm.count("SPL")) SPL = true;
  if (vm.count("OUT")) SPL = false;

  bool LOG = false;
  if (vm.count("LOG")) {
    LOG = true;
  }

  bool Hall = false;
  if (vm.count("Hall")) Hall = true;

  bool verbose = false;
  if (vm.count("verbose")) verbose = true;

  bool debug = false;
  if (vm.count("debug")) debug = true;

  // ==================================================
  // Nice process
  // ==================================================

  if (NICE>0)
    setpriority(PRIO_PROCESS, 0, NICE);

  // ==================================================
  // Phase space
  // ==================================================

  std::string file;

#ifdef DEBUG
  std::cout << "[" << myid << "] Begin phase -space loop" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif	      

  // ==================================================
  // PSP input stream
  // ==================================================

  int iok = 1;
  std::ostringstream s1;
  
  s1.str("");		// Clear stringstream
  if (SPL) s1 << "SPL.";
  else     s1 << "OUT.";
  s1 << runtag << "."<< std::setw(5) << std::setfill('0') << indx;
      
				// Check for existence of next file
  file = dir + s1.str();
  std::ifstream in(file);
  if (!in) {
    std::cerr << "Error opening <" << file << ">" << endl;
    iok = 0;
  }
  
  if (iok==0) {
    MPI_Finalize();
    exit(-1);
  }

  // ==================================================
  // Open output file
  // ==================================================

  std::ofstream out;
  bool ok = true;
  if (myid==0) {
    out.open(prefix + ".out");
    if (!out) {
      std::cerr << "Error opening output file <" << prefix + ".out" << ">" << std::endl;
      ok = false;
    }
  }
  
  {
    int okay = ok ? 1 : 0;
    MPI_Bcast(&okay, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (okay==0) {
      MPI_Finalize();
      exit(-2);
    }
  }

  // ==================================================
  // Open PSP file
  // ==================================================

  PSPptr psp;

  if (SPL) psp = std::make_shared<PSPspl>(s1.str(), dir, true);
  else     psp = std::make_shared<PSPout>(file, true);
  
  tnow = psp->CurrentTime();
  if (myid==0) std::cout << "Beginning partition [time=" << tnow
			 << ", index=" << indx << "] . . . "  << flush;
  
  if (not psp->GetNamed(cname)) {
    if (myid==0) {
      std::cout << "Error finding component named <" << cname << ">" << std::endl;
      psp->PrintSummary(std::cout);
    }
    MPI_Finalize();
    exit(-1);
  }
      
  int nbod = psp->GetNamed(cname)->comp.nbod;

  std::vector<double> KDdens;

  boost::shared_ptr<boost::progress_display> progress;
  if (myid==0) {
    std::cout << std::endl
	      << "Accumulating particle positions . . . "
	      << std::endl;
    progress = boost::make_shared<boost::progress_display>(nbod);
  }

  SParticle *p = psp->GetParticle();
  int icnt = 0;

  // This is the kd- NN density estimate; skipped by default for Ndens=0
  //
  if (Ndens<=0) Ndens = 2;
  if (myid==0) std::cout << "Computing KD density estimate for " << nbod
			 << " points" << std::endl;

  typedef point <double, 3> point3;
  typedef kdtree<double, 3> tree3;

  std::vector<point3> points;

  // Every node needs to make the tree (a parallel share could be
  // implemented in KDtree.H)
  //
  double KDmass = 0.0;
  for (auto part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {
    double ms = part->mass();
    KDmass += ms;
    points.push_back(point3({part->pos(0), part->pos(1), part->pos(2)}, ms));
  }
    
  tree3 tree(points.begin(), points.end());
    
  KDdens.resize(nbod, 0.0);

  int badVol = 0;

  double dR = 2.0*ROUT/(NOUT-1);
  double dZ = 2.0*ZOUT/(NOUT-1);

  for (int j=0; j<NOUT; j++) {
    double Z = -ZOUT + dZ*j;

    for (int i=0; i<NOUT; i++) {
      double R = -ROUT + dR*i;

      auto ret = tree.nearestN({R, 0.0, Z}, Ndens);
      double volume = 4.0*M_PI/3.0*std::pow(std::get<2>(ret), 3.0);
      double densty = 0.0;
      if (volume>0.0) densty = std::get<1>(ret)/volume;

      out << std::setw(12) << R
	  << std::setw(12) << Z
	  << std::setw(12) << densty
	  << std::endl;
    }
    out << std::endl;
  }

  return 0;
}

