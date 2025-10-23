/*
  Separate a psp structure and make a 2d density grid

  MDWeinberg 07/11/20
*/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>

using namespace std;

#include "Progress.H"		// The old boost progress bar
#include "cxxopts.H"		// Option parsing
#include "libvars.H"		// EXP library globals

#include "StringTok.H"
#include "FileUtils.H"
#include "header.H"
#include "Grid2D.H"
#include "PSP.H"

//
// MPI variables
//
int numprocs, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

int
main(int ac, char **av)
{
  // MPI preliminaries
  //
  MPI_Init(&ac, &av);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

  char *prog = av[0];
  bool verbose = false;
  double rmax;
  int numr, ibeg, iend;
  std::string cname, dir, outp, fpre, runtag;

  // Parse command line
  //
  cxxopts::Options options(prog, "Separate a psp structure and make a 2d density grid\n");

  options.add_options()
   ("h,help", "produce help message")
   ("v,verbose", "verbose output")
   ("R,rmax", "maximum component extent",
     cxxopts::value<double>(rmax)->default_value("0.05"))
   ("N,numr", "number of grid points per dimension",
     cxxopts::value<int>(numr)->default_value("40"))
   ("c,comp", "component name",
     cxxopts::value<std::string>(cname)->default_value("star"))
   ("d,dir", "output directory",
     cxxopts::value<std::string>(dir)->default_value("./"))
   ("o,output", "output prefix",
     cxxopts::value<std::string>(outp)->default_value("densgrid"))
   ("p,psptype", "PSP type (OUT or SPL)",
     cxxopts::value<std::string>(fpre)->default_value("OUT"))
   ("r,runtag", "run tag",
     cxxopts::value<std::string>(runtag)->default_value("run001"))
   ("i,beg", "initial snap shot number",
     cxxopts::value<int>(ibeg)->default_value("0"))
   ("f,end", "final snap shot number",
     cxxopts::value<int>(iend)->default_value(std::to_string(std::numeric_limits<int>::max())))
    ;

  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }


  if (vm.count("help")) {
    if (myid==0) {
      std::cout << options.help() << std::endl;
      std::cout << "Example: " << std::endl;
      std::cout << "\t" << av[0]
		<< " --runtag=run001" << std::endl;
    }
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  std::ofstream out;

  int nok = 0;
  if (myid==0) {
    out.open(outp + "." + runtag);
    if (not out) {
      std::cerr << av[0] << ": could not open <"
		<< outp << "." << runtag << ">"
		<< std::endl;
      nok = 1;
    }
  }

  // File open error
  //
  MPI_Bcast(&nok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nok) {
    MPI_Finalize();
    exit(-1);
  }

  // Input file range check
  //
  int n;
  for (n=ibeg; n<=iend; n++) {

    std::ostringstream fname;
    fname << fpre << "." << runtag << "."
	  << std::setw(5) << std::setfill('0') << n;

    std::string file = fname.str();

    if (!FileExists(file)) {
      if (myid==0) {
	std::cerr << "Error opening file <" << file << "> for input"
		  << std::endl;
      }
      break;
    }
  }
  iend = n-1;
  if (myid==0) {
    std::cerr << "Assuming last file has index <" << iend << ">"
	      << std::endl;
  }

  std::shared_ptr<progress::progress_display> progress;
  if (myid==0) {
    progress = std::make_shared<progress::progress_display>(iend - ibeg + 1);
  }

  // Snapshot loop
  //
  for (int n=ibeg; n<=iend; n++) {
    
    std::ostringstream file;
    file << fpre << "." << runtag << "."
	 << std::setw(5) << std::setfill('0') << n;

    // Parse the PSP file
    //
    PSPptr psp;
    if (file.str().find("SPL") != std::string::npos)
      psp = std::make_shared<PSPspl>(file.str(), dir);
    else
      psp = std::make_shared<PSPout>(file.str());


    // Now write a summary
    //
    if (verbose) {

      psp->PrintSummary(cerr);
    
      std::cerr << "\nPSP file named <" << file.str() << "> has time <" 
		<< psp->CurrentTime() << ">\n";
    }


    Grid2D grid(rmax, numr, psp->CurrentTime());

    PSPstanza* stanza;
    SParticle* part;

    unsigned count = 0;

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;
      
      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {
	if (count++ % numprocs != myid) continue;
	grid.addPoint(part->mass(), part->pos(0), part->pos(1));
      }
    }

    grid.sync();

    if (myid==0) {
      grid.write(out);
      ++(*progress);
    }

  }
    
  MPI_Finalize();

  return 0;
}
  
