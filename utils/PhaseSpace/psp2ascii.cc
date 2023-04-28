/*
  Separate a psp structure to ascii components

  MDWeinberg 06/10/02, 11/24/19
*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <random>
#include <list>

using namespace std;

#include <StringTok.H>
#include <header.H>
#include <PSP.H>
#include <cxxopts.H>

int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time=1e20;
  bool verbose = false;
  bool input   = false;
  std::string cname("comp"), new_dir("./"), filename;

  // Parse command line
  //
  cxxopts::Options options("psp2ascii", "Convert PSP output to ascii for analysis");
  
  options.add_options()
    ("h,help", "this help message")
    ("v,verbose", "verbose output")
    ("a,input", "use input format")
    ("t,time", "desired input time slice",
     cxxopts::value<double>(time)->default_value("1.0e+20"))
    ("o,outname", "prefix name for each component (default: comp)",
     cxxopts::value<std::string>(cname)->default_value("comp"))
    ("d,dirname", "replacement SPL file directory",
     cxxopts::value<std::string>(new_dir)->default_value("./"))
    ("f,filename", "input PSP file",
     cxxopts::value<std::string>(filename))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }


  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid==0) {
      std::cout << options.help() << std::endl;
    }
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  if (vm.count("input")) {
    input = true;
  }


  if (vm.count("filename")) {

    std::ifstream in(filename);
    if (!in) {
      cerr << "Error opening file <" << filename << "> for input\n";
      exit(-1);
    }
    
    if (verbose) cerr << "Using filename: " << filename << endl;

  } else {
    std::cout << options.help() << std::endl;
    exit(-1);
  }
				// Parse the PSP file
				// ------------------
  PSPptr psp;
  if (filename.find("SPL") != std::string::npos)
    psp = std::make_shared<PSPspl>(filename, new_dir);
  else
    psp = std::make_shared<PSPout>(filename);

				// Now write a summary
				// -------------------
  if (verbose) {

    psp->PrintSummary(cerr);
    
    cerr << "\nPSP file <" << filename << "> has time <" 
	 << psp->CurrentTime() << ">\n";
  }

				// Dump ascii for each component
				// -----------------------------
  PSPstanza *stanza;
  SParticle* part;

  for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
				// Open an output file for this stanza
				// -----------------------------------
    ostringstream oname;
    oname << cname << "." << stanza->name << '\0';
    ofstream out(oname.str().c_str());
    out.setf(ios::scientific);
    out.precision(10);

    if (!out) {
      cerr << "Couldn't open output name <" << oname.str() << ">\n";
      exit(-1);
    }
				// Print the header
				// ----------------
    out << setw(15) << stanza->comp.nbod 
	<< setw(10) << stanza->comp.niatr 
	<< setw(10) << stanza->comp.ndatr 
	<< endl;

    const int bunchcnt = 16384;
    int cnt = 0;
    std::ostringstream sout;

    for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

      if (not input and stanza->index_size)
	sout << std::setw(18) << part->indx();

      sout << std::setw(18) << part->mass();
      for (int i=0; i<3; i++) sout << std::setw(18) << part->pos(i);
      for (int i=0; i<3; i++) sout << std::setw(18) << part->vel(i);
      if (not input) sout << std::setw(18) << part->phi();
      for (int i=0; i<part->niatr(); i++)
	sout << std::setw(12) << part->iatr(i);
      for (int i=0; i<part->ndatr(); i++)
	sout << std::setw(18) << part->datr(i);

      sout << std::endl;	// End the record

      if (++cnt==bunchcnt) {	// Write and reset the buffer
	out << sout.str();
	sout.str("");
	cnt = 0;
      }
    }
				// Clear the buffer
    if (sout.str().size()>0) out << sout.str();
    
  }
  
  return 0;
}
  
