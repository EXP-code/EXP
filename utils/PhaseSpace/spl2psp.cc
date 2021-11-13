/*
  Put together a PSP file from per-node components

  MDWeinberg 11/12/19
*/

using namespace std;

#include <unistd.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>

#include <StringTok.H>
#include <header.H>
#include <PSP.H>
#include <Particle.H>

#include <yaml-cpp/yaml.h>
#include <cxxopts.H>

int
main(int argc, char **argv)
{
  char *prog = argv[0];
  std::string runtag, inoutdir, outfile, prefix, new_dir("./");
  bool verbose = false, range = false, checkpt = false;
  int seq, bseq=0, fseq=10000;

  // Parse command line
  //
  cxxopts::Options options(prog, "Put together a PSP file from per-node components\n");

  options.add_options()
   ("h,help", "produce help message")
   ("v,verbose", "verbose output")
   ("c,checkpoint", "reassemble a checkpoint PSP file")
   ("d,wd", "working directory for input and output files",
     cxxopts::value<std::string>(inoutdir)->default_value("."))
   ("r,runtag", "EXP runtag name",
     cxxopts::value<std::string>(runtag)->default_value("run0"))
   ("p,prefix", "leading name of PSP output files",
     cxxopts::value<std::string>(prefix)->default_value("OUT"))
   ("s,seq", "SPL sequence counter",
     cxxopts::value<int>(seq)->default_value("0"))
   ("1,first", "initial index in SPL sequence",
     cxxopts::value<int>(bseq))
   ("2,last", "final index in SPL sequence",
     cxxopts::value<int>(fseq))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << argv[0]
	      << "-d outdir -r run -s 0" << std::endl;
    return 1;
  }

  if (vm.count("first")) {
    bseq = vm["first"].as<bool>();
    range = true;
  }

  if (vm.count("last")) {
    fseq = vm["last"].as<bool>();
    range = true;
  }

  if (vm.count("checkpoint")) {
    checkpt = true;
    range   = false;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

				// For magic number checking
				// -------------------------
  const static unsigned long magic = 0xadbfabc0;
  const static unsigned long mmask = 0xf;
  const static unsigned long nmask = ~mmask;

  if (not range) {		// Logic for writing a single PSP file
    bseq = seq;
    fseq = seq;
  }
				// Begin sequence range loop
				// -------------------------
  for (seq=bseq; seq<=fseq; seq++) {

				// Open the master file
				// --------------------
    std::ostringstream mfile;
    if (checkpt)
      mfile << inoutdir << "/SPL." << runtag << ".chkpt";
    else
      mfile << inoutdir << "/SPL." << runtag << "."
	    << std::setw(5) << std::setfill('0') << seq;
    
    std::ifstream master(mfile.str());

    if (master.fail()) {
      std::cerr << "spl2psp: can't open master file <" << mfile.str() 
		<< "> . . . quitting" << std::endl;
      exit(-1);
    }
    
				// Open the output file
				// --------------------
    std::ostringstream ofile;
    if (checkpt)
      ofile << inoutdir << "/" << prefix << "." << runtag << ".chkpt";
    else
      ofile << inoutdir << "/" << prefix << "." << runtag << "."
	    << std::setw(5) << std::setfill('0') << seq;
    
    std::ofstream wholef(ofile.str());
    
    if (wholef.fail()) {
      std::cerr << "spl2psp: can't open output file <" << ofile.str() 
		<< "> . . . quitting" << std::endl;
      exit(-1);
    }
				// Now read and write the header
				// -----------------------------
    struct MasterHeader header;
  
    master.read ((char *)      &header, sizeof(MasterHeader));
    wholef.write((const char *)&header, sizeof(MasterHeader));

    unsigned long cmagic;
    unsigned rsize;

    for (int c=0; c<header.ncomp; c++) {
      
      master.read((char*)&cmagic, sizeof(unsigned long));
      if ( (cmagic & nmask) != magic ) {
	std::string msg("Error identifying new PSP.  Is this an old PSP?");
	std::cerr << msg << std::endl;
      }
      rsize = cmagic & mmask;

      int nprocs;
      master.read((char*)&nprocs, sizeof(int));

      ComponentHeader header;
      if (!header.read(&master)) {
	std::string msg("ComponentHeader: error reading particle header");
	std::cerr << msg << std::endl;
	exit(-2);
      }

      wholef.write((const char*)&cmagic, sizeof(unsigned long));

      if (!header.write(&wholef)) {
	std::string msg("ComponentHeader: error writing particle header");
	std::cerr << msg << std::endl;
	exit(-3);
      }

				// Parse info field to get 
				// id and parameter strings
      YAML::Node config;

      std::istringstream sin(header.info.get());
      config = YAML::Load(sin);
      
      YAML::Node cconf;
      std::string name;

      try {
	name = config["name"].as<std::string>();
	cconf = config["parameters"];
      }
      catch (YAML::Exception & error) {
	if (myid==0) std::cout << "Error parsing YAML in PSP file: "
			       << error.what() << std::endl
			       << std::string(60, '-') << std::endl
			       << "Config node"        << std::endl
			       << std::string(60, '-') << std::endl
			       << config               << std::endl
			       << std::string(60, '-') << std::endl;
	exit(-4);
      }
      
      bool indexing = false;
      if (cconf["indexing"])
	indexing = cconf["indexing"].as<bool>();
      
				// File name buffer from master file
      const size_t PBUF_SIZ = 1024;
      char buf [PBUF_SIZ];
      
      int cnt = 0;
				// Read particles from each component file
      for (int n=0; n<nprocs; n++) {
	
	master.read((char*)buf, PBUF_SIZ);
	
	std::ifstream fcomp(inoutdir + "/" + buf);
	
	if (fcomp.fail()) {
	  std::cerr << "spl2psp: can't open component file <" << buf
		    << "> . . . quitting" << std::endl;
	  exit(-5);
	}

	// Get number of particles
	unsigned int N;
	fcomp.read((char*)&N, sizeof(unsigned int));

	// Read and write particle by particle
	Particle p;
	for (int k=0; k<N; k++) {
	  p.readBinary(rsize, indexing, cnt++, &fcomp);
	  p.writeBinary(rsize, indexing, &wholef);
	}
	// Loop over particles

      }
      // Loop over file chunk

      if (verbose)
	std::cout << cnt << " particles written in component <"
		  << name << ">" << std::endl;
    }
    // Loop over headers

    try {
      wholef.close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "spl2psp: exception closing file <" << ofile.str()
		<< ": " << e.what() << std::endl;
    }

    std::cout << "Successfully wrote PSP file <" << ofile.str() << ">"
	      << std::endl;
  }
  // Loop over sequence indices


  return 0;
}
  
