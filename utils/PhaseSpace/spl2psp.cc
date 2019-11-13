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
#include <boost/program_options.hpp>

namespace po = boost::program_options;

				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;
int multistep=0;

int
main(int argc, char **argv)
{
  char *prog = argv[0];
  std::string runtag, inoutdir, outfile;
  bool verbose = false;
  int seq;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("verbose,v",       "verbose output")
    ("wd,d",	        po::value<std::string>(&inoutdir)->default_value("."),
     "working directory")
    ("runtag,r",	po::value<std::string>(&runtag)->default_value("run0"),
     "component name")
    ("seq,s",		po::value<int>(&seq)->default_value(0),
     "sequence counter")
    ;


  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << argv[0]
	      << "-d outdir -r run -s 0" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

				// For magic number checking
				// -------------------------
  const static unsigned long magic = 0xadbfabc0;
  const static unsigned long mmask = 0xf;
  const static unsigned long nmask = ~mmask;

				// Open the master file
				// --------------------
  std::ostringstream mfile;
  mfile << inoutdir << "/SPL." << runtag << "." << std::setw(5) << std::setfill('0') << seq;

  std::ifstream master(mfile.str());

  if (master.fail()) {
    std::cerr << "psp2whole: can't open master file <" << mfile.str() 
	      << "> . . . quitting" << std::endl;
    exit(-1);
  }

				// Open the output file
				// --------------------
  std::ostringstream ofile;
  ofile << inoutdir << "/OUT." << runtag << "."
	<< std::setw(5) << std::setfill('0') << seq;

  std::ofstream wholef(ofile.str());

  if (wholef.fail()) {
    std::cerr << "psp2whole: can't open output file <" << ofile.str() 
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

    std::cout << "Read component header <" << name << ">" << std::endl;

    bool indexing = false;
    if (cconf["indexing"])
      indexing = cconf["indexing"].as<bool>();

    const size_t PBUF_SIZ = 1024;
    char buf [PBUF_SIZ];

    int seq = 0;
				// Read particles from each component file
    for (int n=0; n<nprocs; n++) {

	master.read((char*)buf, PBUF_SIZ);

	std::ifstream fcomp(buf);

	if (fcomp.fail()) {
	  std::cerr << "psp2whole: can't open component file <" << buf
		    << "> . . . quitting" << std::endl;
	  exit(-5);
	}


	unsigned int N;
	fcomp.read((char*)&N, sizeof(unsigned int));

	Particle p;
	for (int k=0; k<N; k++) {
	  p.readBinary(rsize, indexing, seq++, &fcomp);
	  p.writeBinary(rsize, indexing, &wholef);
	}
	// Loop over particles

	std::cout << seq << " particles written" << std::endl;

    }
    // Loop over file chunk

  }
  // Loop over headers

  return 0;
}
  
