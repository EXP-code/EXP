/*
  Compute min and max for all fields

  MDWeinberg 08/26/11, 11/24/19
*/

using namespace std;

#include <cstdlib>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>

#include <Species.H>

#include <header.H>
#include <PSP.H>
#include <InitContainer.H>

#include <boost/program_options.hpp>

namespace po = boost::program_options;


				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;

int
main(int ac, char **av)
{
  bool verbose = false;
  char *prog = av[0];
  std:: string cname, new_dir("./");
  double time;
  int sindx;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("verbose,v",       "verbose output")
    ("species,s",	po::value<int>(&sindx)->default_value(-1),
     "position of species index")
    ("name,c",	        po::value<std::string>(&cname)->default_value("comp"),
     "component name")
    ("files,f",         po::value< std::vector<std::string> >(), 
     "input files")
    ("dir,d",           po::value<std::string>(&new_dir), 
     "replacement SPL file directory")
    ;


  po::variables_map vm;

  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << av[0]
	      << " -f OUT.run.00001 -s 1 -c gas" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  std::vector<std::string> files = vm["files"].as< std::vector<std::string> >();

  bool first = true;
  
  for (auto file : files ) {

    if (verbose) cerr << "Using filename: " << file << endl;


				// Parse the PSP file
				// ------------------
    PSPptr psp;
    if (file.find("SPL") != std::string::npos)
      psp = std::make_shared<PSPspl>(file, new_dir);
    else
      psp = std::make_shared<PSPout>(file);

				// Now write a summary
				// -------------------
    if (verbose) {
      
      psp->PrintSummary(cerr);
    
      cerr << "\nPSP file <" << file << "> has time <" 
	   << psp->CurrentTime() << ">\n";
    }

				// Dump ascii for each component
				// -----------------------------
    vector<double> pos(3), vel(3);
    int itmp, icnt, iv;

				// Initialize the value maps
				// -------------------------
    PSPstanza *stanza;
    SParticle* part;

    std::set< speciesKey >                       sofar;

    std::map< speciesKey, std::vector<double> >  posmin, posmax;
    std::map< speciesKey, std::vector<double> >  velmin, velmax;
    std::map< speciesKey, float >                phimin, phimax, masmin, masmax;
    std::map< speciesKey, std::vector<int> >     iavmin, iavmax;
    std::map< speciesKey, std::vector<double> >  davmin, davmax;

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;

      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

	double mass = part->mass();

	speciesKey key = {0, 0};

	if (sindx >= 0) {
	  key = KeyConvert(part->iatr(sindx)).getKey();
	}
	
				// Initialize maps?
	if (sofar.find(key) == sofar.end()) {
	  posmin[key].resize(3,  FLT_MAX);
	  posmax[key].resize(3, -FLT_MAX);
	  velmin[key].resize(3,  FLT_MAX);
	  velmax[key].resize(3, -FLT_MAX);
	  masmin[key] =  FLT_MAX;
	  masmax[key] =  0.0;
	  phimin[key] =  FLT_MAX;
	  phimax[key] = -FLT_MAX;
	  iavmin[key].resize(part->niatr(),  INT_MAX);
	  iavmax[key].resize(part->niatr(), -INT_MAX);
	  davmin[key].resize(part->ndatr(),  FLT_MAX);
	  davmax[key].resize(part->ndatr(), -FLT_MAX);
	  sofar.insert(key);
	}
	
	masmin[key] = std::min<double>(part->mass(), masmin[key]);
	masmax[key] = std::max<double>(part->mass(), masmax[key]);

	for (size_t k=0; k<3; k++) {
	  posmin[key][k] = std::min<double>(part->pos(k), posmin[key][k]);
	  posmax[key][k] = std::max<double>(part->pos(k), posmax[key][k]);
	  velmin[key][k] = std::min<double>(part->vel(k), velmin[key][k]);
	  velmax[key][k] = std::max<double>(part->vel(k), velmax[key][k]);
	}
	
	phimin[key] = std::min<double>(part->phi(), phimin[key]);
	phimax[key] = std::max<double>(part->phi(), phimax[key]);

	for (int k=0; k<part->niatr(); k++) {
	  iavmin[key][k] = std::min<int>(part->iatr(k), iavmin[key][k]);
	  iavmax[key][k] = std::max<int>(part->iatr(k), iavmax[key][k]);
	}

	for (int k=0; k<part->ndatr(); k++) {
	  davmin[key][k] = std::min<double>(part->datr(k), davmin[key][k]);
	  davmax[key][k] = std::max<double>(part->datr(k), davmax[key][k]);
	}
      }
    }
    
    if (sofar.size()) {

      //
      // Output
      //
      const size_t fw = 18;
      double Time = psp->CurrentTime();

      std::cout << "Time=" << Time << std::endl << std::endl;

      std::cout << std::string(4*fw, '-') << std::endl;

      std::cout << std::right;
      if (sindx>=0) std::cout << std::setw(fw) << "Key";
      std::cout << std::setw(fw) << "Value"
		<< std::setw(fw) << "Min"
		<< std::setw(fw) << "Max"
		<< std::endl;
      if (sindx>=0) std::cout << std::setw(fw) << "------";
      std::cout << std::setw(fw) << "------"
		<< std::setw(fw) << "------"
		<< std::setw(fw) << "------"
		<< std::endl;
      
      for (auto key : sofar) {
	std::ostringstream skey;
	skey << "<" << key.first << ", " << key.second << ">";
	
	std::cout << std::right;
	if (sindx >= 0) std::cout << std::setw(fw) << skey.str();
	std::cout << std::setw(fw) << "Mass"
		  << std::setw(fw) << masmin[key]
		  << std::setw(fw) << masmax[key]
		  << std::endl;

	for (int k=0; k<3; k++) {
	  std::ostringstream snam;
	  snam << "pos(" << k+1 << ")";
	  std::cout << std::right;
	  if (sindx >= 0) std::cout << std::setw(fw) << skey.str();
	  std::cout << std::setw(fw) << snam.str()
		    << std::setw(fw) << posmin[key][k]
		    << std::setw(fw) << posmax[key][k]
		    << std::endl;
	}

	for (int k=0; k<3; k++) {
	  std::ostringstream snam;
	  snam << "vel(" << k+1 << ")";
	  std::cout << std::right;
	  if (sindx >= 0) std::cout << std::setw(fw) << skey.str();
	  std::cout << std::setw(fw) << snam.str()
		    << std::setw(fw) << velmin[key][k]
		    << std::setw(fw) << velmax[key][k]
		    << std::endl;
	}
	
	std::cout << std::right;
	if (sindx >= 0) std::cout << std::setw(fw) << skey.str();
	std::cout << std::setw(fw) << "Phi"
		  << std::setw(fw) << phimin[key]
		  << std::setw(fw) << phimax[key]
		  << std::endl;
	
	
	for (size_t k=0; k<iavmin[key].size(); k++) {
	  std::ostringstream snam;
	  snam << "iatr(" << k << ")";

	  std::cout << std::right;
	  if (sindx >= 0) std::cout << std::setw(fw) << skey.str();
	  std::cout << std::setw(fw) << snam.str()
		    << std::setw(fw) << iavmin[key][k]
		    << std::setw(fw) << iavmax[key][k]
		    << std::endl;
	}
	
	for (size_t k=0; k<davmin[key].size(); k++) {
	  std::ostringstream snam;
	  snam << "datr(" << k << ")";
	  std::cout << std::right;
	  if (sindx >= 0) std::cout << std::setw(fw) << skey.str();
	  std::cout << std::setw(fw) << snam.str()
		    << std::setw(fw) << davmin[key][k]
		    << std::setw(fw) << davmax[key][k]
		    << std::endl;
	}

	std::cout << std::string(4*fw, '-') << std::endl;
      }

    } else {
      std::cout << "No data . . . does component <" << cname
		<< "> exist?" << std::endl;
    }
    cout << endl;
  }

  return 0;
}
