#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <regex>
#include <map>
#include <set>

#include <mpi.h>

#include <atomic_constants.H>
#include <Badnell.H>

std::string BadnellData::datapath = "./";
std::string BadnellData::coupling = "ic";
bool        BadnellData::reweight = false;

BadnellData::BadnellData()
{
  std::tie(ions, ionQ) = walkDirectory();
}

std::tuple<std::map<int, std::string>, std::vector<BadnellData::lQ>>
BadnellData::walkDirectory()
{
  // Periodic table
  //
  PeriodicTable table;

  // Filter strings
  //
  const std::string user("nrb"), year("20");

  // Set default datapath
  //
  std::filesystem::path dpath{datapath};

  // Data list
  //
  std::set<std::pair<unsigned short, unsigned short>> ionq;
  std::map<int, std::string> ions;

  // Look for atomic data path in environment
  //
  if (const char* env_p = std::getenv("ATOMIC_DATA_PATH")) {
    dpath = std::filesystem::path( env_p ) / "RR";

    const std::regex txt_regex("[a-z]+");
    const std::regex num_regex("[0-9]+");
    std::smatch matches;
    std::string element, charge, couple;

    for (auto const& dir_entry : std::filesystem::recursive_directory_iterator{dpath}) {
      std::string fname = dir_entry.path().string();
      if (fname.find(user+year) != std::string::npos) {
	// Break into data elements
	int beg = fname.find('#') + 1;
	int end = fname.find('_');

	// Get sequence
	std::string seq = fname.substr(beg, end - beg);

	// Get the element
	//
	fname = fname.substr(end+1, std::string::npos);
	if (std::regex_search(fname, matches, txt_regex))
	  element = matches[0];
	else {
	  std::runtime_error("No element match in: " + fname);
	}

	// Get the charge
	//
	fname = fname.substr(element.size(), std::string::npos);
	if (std::regex_search(fname, matches, num_regex))
	  charge = matches[0];
	else {
	  std::runtime_error("No charge match in: " + fname);
	}

	// Get the couple
	//
	fname = fname.substr(charge.size(), std::string::npos);
	if (std::regex_search(fname, matches, txt_regex))
	  couple = matches[0];
	else {
	  std::runtime_error("No couple match in: " + fname);
	}

	// Ignore unwanted coupling
	//
	if (couple != coupling) continue;

	std::string abbrev(element);
	abbrev[0] = std::toupper(abbrev[0]);

	auto ret = table[abbrev];

	std::istringstream sin(charge);
	unsigned short chg;
	sin >> chg;

	/*
	std::cout << std::setw( 8) << seq
		  << std::setw( 8) << element
		  << std::setw(20) << std::get<0>(*ret)
		  << std::setw( 8) << std::get<2>(*ret)
		  << std::setw( 8) << chg
		  << std::setw( 8) << couple
		  << std::endl;
	*/

	ions[std::get<2>(*ret)] = element;
	ionq.insert({std::get<2>(*ret), chg+1});
      }
    }
  }
  
  std::vector<std::pair<unsigned short, unsigned short>> ionQ;
  for (auto s : ionq) {
    ionQ.push_back(s);
  }

  return {ions, ionQ};
}


void BadnellData::initialize(atomicData* ad)
{
  // Attempt to open files for ions in ionQ
  //

  /** From Nigel

      The naming convention is:
      user,year-code,sequence,element,charge,coupling scheme,core
      excitation (n->n').

      The year code is when the approach was introduced (NOT when
      calculated).  So, "20" is a new year code for the binned cross
      sections.  And "00" is the old year code for j-resolved rate
      coeffs.
  */

  // Filter strings
  //
  const std::string user("nrb"), year("20");

  // Set default datapath
  //
  std::filesystem::path dpath{datapath};

  // Look for atomic data path in environment
  //
  if (const char* env_p = std::getenv("ATOMIC_DATA_PATH")) {
    dpath = std::filesystem::path( env_p );
  }
  
  // Loop through ion list and read data
  //
  for (auto ZC : ionQ) {

    // File regular expression
    //
    std::ostringstream rstr;

    // Path--+               Sequence-----+     +---Separator
    //       |                            |     |
    //       |               Separator--+ |     |           +--Core excitation
    //       |                          | |     |           |  n->n' for
    //       v                          v v     v           |  dielectronic
    rstr << ".*[/]" << user << year << "#[a-z]+[_]" //      v         
	 << ions[ZC.first] << ZC.second - 1 << coupling << "[0-9]*[.]dat";
    //        ^               ^                ^                     ^
    //        |               |                |                     |
    // Ion name               Charge           Coupling     Suffix---+

    std::regex species(rstr.str());
    std::smatch matches;

    // Will contain matched file(s)
    std::vector<std::string> names;

    try {
      for (auto const& dir_entry : std::filesystem::recursive_directory_iterator{dpath / "RR"}) {
	std::string name = dir_entry.path().string();
	if (std::regex_search(name, matches, species))
	  names.push_back(name);
      }
    }
    catch(std::filesystem::filesystem_error const& ex) {
      std::cout
	<< "Directory scan error for RR..."
	<< "what():  " << ex.what() << '\n'
	<< "path1(): " << ex.path1() << '\n'
	<< "path2(): " << ex.path2() << '\n'
	<< "code().value():    " << ex.code().value() << '\n'
	<< "code().message():  " << ex.code().message() << '\n'
	<< "code().category(): " << ex.code().category().name() << '\n';
      
      MPI_Finalize();
      exit(32);
    }
    
    if (names.size() != 1) {
      std::cout << "Unexpected RR matches: " << names.size() << std::endl;
      for (auto s : names) std::cout << s << " ";
      std::cout << std::endl;
      MPI_Finalize();
      exit(33);
    }

    std::ifstream in{dpath / names[0]};
    
    std::string line;
    std::getline(in, line);

    double E, X;
    
    // Create a data record
    bdPtr d = std::make_shared<BadnellRec>();
    
    bool looking = true;
    while (in) {
      if (looking) {
	// Keep reading lines until we reach the totals
	if (line.find("E(RYD)") != std::string::npos) {
	  looking = false;
	  std::getline(in, line); // Read header
	}
      } else {
	// We are done!
	if (line[0] == 'C') break;
	
	// Read the line
	std::istringstream ins(line);
	
	// Values (Energy in Rydbergs and Energy*Sigma in Mbarn*Rydberg)
	ins >> E >> X;
	d->E_rr.push_back(E*RydtoeV);
	d->X_rr.push_back(X*RydtoeV);
      }
      
      // Read next line
      std::getline(in, line);
    }
    // END: RR file read
    
    // Now read DR
    //

    names.clear();
    try {
      for (auto const& dir_entry : std::filesystem::directory_iterator{dpath / "DR"}) {
	std::string name = dir_entry.path().string();
	if (std::regex_search(name, matches, species)) {
	  names.push_back(name);
	}
      }
    }
    catch(std::filesystem::filesystem_error const& ex) {
      std::cout
	<< "Directory scan error for DR..."
	<< "what():  " << ex.what() << '\n'
	<< "path1(): " << ex.path1() << '\n'
	<< "path2(): " << ex.path2() << '\n'
	<< "code().value():    " << ex.code().value() << '\n'
	<< "code().message():  " << ex.code().message() << '\n'
	<< "code().category(): " << ex.code().category().name() << '\n';
      
      MPI_Finalize();
      exit(34);
    }

    if (names.size()) {

      if (names.size() != 1) {
	std::cout << "Unexpected DR matches: " << names.size() << std::endl;
	for (auto s : names) std::cout << s << " ";
	std::cout << std::endl;
	MPI_Finalize();
	exit(35);
      }

      in = std::ifstream{dpath / "DR" / names[0]};
      
      std::getline(in, line);
    
      looking = true;
      while (in) {
	if (looking) {
	  // Keep reading lines until we reach the totals
	  if (line.find("E(RYD)") != std::string::npos) {
	    looking = false;
	    std::getline(in, line); // Read header
	  }
	} else {
	  // We are done!
	  if (line[0] == 'C') break;
	  
	  // Read the line
	  std::istringstream ins(line);
	  
	  // Values (Energy in Rydberg and energy averaged cross
	  // section in Mbarn)
	  ins >> E >> X;
	  d->E_dr.push_back(E*RydtoeV); // Energy is now in eV
	  d->X_dr.push_back(X);
	}
	
	// Read next line
	std::getline(in, line);
      }
      // END: DR file read

      // Reweight energy average
      //
      if (reweight) {
	for (int i=0; i<d->E_dr.size()-1; i++) {
	  if (i<d->X_dr.size())
	    d->X_dr[i] *= 2.0/(1.0 + d->E_dr[i]/d->E_dr[i+1]);
	  else
	    std::cout << "DR error for [" << ZC.first << ", " << ZC.second
		      << "]: i=" << i << " !< " << d->X_dr.size() << std::endl;
	}
      } // END: reweight
      
    } // END: DR will not exist for hydrogenic case

    data[ZC] = d;
  }

  // DONE
}
