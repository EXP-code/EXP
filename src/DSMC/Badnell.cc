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
  // Periodic table
  //
  PeriodicTable table;

  /** From Nigel

      The naming convention is:
      user,year-code,sequence,element,charge,coupling scheme,core
      excitation (n->n').

      The year code is when the approach was introduced (NOT when
      calculated).  So, "20" is a new year code for the binned cross
      sections.  And "00" is the old year code for j-resolved rate
      coeffs.
  */

  // Set default datapath
  //
  std::filesystem::path dpath{datapath};

  // Look for atomic data path in environment
  //
  if (const char* env_p = std::getenv("ATOMIC_DATA_PATH")) {
    dpath = std::filesystem::path( env_p );
  }
  
  // Filter strings
  //
  const std::string user("nrb"), year("20"); // RR 


  // File regular expression for RR
  //
  std::ostringstream rstr;

  // Path--+               Sequence-----+     +---Separator
  //       |                            |     |
  //       |               Separator--+ |     |
  //       |                          | |     |
  //       v                          v v     v
  rstr << ".*[/]" << user << year << "#[a-z]+[_]"
       << "([a-z]+)([0-9]+)" << coupling << "[0-9]*[.]dat";
  //        ^       ^           ^            ^        ^
  //        |       |           |            |        |
  //        |       |           |            +--------|----- Core excitation
  //        |       |           |                     |      n->n' for
  // Ion name       Charge      Coupling     Suffix---+      dielectronic
  // (Group 1)      (Group 2)

  std::regex species(rstr.str());
  std::smatch matches;

  // Will contain matched file(s).  Should only be one.
  //
  std::map<lQ, std::string> names;

  try {
    for (auto const& dir_entry : std::filesystem::recursive_directory_iterator{dpath / "RR"}) {
      std::string name = dir_entry.path().string();
      if (std::regex_search(name, matches, species)) {
	auto it = table[matches[1].str()];
	if (it == 0) {
	  std::ostringstream sout;
	  sout << "Badnell::initialize: [0] could not locate element <"
	       << matches[1].str() << ">";
	  throw std::runtime_error(sout.str());
	}

	// The element and and charge
	//
	lQ ZC = {it->number(), stoi(matches[2].str()) + 1};

	// Add the file name
	//
	names[ZC] = name;
      }
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
    
  // File regular expressions for DR
  //
  std::ostringstream rdr1, rdr2;
  std::string yr1("00"), yr2("12");

  // Path--+              Sequence-----+     +---Separator
  //       |                           |     |
  //       |              Separator--+ |     |
  //       |                         | |     |
  //       v                         v v     v
  rdr1 << ".*[/]" << user << yr1 << "#[a-z]+[_]"
       << "([a-z]+)([0-9]+)" << coupling << "([0-9]*)_t5.dat";
  rdr2 << ".*[/]" << user << yr2 << "#[a-z]+[_]"
       << "([a-z]+)([0-9]+)" << coupling << "([0-9]*)_t5.dat";
  //        ^       ^           ^             ^          ^
  //        |       |           |             |          |
  //        |       |           |             +----------|-- Core excitation
  //        |       |           |                        |   n->n' for
  // Ion name       Charge      Coupling        Suffix---+   dielectronic
  // (Group 1)      (Group 2)                                (Group 3)

  std::regex spc1(rdr1.str());
  std::regex spc2(rdr2.str());

  std::smatch match1, match2;

  std::map<lQ, std::vector<std::pair<std::string, std::string>>> names2;

  try {
    for (auto const& dir_entry : std::filesystem::recursive_directory_iterator{dpath / "DR"}) {
      std::string name = dir_entry.path().string();
      if (std::regex_search(name, match1, spc1)) {
	auto it = table[match1[1].str()];
	if (it == 0) {
	  std::ostringstream sout;
	  sout << "Badnell::initialize: [1] could not locate element <"
	       << match1[1].str() << ">";
	  throw std::runtime_error(sout.str());
	}

	// The element and and charge
	//
	lQ ZC = {it->number(), stoi(match1[2].str()) + 1};

	// Add this excitaton to the db
	//
	names2[ZC].push_back({match1[3].str(), name});
      }
      if (std::regex_search(name, match2, spc2)) {
	auto it = table[match2[1].str()];
	if (it == 0) {
	  std::ostringstream sout;
	  sout << "Badnell::initialize: [2] could not locate element <"
	       << match2[1].str() << ">";
	  throw std::runtime_error(sout.str());
	}

	// The element and and charge
	//
	lQ ZC = {it->number(), stoi(match2[2].str()) + 1};

	// Add this excitaton to the db
	//
	names2[ZC].push_back({match1[3].str(), name});
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

  // Loop through each desired ion and process the files
  //
  for (auto ZC : ionQ) {

    // Create a data record
    //
    bdPtr d = std::make_shared<BadnellRec>();

    // Begin with RR files
    //
    auto it = names.find(ZC);

    // Have a file
    //
    if (it != names.end()) {

      std::ifstream in{dpath / it->second};
    
      std::string line;
      std::getline(in, line);

      double E, X;
    
      bool looking = true;
      while (in) {
	if (looking) {
	  // Keep reading lines until we reach the totals
	  //
	  if (line.find("E(RYD)") != std::string::npos) {
	    looking = false;
	    std::getline(in, line); // Read header
	  }
	} else {
	  // We are done!
	  //
	  if (line[0] == 'C') break;
	  
	  // Read the line
	  //
	  std::istringstream ins(line);
	  
	  // Values (Energy in Rydbergs and Energy*Sigma in Mbarn*Rydberg)
	  //
	  ins >> E >> X;
	  d->E_rr.push_back(E*RydtoeV);
	  d->X_rr.push_back(X*RydtoeV);
	}
	
	// Read next line
	std::getline(in, line);
      }
      // END: RR file read

      //       +--- Verbose reporting on RR files found
      //       |    (not for production mode)
      //       v
    } else  if (false) {
      std::ostringstream sout;
      sout << "Badnell::initialize: could not locate RR file for ion (Z, C) = ("
	   << ZC.first << ", " << ZC.second << ")";
      std::cout << sout.str() << std::endl;
    }
    
    // Now, on to DR files
    //
    auto it2 = names2.find(ZC);

    if (it2 != names2.end()) {
      
      for (auto v : it2->second) {

	auto cr = v.first;
	// Test for uncertain parent-line files
	if (cr.size()>2) continue;

	auto in = std::ifstream{dpath / "DR" / v.second};
	
	std::string line;
	std::getline(in, line);
	
	double E, X;
	bool looking = true;
	while (in) {
	  if (looking) {
	    // Keep reading lines until we reach the totals
	    //
	    if (line.find("E(RYD)") != std::string::npos) {
	      looking = false;
	      std::getline(in, line); // Read header
	    }
	  } else {
	    // We are done!
	    //
	    if (line[0] == 'C') break;
	    
	    // Read the line
	    //
	    std::istringstream ins(line);
	    
	    // Default values
	    //
	    E = X = 0.0;
	    
	    // Values (Energy in Rydberg and energy averaged cross
	    // section in Mbarn)
	    ins >> E >> X;
	    d->E_dr[cr].push_back(E*RydtoeV); // Energy is now in eV
	    d->X_dr[cr].push_back(X);
	  }
	
	  // Read next line
	  std::getline(in, line);
	}
	// END: DR file read
	
	// Reweight energy average
	//
	if (reweight) {
	  for (int i=0; i<d->E_dr[cr].size()-1; i++) {
	    if (i<d->X_dr[cr].size())
	      d->X_dr[cr][i] *= 2.0/(1.0 + d->E_dr[cr][i]/d->E_dr[cr][i+1]);
	    else
	      std::cout << "DR error for [" << ZC.first << ", " << ZC.second
			<< "]: i=" << i << " !< " << d->X_dr[cr].size() << std::endl;
	  }
	}
	// END: reweight
	
      }
      // END: DR stanza

      //       +--- Verbose reporting on DR files found
      //       |    (not for production mode)
      //       v
    } else if (false) {
      std::ostringstream sout;
      sout << "Badnell::initialize: could not locate DR file for element "
	   << table[ZC.first]->name() << ", ion (Z, C) = ("
	   << ZC.first << ", " << ZC.second << ")";
      std::cout << sout.str() << std::endl;
    }
    // END: file loop

    data[ZC] = d;
  }

  // DONE
}

void BadnellData::report()
{
  std::cout << std::string(64, '-') << std::endl
	    << "ADAS cross section check" << std::endl
	    << std::string(64, '-') << std::endl
	    << std::setw(8) << "Z"
	    << std::setw(8) << "C"
	    << std::setw(8) << "N(RR)"
	    << std::setw(8) << "core"
	    << std::setw(8) << "N(DR)"
	    << std::endl
	    << std::setw(8) << "-----"
	    << std::setw(8) << "-----"
	    << std::setw(8) << "-----"
	    << std::setw(8) << "-----"
	    << std::setw(8) << "-----"
	    << std::endl;

  for (auto d : data) {
    std::cout << std::setw(8) << d.first.first
	      << std::setw(8) << d.first.second
	      << std::setw(8) << d.second->E_rr.size()
	      << std::endl;
    for (auto s : d.second->E_dr) {
      std::cout << std::setw(8) << ""
		<< std::setw(8) << ""
		<< std::setw(8) << ""
		<< std::setw(8) << s.first
		<< std::setw(8) << s.second.size() << std::endl;
    }
  }
  std::cout << std::string(32, '-') << std::endl;
}

