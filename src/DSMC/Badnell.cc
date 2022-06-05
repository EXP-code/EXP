#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include <mpi.h>

#include <Badnell.H>

std::string BadnellData::datapath = "./";
bool        BadnellData::reweight = false;

BadnellData::BadnellData()
{
  // Ion registry (just He so far)
  //
  ions = { {2, "he"} };		// map of mnemonics
  ionQ = { {2, 2}, {2, 3} };	// vector of ions
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

    // File filters
    std::ostringstream m1, m2;
    m1 << user << year;
    m2 << ions[ZC.first] << ZC.second - 1;

    auto filter = [&m1, &m2](const std::string& name) {
		    return (name.find(m1.str()) != std::string::npos &&
			    name.find(m2.str()) != std::string::npos );
		  };

    // Will contain matched file(s)
    std::vector<std::string> names;

    try {
      for (auto const& dir_entry : std::filesystem::directory_iterator{dpath / "RR"}) {
	std::string name = dir_entry.path().string();
	if (filter(name)) names.push_back(name);
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
	if (filter(name)) names.push_back(name);
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
