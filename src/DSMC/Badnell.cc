#include <iostream>
#include <fstream>
#include <sstream>

#include <mpi.h>

#include <Badnell.H>
#include <scandir.H>

#include <boost/make_shared.hpp>

std::string BadnellData::RRdir = "./RR";
std::string BadnellData::DRdir = "./DR";
bool BadnellData::reweight = true;

BadnellData::BadnellData()
{
  // Ion registry (just He so far)
  //
  ions = { {2, "he"} };		// map of mnemonics
  ionQ = { {2, 2}, {2, 3} };	// vector of ions
}

void BadnellData::initialize(chdata* ch)
{
  // Attempt to open files for existing
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

  const std::string user("nrb"), year("20");

  for (auto ZC : ionQ) {

    // File filters
    std::ostringstream m1, m2;
    m1 << user << year;
    m2 << ions[ZC.first] << ZC.second - 1;

    auto filter = [&m1, &m2](const std::string& name) {
		    return (name.find(m1.str()) != std::string::npos &&
			    name.find(m2.str()) != std::string::npos );
		  };

    std::cout << "[before] Filter=" << m1.str() << ", " << m2.str() << std::endl;

    std::vector<std::string> names = scandirpp::get_names(RRdir, filter);

    std::cout << "[RR] Filter=" << m1.str() << ", " << m2.str()
	      << " :: " << RRdir + "/" + names[0] << std::endl;

    if (names.size() != 1) {
      std::cout << "Unexpected RR matches: " << names.size() << std::endl;
      MPI_Finalize();
      exit(33);
    }

    std::ifstream in(RRdir + "/" + names[0]);
    
    std::string line;
    std::getline(in, line);

    double E, X;
    
    // Create a data record
    bdPtr d = boost::make_shared<BadnellRec>();

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
	
	// Values
	ins >> E >> X;
	d->E_rr.push_back(E);
	d->X_rr.push_back(X);
      }
      
      // Read next line
      std::getline(in, line);
    }
    
    // Reweight energies
    //
    if (reweight) {
      for (int i=0; i<d->E_rr.size()-1; i++) {
	if (i<d->X_rr.size())
	  d->X_rr[i] *= 2.0/(1.0 + d->E_rr[i]/d->E_rr[i+1]);
	else {
	  std::cout << "RR error for [" << ZC.first << ", " << ZC.second
		    << "]: i=" << i << " !< " << d->X_rr.size() << std::endl;
	}
      }
    }

    // Now read DR

    names = scandirpp::get_names(DRdir, filter);

    if (names.size()) {

      if (names.size() != 1) {
	std::cout << "Unexpected DR matches: " << names.size() << std::endl;
	MPI_Finalize();
	exit(34);
      }

      std::cout << "[DR] Filter=" << m1.str() << ", " << m2.str()
		<< " :: " << DRdir + "/" + names[0] << std::endl;
      in.close();
      in.open(DRdir + "/" + names[0]);
      
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
	  
	  // Values
	  ins >> E >> X;
	  d->E_dr.push_back(E);
	  d->X_dr.push_back(X);
	}
	
	// Read next line
	std::getline(in, line);
      }

      if (reweight) {
	for (int i=0; i<d->E_dr.size()-1; i++) {
	  if (i<d->X_dr.size())
	    d->X_dr[i] *= 2.0/(1.0 + d->E_dr[i]/d->E_dr[i+1]);
	  else
	    std::cout << "DR error for [" << ZC.first << ", " << ZC.second
		      << "]: i=" << i << " !< " << d->X_dr.size() << std::endl;
	}
      }

    } // END: DR will not exist for hydrogenic case

    data[ZC] = d;
  }

}
