#include <iostream>
#include <fstream>
#include <sstream>

#include <Badnell.H>

#include <boost/make_shared.hpp>

std::string BadnellData::RRdir = "./RR";
std::string BadnellData::DRdir = "./DR";

BadnellData::BadnellData()
{
  // Ion registry
  //
  ions = { {2, "he"} };		// map of mnemonics
  ionQ = { {2, 2} };		// vector of ions
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

  const std::string user("nrb"), year("20"), seq("#h_"), coupling("ic"), core("12");

  for (auto ZC : ionQ) {

    // First read RR
    std::ostringstream file;
    file << RRdir << "/" << user << year << seq << ions[ZC.first] << ZC.second - 1 << coupling << ".dat";
    std::cout << "File=" << file.str() << std::endl;
    std::ifstream in(file.str());

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
    
    // Now read DR
    file.str("");
    file << DRdir << "/" << user << year << seq << ions[ZC.first] << ZC.second - 1 << coupling << core << ".dat";
    std::cout << "File=" << file.str() << std::endl;
    in.close();
    in.open(file.str());

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

    data[ZC] = d;
  }

}
