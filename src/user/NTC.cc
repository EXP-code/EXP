#include <iostream>
#include <iomanip>
#include <sstream>

#include "Particle.H"
#include "pCell.H"
#include "NTC.H"

using namespace NTC;

// For verbose debugging
static const bool DEBUG_V   = false;

// Chatty output for debugging
bool NTCdb::chatty = false;

// Save interval
unsigned NTCdb::intvl = 10;

// Skip this number of calls between reports
unsigned NTCitem::skprpt    = 5;

// Maximum number of reports for each cell item
unsigned NTCitem::maxrpt    = 5;

// Minimum CrossSection x Velocity value
double   NTCitem::Min = 1.0e-24;

// Minimum CrossSection x Velocity value
double   NTCitem::Def = 10.0;

// Count live instances
unsigned NTCitem::instance  = 0;

void NTCitem::Test()
{
  inTest = true;		// To prevent recursive calling in Prob

				// Size of sample and separator bar
  const size_t number = 10;
  size_t wid = number * 12 + 10 + 14;

  if (db.size()==0)
    std::cout << "Empty caller: " << caller << " : " 
	      << " p: "  << &db
	      << std::endl;

  if (db.size()>0) {

    std::cout << std::string(70, '-') << std::endl
	      << "VelCrs structure: "  << this << " caller: " << caller 
	      << std::endl << std::string(70, '-') << std::endl;

    for (auto k : db) {
      sKeyPair p = k.first;
      std::ostringstream sout;
      sout << "<" << p.first.first  << "," << p.first.second
	   << "|" << p.second.first << "," << p.second.second << ">";

      std::cout << std::setw(14) << sout.str()
		<< std::setw(10) << k.second.count();

      size_t prc = std::cout.precision(2);
      for (size_t i=0; i<number; i++) {
	double P = (0.5 + i)/number;
	std::cout << std::setw(10) << k.second(P);
      }
      std::cout.precision(prc);
    }
    std::cout << std::string(wid, '-') << std::endl;
  } else {
    std::cout << std::string(wid, '-') << std::endl
	      << "VelCrs no points: " << this << std::endl
	      << std::string(wid, '-') << std::endl;
  }

  inTest = false;
}

NTCitem::vcMap NTCitem::Prob(double x)
{
  vcMap ret;
  
  // Deal with round off issues on resume
  //
  for (auto v : db) ret[v.first] = v.second.inverse(x);

  return ret;
}

double NTCitem::Prob(sKeyPair indx, double x)
{
  // Get stanza in db
  qpMap::iterator it = db.find(indx);

  // Default value
  if (it == db.end()) return 0.5;

  // Debug reporting
  //
  if (DEBUG_V and !inTest) {
    report++;			// Make sure the first one is skipped
    if (report % skprpt == 0 and report < maxrpt*skprpt)  Test();
  }

  // Return quantile
  return it->second.inverse(x);
}

double NTCitem::CrsVel(sKeyPair indx, double p)
{
  // Get stanza in db
  qpMap::iterator it = db.find(indx);

  // Default value
  if (it == db.end()) return Def;

  // Return value
  return it->second(p);
}

void NTCitem::Add(sKeyPair indx, double val)
{
  // value below threshold
  //
  if (val<=Min) return;

  // Check for initialization of Quantile
  //
  if (db.find(indx) == db.end()) {
    db[indx].histogram(Nequal);
    for (auto v : qs) db[indx].newQ(v);
  }

  // Add new element
  //
  db[indx].add(val);
}


NTCitem& NTCdb::operator[](const key_type& k)
{
  NTCdata::iterator it = data.find(k);

  if (it == data.end()) {
    
    // Look for an initialization parent
    //
    key_type check = k >> 3;
    it = data.find(check);
    while (it == data.end() and check != 0) {
      check = check >> 3;
      it = data.find(check);
    }

    // This is only informational, not a problem
    if (0) {
      if (DEBUG_V && it == data.end()) {
	std::cout << "NTC: no parent found for " << k << std::endl;
      }
    }

    // If none, create an empty item.  Otherwise, initialize from the
    // parent
    //
    if (it == data.end()) {
      data[k] = NTCitem();
      it = data.find(k);
    }
  }

  if (DEBUG_V) it->second.setKey(k);

  return it->second;
}
