#include <iostream>
#include <iomanip>
#include <sstream>

#include "Particle.H"
#include "pCell.H"
#include "NTC.H"

// For verbose debugging
static const bool DEBUG_V = true;

// Skip this number of calls between reports
unsigned NTCitem::skprpt    = 5;

// Maximum number of reports for each cell item
unsigned NTCitem::maxrpt    = 5;

// Maximum number of cached relative velocities
size_t   NTCitem::VelCrsSZ  = 4096; 

// Minimum CrossSection x Velocity value
double   NTCitem::VelCrsMin = 1.0e-24;

// Minimum CrossSection x Velocity value
double   NTCitem::VelCrsDef = 10.0;

void NTCitem::VelCrsTest()
{
  inTest = true;

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

      vcTup q1(VelCrsAvg(p, 0.25));
      vcTup q2(VelCrsAvg(p, 0.50));
      vcTup q3(VelCrsAvg(p, 0.75));

      vcTup p1( std::get<0>(q1),  std::get<0>(q2), std::get<0>(q3) );
      vcTup p2( std::get<1>(q1),  std::get<1>(q2), std::get<1>(q3) );
      vcTup p3( std::get<2>(q1),  std::get<2>(q2), std::get<2>(q3) );

      std::cout << std::setw(14) << sout.str()
		<< std::setw(10) << k.second.size()
		<< "  " << p1 << "  product" << std::endl 
		<< std::string(70, '-')      << std::endl
		<< std::setw(24) << ""
		<< "  " << p2 << "  cross"   << std::endl
		<< std::string(70, '-')      << std::endl
		<< std::setw(24) << ""
		<< "  " << p3 << "  ratio"   << std::endl
		<< std::string(70, '-')      << std::endl;
    }
    std::cout << std::string(70, '-') << std::endl;
  } else {
    std::cout << std::string(70, '-') << std::endl
	      << "VelCrs no points: " << this << std::endl
	      << std::string(70, '-') << std::endl;
  }

  inTest = false;
}

NTCitem::vcMap NTCitem::VelCrsAvg(double quant)
{
  vcMap ret;

  // Sanity
  quant = std::min(std::max(0.0, quant), 0.9999);

  for (auto v : db) {
    // Sort for quantile
    std::sort(v.second.begin(), v.second.end());
    //        ^
    //        |
    //        +--- This is a copy

    // Assign quantile
    //
    ret[v.first] = v.second[floor(v.second.size()*quant)];
  }

  return ret;
}

NTCitem::vcTup NTCitem::VelCrsAvg(sKeyPair indx, double quant)
{
  // Get stanza in db
  dqMap::iterator it = db.find(indx);

  // Default value
  if (it == db.end()) return vcTup(VelCrsDef, 0, 0);

  // Sanity
  quant = std::min(std::max(0.0, quant), 0.9999);

  // Sort for quantile
  dqTup tmp(it->second); std::sort(tmp.begin(), tmp.end());
  //        ^
  //        |
  //        +--- This is NOT a copy
  
  // Debug reporting
  //
  if (DEBUG_V and !inTest) {
    report++;			// Make sure the first one is skipped
    if (report % skprpt == 0 and report < maxrpt*skprpt)  VelCrsTest();
  }

  // Return quantile
  return tmp[floor(tmp.size()*quant)];
}

void NTCitem::VelCrsAdd(dqMap& vals)
{
  for (auto it : vals) {
    db[it.first].insert(db[it.first].begin(), it.second.begin(), it.second.end());
    if (db[it.first].size() > VelCrsSZ) db[it.first].resize(VelCrsSZ);
  }
}

void NTCitem::VelCrsAdd(sKeyPair indx, const vcTup& val)
{
  // value below threshold
  //
  if (std::get<0>(val)<=VelCrsMin) return;

  // Add new element
  //
  db[indx].push_front(val);

  // Push out the oldest element if deque is full
  //
  if (db[indx].size()>VelCrsSZ) db[indx].resize(VelCrsSZ);

}


NTCptr NTCdb::operator[](const key_type& k)
{
  NTCdata::iterator it = data.find(k);
  NTCptr ret;

  if (it == data.end()) {
    
    // Look for an initialization parent
    //
    key_type check = k >> 3;
    it = data.find(check);
    while (it == data.end() and check != 0) {
      check = check >> 3;
      it = data.find(check);
    }

    if (DEBUG_V && it == data.end()) {
      std::cout << "NTC: no parent found for " << k << std::endl;
    }

    // If none, create an empty item.  Otherwise, initialize from the
    // parent
    //
    if (it == data.end()) {
      data[k] = ret = NTCptr(new NTCitem);
    } else {
      data[k] = ret = NTCptr(new NTCitem(it->second));
    }

  } else {
    ret = it->second;
  }

  if (DEBUG_V) ret->setKey(k);

  return ret;
}
