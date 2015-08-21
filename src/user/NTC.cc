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
double   NTCitem::VelCrsMin = 1.0e-24;

// Minimum CrossSection x Velocity value
double   NTCitem::VelCrsDef = 10.0;

// Count live instances
unsigned NTCitem::instance  = 0;

void NTCitem::VelCrsTest()
{
  inTest = true;		// To prevent recursive calling in VelCrsAvg

  if (db.size()==0)
    std::cout << "Empty caller: " << caller << " : " 
	      << " p: "  << &db
	      << std::endl;

  if (db.size()>0) {

    size_t wid = qs.size() * 12;

    std::cout << std::string(70, '-') << std::endl
	      << "VelCrs structure: "  << this << " caller: " << caller 
	      << std::endl << std::string(70, '-') << std::endl;

    for (auto k : db) {
      sKeyPair p = k.first;
      std::ostringstream sout;
      sout << "<" << p.first.first  << "," << p.first.second
	   << "|" << p.second.first << "," << p.second.second << ">";

      std::vector<double> p1, p2, p3;
      for (auto P : qs) {
	p1.push_back(k.second[P][0]());
	p2.push_back(k.second[P][1]());
	p3.push_back(k.second[P][2]());
      }

      std::cout << std::setw(14) << sout.str()
		<< std::setw(10) << k.second[u3][0].count()
		<< "  " << vbkts(p1) << "  product" << std::endl 
		<< std::setw(26) << "" << std::string(wid, '-') << std::endl
		<< std::setw(24) << ""
		<< "  " << vbkts(p2) << "  cross"   << std::endl
		<< std::setw(26) << "" << std::string(wid, '-') << std::endl
		<< std::setw(24) << ""
		<< "  " << vbkts(p3) << "  ratio"   << std::endl
		<< std::setw(26) << "" << std::string(wid, '-') << std::endl;
    }
    std::cout << std::string(70, '-') << std::endl;
  } else {
    std::cout << std::string(70, '-') << std::endl
	      << "VelCrs no points: " << this << std::endl
	      << std::string(70, '-') << std::endl;
  }

  inTest = false;
}

NTCitem::vcMap NTCitem::VelCrsAvg(double _quant)
{
  vcMap ret;
  
  // Check for quant
  //
  const double slop = 1.0e-3;
  double qmax = 1.0e20, quant = 0.5;
  for (auto v : qs) {
    double qtst = fabs(v - _quant);
    if (qmax > qtst) {
      quant = v;
      qmax  = qtst;
    }
  }

  bool ok = (qmax < slop);

  // Deal with round off issues on resume
  //
  for (auto v : db) {
    if (ok)
      ret[v.first] = vcTup(v.second[quant][0](), 
			   v.second[quant][1](), 
			   v.second[quant][2]());
    else
      ret[v.first] = vcTup(NTCitem::VelCrsMin, 0.0, 0.0);
  }

  return ret;
}

NTCitem::vcTup NTCitem::VelCrsAvg(sKeyPair indx, double _quant)
{
  // Check for quant
  //
  const double slop = 1.0e-3;
  double qmax = 1.0e20, quant = 0.5;
  for (auto v : qs) {
    double qtst = fabs(v - _quant);
    if (qmax > qtst) {
      quant = v;
      qmax  = qtst;
    }
  }

  bool ok = (qmax < slop);


  // Get stanza in db
  qpMap::iterator it = db.find(indx);

  // Default value
  if (it == db.end()) return vcTup(VelCrsDef, 0, 0);

  // Debug reporting
  //
  if (DEBUG_V and !inTest) {
    report++;			// Make sure the first one is skipped
    if (report % skprpt == 0 and report < maxrpt*skprpt)  VelCrsTest();
  }

  // Return quantile
  if (ok)
    return vcTup(it->second[quant][0](), 
		 it->second[quant][1](), 
		 it->second[quant][2]());
  else
    return vcTup(NTCitem::VelCrsMin, 0.0, 0.0);
}

void NTCitem::VelCrsAdd(sKeyPair indx, const vcTup& val)
{
  // value below threshold
  //
  if (std::get<0>(val)<=VelCrsMin) return;

  // Check for initialization of Quantile
  //
  if (db.find(indx) == db.end()) db[indx].initialize(qs);

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
