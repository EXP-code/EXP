#include <iostream>
#include <iomanip>
#include <sstream>

#include "Particle.H"
#include "pCell.H"
#include "NTC.H"

using namespace NTC;

// For verbose debugging
static const bool DEBUG_V   = false;

// Skip this number of calls between reports
unsigned NTCitem::skprpt    = 5;

// Maximum number of reports for each cell item
unsigned NTCitem::maxrpt    = 5;

// Minimum CrossSection x Velocity value
double   NTCitem::VelCrsMin = 1.0e-24;

// Minimum CrossSection x Velocity value
double   NTCitem::VelCrsDef = 10.0;

void NTCitem::VelCrsTest()
{
  inTest = true;		// To prevent recursive calling in VelCrsAvg

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

      vcTup q1(k.second[u1].get());
      vcTup q2(k.second[u2].get());
      vcTup q3(k.second[u3].get());

      vcTup p1( std::get<0>(q1),  std::get<0>(q2), std::get<0>(q3) );
      vcTup p2( std::get<1>(q1),  std::get<1>(q2), std::get<1>(q3) );
      vcTup p3( std::get<2>(q1),  std::get<2>(q2), std::get<2>(q3) );

      std::cout << std::setw(14) << sout.str()
		<< std::setw(10) << k.second[u3].count()
		<< "  " << p1 << "  product" << std::endl 
		<< std::setw(26) << "" << std::string(42, '-') << std::endl
		<< std::setw(24) << ""
		<< "  " << p2 << "  cross"   << std::endl
		<< std::setw(26) << "" << std::string(42, '-') << std::endl
		<< std::setw(24) << ""
		<< "  " << p3 << "  ratio"   << std::endl
		<< std::setw(26) << "" << std::string(42, '-') << std::endl;
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
      ret[v.first] = v.second[quant].get();
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
    return it->second[quant].get();
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
  if (db.find(indx) == db.end()) db[indx].initialize(ic);

  // Add new element
  //
  db[indx].add(val);
}


NTCitem NTCdb::operator[](const key_type& k)
{
  NTCdata::iterator it = data.find(k);
  NTCitem ret;

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
      data[k] = ret;
    } else {
      ret = it->second;
    }

  } else {
    ret = it->second;
  }

  if (DEBUG_V) ret.setKey(k);

  return ret;
}
