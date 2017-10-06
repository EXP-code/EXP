#include <iostream>
#include <iomanip>
#include <sstream>

#include "global.H"
#include "Particle.H"
#include "pCell.H"
#include "NTC.H"

using namespace NTC;

// For verbose debugging
static const bool DEBUG_V = false;

// Copy data to ascii file
bool NTCdb::debug_output = false;

// Ascii version number
unsigned NTCdb::debug_count = 1;

// Maximum age for uncalled cell entries in number of time steps
unsigned NTCdb::maxAge = 20;

// Chatty output for debugging
bool NTCdb::chatty = false;

// Save interval
unsigned NTCdb::intvl = 10;

// Skip this number of calls between reports
double   NTCitem::def_time  = -1.0e32;

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

// Label database
std::map<Interact::pType, std::string>
Interact::pLabs = {
  {Interact::simple,   "Simple"},
  {Interact::neutral,  "Neutral"},
  {Interact::ion,      "Ion"},
  {Interact::electron, "Electron"}
};

// Default null pElem
constexpr Interact::pElem const Interact::pdef {Interact::simple, speciesKey(0, 0)};

// Default electron pElem
constexpr Interact::pElem const Interact::edef {Interact::electron, speciesKey(0, 0)};

// Default interaction key for singleton
Interact::T Interact::single(0, Interact::pdef, Interact::pdef);

// Default interaction key for singleton
constexpr static unsigned short singleVal = std::numeric_limits<unsigned short>::max();
Interact::T NTCitem::single {singleVal, Interact::pdef, Interact::pdef};

std::ostream& operator<< (std::ostream& out, const Interact::pElem& e)
{
  std::ostringstream sout;
  sout << Interact::label(e)
       << " [" << e.second.first << ", " << e.second.second << "]";
  out << sout.str();
  return out;
}

template<typename T, typename U, typename V>
std::ostream& operator<< (std::ostream& out, const std::tuple<T, U, V>& t)
{
  std::ostringstream sout;
  sout << '(' << std::get<0>(t)
       << ',' << std::get<1>(t)
       << ',' << std::get<2>(t)
       << ')';
  out << sout.str();
  return out;
}

std::ostream& operator<< (std::ostream& out, const sKeyPair& t)
{
  std::ostringstream sout;
  sout << "[(" << t.first.first
       << ',' << t.first.second
       << ")(" << t.second.first
       << ',' << t.second.second
       << ")]";
  out << sout.str();
  return out;
}

void NTCdb::dump()
{
  if (myid>0 or not debug_output) return;
  
  std::ostringstream sout;
  sout << "ntcdb_debug." << debug_count++;

  std::ofstream out(sout.str());
  if (out) {
    out << std::string(72, '-') << std::endl
	<< std::string( 5, '=') << " Time=" << tnow << std::endl;

    for (auto v : data) {
      out << std::string(72, '-') << std::endl
	  << std::string(10, '=') << " Key " << v.first << std::endl;
      v.second.dump(out);
    }
  }
}  

void NTCitem::dump(std::ostream& out)
{
  for (auto v : db) {
    out << std::string(72, '-') << std::endl
	<< std::string(15, '=') << " Key pair " << v.first << std::endl;
    v.second.dump(out);
  }
}


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
    
    for (auto v : db) {
      sKeyPair p = v.first;
      std::ostringstream sout;
      sout << "<" << p.first.first  << "," << p.first.second
	   << "|" << p.second.first << "," << p.second.second << ">";
      
      size_t prc = std::cout.precision(2);
      
      for (auto k : v.second) {
	
	std::cout << std::setw(14) << sout.str()
		  << std::setw(10) << k.first
		  << std::setw(10) << k.second.count();
	
	for (size_t i=0; i<number; i++) {
	  double P = (0.5 + i)/number;
	  std::cout << std::setw(10) << k.second(P);
	}
	std::cout << std::endl;
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
  for (auto v : db) {
    for (auto u : v.second)
      ret[v.first][u.first] = u.second.inverse(x);
  }
  return ret;
}

double NTCitem::Prob(sKeyPair indx, const T& intr, double x)
{
  oddBall(indx, intr);
  
  // Get stanza in db
  qpMap::iterator jt = db.find(indx);
  
  // Default value
  if (jt == db.end()) return 0.5;
  
  // Get next stanza 
  uqMap::iterator it = jt->second.find(intr);
  
  // Default value
  if (it == jt->second.end()) return 0.5;
  
  // Debug reporting
  //
  if (DEBUG_V and !inTest) {
    report++;			// Make sure the first one is skipped
    if (report % skprpt == 0 and report < maxrpt*skprpt)  Test();
  }
  
  // Return quantile
  return it->second.inverse(x);
}

void NTCitem::debug()
{
  std::cout << std::string(94, '-') << std::endl;
  std::cout << std::left << std::setw(30) << "Species pair"
	    << std::left << std::setw( 6) << "Inter"
	    << std::left << std::setw(10) << "Count"
	    << std::left << std::setw(18) << "Min value"
	    << std::left << std::setw(18) << "Max value"
	    << std::left << std::setw(18) << "Quantile"
	    << std::endl;
  
  for (auto u : db) {
    for (auto v : u.second) {
      std::cout << std::left << std::setw(30) << u.first 
		<< std::left << std::setw(10) << v.first
		<< std::left << std::setw(10) << v.second.count() 
		<< std::left << std::setw(18) << v.second.xmin()
		<< std::left << std::setw(18) << v.second.xmax()
		<< std::left << std::setw(18) << v.second(0.95)
		<< std::endl;
      v.second.dump(std::cout);
    }    
  }
  std::cout << std::string(94, '-') << std::endl;
}

double NTCitem::CrsVel(sKeyPair indx, const T& intr, double p)
{
  oddBall(indx, intr);
  
  // Get stanza in db
  qpMap::iterator it = db.find(indx);
  
  // Default value
  if (it == db.end()) return Def;
  
  // Get next stanza
  uqMap::iterator jt = it->second.find(intr);
  
  // Default value
  if (jt == it->second.end()) return Def;
  
  
  // Default value
  if (! jt->second.full() ) return Def;
  
  if (0) debug();
  
  // Set time stamp
  //
  ts = tnow;

  // Return value
  return jt->second(p);
}

NTCitem::vcMap NTCitem::CrsVel(double p)
{
  vcMap ret;
  
  for (auto v : db) {
    for (auto u : v.second)
      ret[v.first][u.first] = u.second(p);
  }
  
  // Set time stamp
  //
  ts = tnow;

  // Return value
  return ret;
}

void NTCitem::Add(sKeyPair indx, const T& intr, double val)
{
  oddBall(indx, intr);
  
  // value below threshold
  //
  if (val<=Min) return;
  
  // Add new element
  //
  db[indx][intr].add(val);

  // Set time stamp
  //
  ts = tnow;
}

NTCitem& NTCdb::operator[](const key_type& k)
{
  NTCdata::iterator it = data.find(k);
  
  if (it == data.end()) {
    //  +-- Debugging
    //  |
    //  v
    if (true and data.size()>1000)
      std::cout << "NTCdb [" << std::setw(2) << myid << "]: "
		<< std::hex << k << " not found, tnow="
		<< std::setw(10) << tnow << " live=" << std::dec
		<< std::setw(9) << std::left << NTCitem::instance
		<< " size=" << std::setw(9) << data.size() << std::endl;

    // Look for an initialization parent
    //
    key_type check = k >> 3;
    it = data.find(check);
    while (it == data.end() and check != 0) {
      check = check >> 3;
      it = data.find(check);
    }
    
    // If none, create an empty item.  Otherwise, initialize from the
    // parent
    //
    if (it == data.end()) data[k] = NTCitem();
    else                  data[k] = NTCitem(it->second);
    it = data.find(k);
  }
  
  if (DEBUG_V) it->second.setKey(k);
  
  // Set time stamp
  //
  it->second.ts = tnow;

  return it->second;
}
