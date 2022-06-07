#include <filesystem>
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

// Count live instances
unsigned NTCitem::instance  = 0;

// Default null
speciesKey const NTC::nullKey {speciesKey(0, 0)};

// Default interaction key for electron
static unsigned short elecVal = std::numeric_limits<unsigned short>::max();

// Default electron
speciesKey const NTC::electron {speciesKey(elecVal, elecVal)};

// Default electron
speciesKey const NTC::proton {speciesKey(1, 2)};

// Default interaction key for singleton
T const NTC::single {0, nullKey, nullKey};

std::ostream& operator<< (std::ostream& out, const speciesKey& k)
{
  std::ostringstream sout;
  sout << " [" << k.first << ", " << k.second << "]";
  out << sout.str();
  return out;
}

std::ostream& operator<< (std::ostream& out, const NTC::T& e)
{
  std::ostringstream sout;
  sout << std::get<0>(e)
       << std::get<1>(e)
       << std::get<2>(e);
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

/*
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
*/

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
  if (it == db.end()) throw Error();
  
  // Get next stanza
  uqMap::iterator jt = it->second.find(intr);
  
  // Default value
  if (jt == it->second.end()) throw Error();
  
  
  // Default value
  if (! jt->second.full() ) throw Error();
  
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
    if (false and data.size()>1000)
      std::cout << "NTCdb [" << std::setw(2) << myid << "]: "
		<< std::hex << k << " not found, tnow="
		<< std::setw(10) << tnow << " live=" << std::dec
		<< std::setw(9) << std::left << NTCitem::instance
		<< " size=" << std::setw(9) << data.size() << std::endl;

    // Look for an initialization parent
    //
    key_type check = k >> 3;
    it = data.find(check);
    while (it == data.end() and check != 0ul) {
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


//! Save contents to a file
void NTCdb::save_myself(const std::string& filename)
{
  namespace fs = std::filesystem;
      
  // Synchronize
  sync();
      
  // Only the root process makes an archive
  if (myid==0) {
    std::ofstream ofs;
    bool exist = fs::exists(filename);
    std::string tfile(".tmp." + filename);
    
    if (exist) {
      ofs.open(tfile.c_str());
    } else {
      ofs.open(filename.c_str());
    }
    
#ifdef BINARY_ARCHIVE
    boost::archive::binary_oarchive oa(ofs);
#else
    boost::archive::xml_oarchive oa(ofs);
#endif
    oa << BOOST_SERIALIZATION_NVP(pdata);
    oa << BOOST_SERIALIZATION_NVP(curr_time);
    oa << BOOST_SERIALIZATION_NVP(next_time);
    oa << BOOST_SERIALIZATION_NVP(debug_count);
    
    if (exist) {
      fs::rename(tfile, filename);
    }
    
    if (chatty and NTCitem::live()) {
      std::cout << std::endl << std::string(70, '-') << std::endl
		<< "-- Took NTC checkpoint at time <" << curr_time << ">" << std::endl
		<< "-- Next NTC checkpoint at time <" << next_time << ">" << std::endl
		<< std::string(70, '-') << std::endl;
    }
  }
      
  // Clean up node-aggregated db
  pdata.clear();

  // Debugging
  if (myid==0 and NTCitem::live()) {
    std::cout << std::string(40, '-') << std::endl
	      << "--- QuantileBag instances: " << QuantileBag::live()
	      << std::endl
	      << "--- NTCitem instances:  " << NTCitem::live()
	      << std::endl
	      << std::string(40, '-') << std::endl;
  }
  dump();
}


void NTCdb::restore_myself(const std::string& filename)
{
  // open the archive; all processes
  std::ifstream ifs(filename.c_str());
#ifdef BINARY_ARCHIVE
  boost::archive::binary_iarchive ia(ifs);
#else
  boost::archive::xml_iarchive ia(ifs);
#endif
  // restore the data from the archive
  ia >> BOOST_SERIALIZATION_NVP(data);
  ia >> BOOST_SERIALIZATION_NVP(curr_time);
  ia >> BOOST_SERIALIZATION_NVP(next_time);
  ia >> BOOST_SERIALIZATION_NVP(debug_count);
  
  if (myid==0 and chatty) {
    std::cout << std::endl << std::string(70, '-') << std::endl
	      << "-- Restoring NTC checkpoint from time <" << curr_time << ">" << std::endl
	      << "-- Current simulation time is <" << tnow << ">" << std::endl
	      << "-- Next NTC checkpoint scheduled for <" << next_time << ">" << std::endl
	      << std::string(70, '-') << std::endl;
  }
  
  // For debugging . . . 
  dump();
}

// Constructor
NTCdb::NTCdb()
{
  namespace fs = std::filesystem;

  // Attempt to restore from save set
  try {
    std::string filename = ".ntcdb." + runtag;
    if ( fs::exists(filename) ) {
      restore_myself(filename);
    }
  }
  catch (fs::filesystem_error & error) {
    if (myid==0) {
      std::cout << "Failure to open archive <.ntcdb." << runtag
		<< "> with error: " << error.what() << std::endl
		<< "We will make a new archive.  Continuing . . . "
		<< std::endl;
    }
  }	
  //
  curr_time = tnow;
  next_time = curr_time + dtime*(static_cast<double>(intvl)-eps);
}
