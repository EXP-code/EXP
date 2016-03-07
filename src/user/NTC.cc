#include <iostream>
#include <iomanip>
#include <sstream>

#include "Particle.H"
#include "pCell.H"
#include "NTC.H"

using namespace NTC;

// For verbose debugging
static const bool DEBUG_V = false;

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

// Default interaction key for singleton
NTC::T NTCitem::single {0, 0, 0};

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::tuple<T, T, T>& t)
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

   // Return value
   return ret;
 }

 void NTCitem::Add(sKeyPair indx, const T& intr, double val)
 {
   oddBall(indx, intr);

   // value below threshold
   //
   if (val<=Min) return;

   // Check for initialization of Quantile
   //
   if (db.find(indx) == db.end()) {
     db[indx][intr].histogram(Nequal);
     for (auto v : qs) db[indx][intr].newQ(v);
   }

   if (0) {
     std::ostringstream sout; sout << "<" << val << ">";

     std::cout << "Adding " << std::setw(16) << sout.str() << " to " << indx 
	       << " [" << db[indx][intr].datums() 
	       << "/"  << db[indx][intr].target() << "]";
     if (db[indx][intr].full()) 
       std::cout << " P(0.5) =" 
		 << std::setw(16) << db[indx][intr](0.5)
		 << " P(0.95) =" 
		 << std::setw(16) << db[indx][intr](0.95)
		 << " xmax =" 
		 << std::setw(16) << db[indx][intr].xmax();
     std::cout << std::endl;
   }

   // Add new element
   //
   db[indx][intr].add(val);
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
