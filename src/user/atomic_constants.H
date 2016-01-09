#ifndef _atomic_constants_H
#define _atomic_constants_H

// Proton mass (g)
const double mp        = 1.67262158e-24;

// Electron mass (g)
const double me        =  9.10938291e-28;

// Speed of light (cm/s)
const double light     =  2.99792458e+10;

// electron volt in (cgs)
const double eV        =  1.60217653e-12;

// Boltzmann constant (cgs)
const double boltz     = 1.3806504e-16;

// Boltzmann constant (eV)
const double boltzEv   = 8.6173324e-5;

// Planck's constant (cgs)
const double planck    = 6.626075540e-27;

// Electron charge (cgs)
const double esu       = 4.80320427e-10;

// Atomic mass unit in grams
const double amu       = 1.660539e-24;

// Electron rest mass in MeV
const double mec2      = 0.51099891013;

// Parsec (cm)
const double pc        = 3.08567758e18;

// Bohr radius (cm)
const double a0        = 0.052917721092e-7;

// Year (sec)
const double year      = 365.242*24*3600;

// Solar mass (g)
const double msun      = 1.9891e33;

#include <tuple>
#include <map>

#include <boost/shared_ptr.hpp>

/** Store and retrieve periodic table info.  This can be extended to
    include additional data fields, etc.  Would be nice to add atomic
    weights and possibly family info in the future.
*/
class AtomicElement : public std::tuple<std::string, std::string, unsigned, double>
{
public:
  
  //! Null constructor: returns initialized but empty record
  AtomicElement() 
  {
    std::get<2>(*this) = 0;
    std::get<3>(*this) = 1.0;
  }

  //! Constructor
  AtomicElement(const std::string& name, 
                const std::string& abbrev, 
                const unsigned z)
  {
    std::get<0>(*this) = name;
    std::get<1>(*this) = abbrev;
    std::get<2>(*this) = z;
    std::get<3>(*this) = 1.0;
  }

  //! Copy constructor
  AtomicElement(const AtomicElement & p)
  {
    std::get<0>(*this) = std::get<0>(p);
    std::get<1>(*this) = std::get<1>(p);
    std::get<2>(*this) = std::get<2>(p);
    std::get<3>(*this) = std::get<3>(p);
  }

  //! Element name
  const std::string&   name()   { return std::get<0>(*this); }

  //! Element abbreviation
  const std::string&   abbrev() { return std::get<1>(*this); }

  //! Atomic number
  const unsigned       number() { return std::get<2>(*this); }

  //! Data (used for cross-section scale factor by CollideIon)
  const double         scale()  { return std::get<3>(*this); }

  //! Set double value (currently used by CollideIon as above)
  void set(double s) { std::get<3>(*this) = s; }
};

//! A full periodic table registry
class PeriodicTable
{
private:

  typedef boost::shared_ptr<AtomicElement> AEptr;

  //! The atomic data
  std::map<unsigned,    AEptr> dataZ;
  std::map<std::string, AEptr> dataS;

  //! Add data to database (only used for initialization)
  void add(const std::string& n, const std::string& a, const unsigned z) 
  {
    AEptr p(new AtomicElement(n, a, z));
    dataZ[z] = p;
    dataS[n] = p;
  }

public:

  //! Initialize data
  PeriodicTable();

  //! Retrieve the data tuple based on atomic number
  AtomicElement* operator[](unsigned z) { 
    std::map<unsigned, AEptr>::iterator it = dataZ.find(z);
    if (it == dataZ.end()) return 0;
    return it->second.get();
  }

  //! Retrieve the data tuple based on abbrev
  AtomicElement* operator[](std::string s) { 
    std::map<std::string, AEptr>::iterator it = dataS.find(s);
    if (it == dataS.end()) return 0;
    return it->second.get();
  }
};


#endif
