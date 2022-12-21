#ifndef _atomic_constants_H
#define _atomic_constants_H

// Fine structure constant
constexpr double afs       = 0.0072973525693;

// Proton mass (g)
constexpr double mp        = 1.67262158e-24;

// Electron mass (g)
constexpr double me        =  9.10938291e-28;

// Speed of light (cm/s)
constexpr double light     =  2.99792458e+10;

// electron volt in (cgs)
constexpr double eV        =  1.60217653e-12;

// Boltzmann constant (cgs)
constexpr double boltz     = 1.3806504e-16;

// Boltzmann constant (eV)
constexpr double boltzEv   = 8.6173324e-5;

// Planck's constant (cgs)
constexpr double planck    = 6.626075540e-27;

// Planck's constant h-bar (eV-s)
constexpr double hbar      = 6.58211951440e-16;

// Electron charge (cgs)
constexpr double esu       = 4.80320427e-10;

// Atomic mass unit in grams
constexpr double amu       = 1.660539e-24;

// Electron rest mass in MeV
constexpr double mec2      = 0.51099891013;

// Parsec (cm)
constexpr double pc        = 3.08567758e18;

// Bohr radius (cm)
constexpr double a0        = 0.052917721092e-7;

// Year (sec)
constexpr double year      = 365.242*24*3600;

// Solar mass (g)
constexpr double msun      = 1.9891e33;

#include <algorithm>
#include <iostream>
#include <cctype>
#include <string>
#include <tuple>
#include <map>

#include <memory>

//! Case insenstive comparison
struct StrMapComparator
{
  bool operator()( const std::string& a, const std::string& b ) const 
  {
    return std::lexicographical_compare
      (
       a.begin(), a.end(),
       b.begin(), b.end(),
       [](auto ch1, auto ch2) {
	 return std::tolower(ch1) < std::tolower(ch2);
       });
  }
};

/** Store and retrieve periodic table info.  This can be extended to
    include additional data fields, etc.  Would be nice to add atomic
    weights and possibly family info in the future.
*/
class AtomicElement : public std::tuple<std::string, std::string, 
					unsigned, double, double>
{
public:
  
  //! Null constructor: returns initialized but empty record
  AtomicElement() 
  {
    std::get<2>(*this) = 0;
    std::get<3>(*this) = 0.0;
    std::get<4>(*this) = 1.0;
  }

  //! Constructor
  AtomicElement(const std::string& name, 
                const std::string& abbrev, 
                const unsigned z, const double w)
  {
    std::get<0>(*this) = name;
    std::get<1>(*this) = abbrev;
    std::get<2>(*this) = z;
    std::get<3>(*this) = w;
    std::get<4>(*this) = 1.0;
  }

  //! Copy constructor
  AtomicElement(const AtomicElement & p)
  {
    std::get<0>(*this) = std::get<0>(p);
    std::get<1>(*this) = std::get<1>(p);
    std::get<2>(*this) = std::get<2>(p);
    std::get<3>(*this) = std::get<3>(p);
    std::get<4>(*this) = std::get<4>(p);
  }

  //! Element name
  const std::string&   name()   { return std::get<0>(*this); }

  //! Element abbreviation
  const std::string&   abbrev() { return std::get<1>(*this); }

  //! Atomic number
  unsigned             number() { return std::get<2>(*this); }

  //! Atomic weight
  double               weight() { return std::get<3>(*this); }

  //! Data (used for cross-section scale factor by CollideIon)
  double               scale()  { return std::get<4>(*this); }

  //! Set double value (currently used by CollideIon as above)
  void set(double s) { std::get<4>(*this) = s; }
};

//! A full periodic table registry
class PeriodicTable
{
private:

  typedef std::shared_ptr<AtomicElement> AEptr;

  //! The atomic data
  std::map<unsigned,    AEptr> dataZ;
  std::map<std::string, AEptr, StrMapComparator> dataS;

  //! Add data to database (only used for initialization)
  void add(const std::string& n, const std::string& a, const unsigned z,
	   const double w) 
  {
    AEptr p(new AtomicElement(n, a, z, w));

    dataS[a] = p;
    dataZ[z] = p;
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

  //! Print data to an output stream, defaulting to stdout
  void print(std::ostream& out=std::cout);

};


#endif
