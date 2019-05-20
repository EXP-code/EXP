#ifndef _Elastic_H
#define _Elastic_H

#include <map>
#include <cmath>
#include <vector>

class Elastic 
{
public:

  //@{
  //! Database structures
  typedef std::map<double, double> InterpPair;
  typedef std::map<unsigned short, InterpPair> xSection;
  //@}

  static bool extrapolate;
  static bool pin;

protected:

  //! The cross-section database
  xSection atomicdata, iondata;

  //! Internal interpolation
  double interpolate(const std::map<double, double>& data, double x);

public:
  
  //! Charged particle type
  enum CPtype { electron, proton };

  //! Constructor
  Elastic();

  /** Get the cross section for atomic number Z and energy in eV.  The
      cross section is units of nm^2 for consistency with CHIANTI.  

      If you ask for data beyond the grid for a given Z, the closest
      value is used.

      Similarly, if you ask for a set that does not exist, data for
      the closest Z is used.
  */
  double operator()(unsigned short Z, double x, CPtype type=electron)
  {
    // Bohr cross section (pi*a_0^2) in nm
    const double b_cross = 0.00879735542978;
    Z = std::max<unsigned short>(Z, 1);
    if (type == electron) {
      if (atomicdata.find(Z)==atomicdata.end()) return 0.0;
      xSection::const_iterator i = atomicdata.lower_bound(Z);
      if (i==atomicdata.end()) --i;
      return b_cross * interpolate(i->second, x);
    } else {
      if (iondata.find(Z)==iondata.end()) return 0.0;
      // Make sure energy is on grid
      x = log10(x);
      if (x > iondata[Z].rbegin()->first) return 0.0;
      // Ok . . . do the interpolation
      xSection::const_iterator i = iondata.lower_bound(Z);
      if (i==iondata.end()) --i;
      return b_cross * pow(10.0, interpolate(i->second, x));
    }
  }

};

class Geometric
{
public:
  typedef std::map<unsigned short, double> radData;

private:
  radData radii;

public:
  
  //! Constructor
  Geometric();

  /** Get the effective radius for the geometric cross section for
      atomic number Z based on computed atomic radii from Clementi et
      al. (see Elastic.cc). The units are in nm for consistency with
      CHIANTI.

      If you ask for data not in the data base, the closest Z value is
      used.
  */
  double operator()(unsigned short Z)
  {
    radData::const_iterator i = radii.lower_bound(Z);

    // Convert radii in picometers to nanometers
    double radius = (i->second) * 1.0e-3;

    return radius;
  }

};

#endif
