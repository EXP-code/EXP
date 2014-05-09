#ifndef _Elastic_H
#define _Elastic_H

#include <map>
#include <vector>

class Elastic 
{
public:

  //@{
  //! Database structures
  typedef std::map<double, double> InterpPair;
  typedef std::map<unsigned short, InterpPair> xSection;
  //@}

private:

  //! The cross-section database
  xSection atomicdata;

  //! Internal interpolation
  double interpolate(const std::map<double, double>& data, double x);

public:
  
  //! Constructor
  Elastic();

  /** Get the cross section for atomic number Z and energy in eV.  The
      cross section is units of nm^2 for consistency with CHIANTI.  

      If you ask for data beyond the grid for a given Z, the closest
      value is used.

      Similarly, if you ask for a set that does not exist, data for
      the closest Z is used.
  */
  double operator()(unsigned short Z, double x)
  {
    // Bohr cross section (pi*a_0^2) in nm
    const double b_cross = 0.00879735542978;
    Z = std::max<unsigned short>(Z, 1);
    xSection::const_iterator i = atomicdata.lower_bound(Z);
    if (i==atomicdata.end()) --i;
    return b_cross * interpolate(i->second, x);
  }

};

#endif
