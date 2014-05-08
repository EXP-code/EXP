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

  //! Get the cross section for atomic number Z and energy x in ergs.
  //! The cross section is in cgs units (cm^2 * 1.0e14)
  double operator()(unsigned short Z, double x)
  {
    Z = std::max<unsigned short>(Z, 1);
    xSection::const_iterator i = atomicdata.lower_bound(Z);
    if (i==atomicdata.end()) --i;
    return interpolate(i->second, x);
  }

};

#endif
