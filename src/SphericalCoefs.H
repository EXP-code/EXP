#ifndef _SPHERICAL_COEFS_H
#define _SPHERICAL_COEFS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>

#include "coef.H"

class SphericalCoefs
{
private:

  int ndigits;

  double to_ndigits(double x)
  { 
    double scale = std::pow(10.0, ndigits);
    return std::round(x * scale) / scale;
  }

public:
  using LMkey = std::tuple<unsigned, unsigned>;

  struct Coefs
  {
    double time;
    int nmax, lmax;
    std::map<LMkey, std::vector<double> > cos_c, sin_c;
    
    bool read(std::istream& in);
  };

  using CoefPtr = std::shared_ptr<Coefs>;

public:
  int lmax, nmax, ntimes;
  
  using Dvector  = std::vector<double>;
  using D2vector = std::vector<Dvector>;

  std::map<double, CoefPtr>  data;
  std::map<double, D2vector> coefs;
  
  std::vector<double> times;

  SphericalCoefs(const std::string& file, unsigned stride=1);

  D2vector interpolate(const double time);
};

#endif
