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
public:

  using Dvector = std::vector<std::vector<double>>;
  using DataPtr  = std::shared_ptr<Dvector>;

private:

  int ndigits;

  double to_ndigits(double x)
  { 
    double scale = std::pow(10.0, ndigits);
    return std::round(x * scale) / scale;
  }

  DataPtr ret;
  DataPtr zero_ret();

public:
  using LMkey = std::tuple<unsigned, unsigned>;

  struct Coefs
  {
    double time, rscl;
    int nmax, lmax;
    std::vector<std::vector<double>> data;
    
    bool read(std::istream& in);
  };

  using CoefPtr = std::shared_ptr<Coefs>;

public:
  int lmax, nmax, ntimes;
  
  std::map<double, Dvector> coefs;
  
  std::vector<double> times;

  SphericalCoefs(const std::string& file, unsigned stride=1);

  DataPtr operator() (const double time);
  DataPtr interpolate(const double time);
};

#endif
