#ifndef _CYLINDRICAL_COEFS_H
#define _CYLINDRICAL_COEFS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
#include <tuple>
#include <map>

#include "coef.H"

class CylindricalCoefs
{
public:
  
  using Dvector  = std::vector<std::vector<double>>;
  using Dpair    = std::tuple<Dvector, Dvector>;
  using DataPtr  = std::shared_ptr<Dpair>;

private:

  int ndigits;

  double to_ndigits(double x)
  { 
    double scale = std::pow(10.0, ndigits);
    return std::round(x * scale) / scale;
  }

  DataPtr ret;

public:

  struct Coefs
  {
    double time;
    int nmax, mmax;
    std::vector< std::vector<double> > cos_c, sin_c;
    
    bool read(std::istream& in);
  };

  using CoefPtr= std::shared_ptr<Coefs>;

public:
  int mmax, nmax, ntimes;
  
  std::map<double, CoefPtr> data;
  
  std::vector<double> times;

  CylindricalCoefs(const std::string& file, int ndigits=10, unsigned stride=1);

  DataPtr operator() (const double time);
  DataPtr interpolate(const double time);
};

#endif
