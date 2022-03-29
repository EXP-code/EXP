#ifndef _TWOD_COEFS_H
#define _TWOD_COEFS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>

#include <Eigen/Eigen>

class TwoDCoefs
{
public:

  using DataPtr  = std::shared_ptr<Eigen::MatrixXd>;

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
    std::string id;
    double time, scale, rmax;
    int nmax, Lmax;
    Eigen::MatrixXd data;
    
    bool read(std::istream& in);
  };

  using CoefPtr = std::shared_ptr<Coefs>;

public:
  int Lmax, nmax, ntimes;
  
  std::map<double, Eigen::MatrixXd> coefs;
  
  std::vector<double> times;

  TwoDCoefs(const std::string& file, unsigned stride=1);

  DataPtr operator() (const double time);
  DataPtr interpolate(const double time);
};

#endif
