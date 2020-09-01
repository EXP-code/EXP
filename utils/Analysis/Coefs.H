#ifndef _COEFS_H
#define _COEFS_H

#include <boost/make_unique.hpp>

#include <yaml-cpp/yaml.h>	// YAML support

#include "../../src/coef.H"

struct Coefs
{
  double time;
  int nmax, mmax;
  std::vector< std::vector<double> > cos_c, sin_c;

  bool read(std::istream& in, bool verbose=false);
};

typedef std::shared_ptr<Coefs> CoefPtr;

struct SphCoefs
{
  SphCoefHeader header;
  std::vector< std::vector<double> > coefs;

  bool read(std::istream& in, bool exp_type);
};

typedef std::shared_ptr<SphCoefs> SphCoefsPtr;

#endif
