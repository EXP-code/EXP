// This may look like C code, but it is really -*- C++ -*-

#ifndef _Gamma_h
#define _Gamma_h

#include <Random.h>

/**
   Generate Gamma variates
   Adapted from TOMS 599 by MDW
*/
class Gamma : public Random 
{
 protected:
  double pAlpha;
  char haveCachedNormal;
  double cachedNormal;
  double getNormal();
    
 public:
  Gamma(double alpha, RNG *gen);

  double alpha();
  double alpha(double x);

  virtual double operator()();
};

    
inline Gamma::Gamma(double alpha, RNG *gen) : Random(gen)
{
  pAlpha = alpha;
  haveCachedNormal = 0;
}

inline double Gamma::alpha() { return pAlpha; }

inline double Gamma::alpha(double x) {
  double tmp = pAlpha;
  pAlpha = x;
  return tmp;
}

#endif
