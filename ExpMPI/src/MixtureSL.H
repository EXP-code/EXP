#ifndef _MixtureSL_H
#define _MixtureSL_H

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

#include <SphereTwoCenter.H>

class MixtureSL
{

  typedef double (SphereTwoCenter::*mixFunc)(double *p);

private:

  vector<double> ctr;
  SphereTwoCenter *p;
  mixFunc f;

  // For debugging
  string id;

public:

  MixtureSL(SphereTwoCenter& instance, string ID, mixFunc func) : 
    p(&instance), id(ID), f(func) {}

  void setCenter(vector<double>& c) { ctr = c; }
  void getCenter(vector<double>& c) { c = ctr; }
  double Mixture(double* pos) { return CALL_MEMBER_FN(*p, f)(pos); }
  
};

#endif
