#ifndef _MixtureSL_H
#define _MixtureSL_H

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

#include <SphereTwoCenter.H>

class MixtureSL
{

  typedef double (SphereTwoCenter::*mixFunc)(double *p);

private:

  vector<double>* ctr;
  SphereTwoCenter *p;
  mixFunc f;

  // For debugging
  string id;

public:

  MixtureSL(SphereTwoCenter& instance, vector<double> *c,
	    string ID, mixFunc func) : 
    p(&instance), ctr(c), id(ID), f(func) {}

  void getCenter(vector<double>& c) { c = *ctr; }
  double Mixture(double* pos) { return CALL_MEMBER_FN(*p, f)(pos); }
  
};

#endif
