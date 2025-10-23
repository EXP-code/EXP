#ifndef _MixtureBasis_H
#define _MixtureBasis_H

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

#include "TwoCenter.H"

class MixtureBasis
{

  typedef double (TwoCenter::*mixFunc)(double *p);

private:

  vector<double>* ctr;
  TwoCenter *p;
  mixFunc f;

  // For debugging
  string id;

public:

  MixtureBasis(TwoCenter& instance, vector<double> *c,
	       string ID, mixFunc func) : 
    p(&instance), ctr(c), id(ID), f(func) {}

  void getCenter(vector<double>& c) { c = *ctr; }
  double Mixture(double* pos) { return CALL_MEMBER_FN(*p, f)(pos); }
  
};

#endif
