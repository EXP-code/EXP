#ifndef _Collide_H
#define _Collide_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <ACG.h>
#include <Uniform.h>

#include "Particle.H"
#include "pHOT.H"
#include "Timer.h"

class Collide;

//! Helper class to pass info to threaded member
struct thrd_pass_arguments
{
  //! Pointer to iterator
  pHOT *tree;
  //! Effective number of particles
  double fn;
  //! Time step
  double tau;
  //! Thread counter id
  int id;
};
  
struct thrd_pass_Collide 
{
  Collide *p;
  //! Arguments
  thrd_pass_arguments arg;
};


class Collide
{
private:

  // Threading stuff

  thrd_pass_Collide *td;
  pthread_t *t;
  Timer forkTime, snglTime, waitTime, joinTime;
  TimeElapsed forkSoFar, snglSoFar, waitSoFar, joinSoFar;

  vector<Timer> collTime;
  vector<int> collCnt;
  vector<TimeElapsed> collSoFar;

  ACG *gen;
  Uniform *unit;

protected:

  int nthrds;

  std::vector<unsigned> numcnt;
  std::vector<unsigned> colcnt;

  std::vector< vector<unsigned> > numcntT;
  std::vector< vector<unsigned> > colcntT;
  std::vector<unsigned> error1T;
  std::vector<unsigned> col1T;
  std::vector<double> KEtotT, KElostT;

  const static unsigned numdiag = 9;
  std::vector<unsigned> tdiag, tdiag1, tdiag0;
  std::vector< vector<unsigned> > tdiagT;

  std::vector< vector<pCell*> > cellist;

  double diam0, diam;
  unsigned coltot;
  unsigned errtot;
  unsigned numtot;

  virtual void initialize_cell(pCell* cell, double vrelmax,
			       double dt, double number, int id=0) {}

  virtual int inelastic(pHOT* tree, Partstruct* p1, Partstruct* p2,
			double *cr, int id=0) { return 0; }

  void debug_list(pHOT& tree);

public:

  static int CNUM;		// Target CNUM collisions if nonzero
				// default: 0

  static bool CBA;		// Use the Consistent Boltzmann Algorithm
				// Alexander, Garcia, Alder 1995, PRL, 74, 5212

  Collide(double diameter, int nth=1);
  ~Collide();

  unsigned total()  { return coltot; }
  unsigned errors() { return errtot; }
  unsigned cellN()  { return numtot; }
  unsigned medianNumber();
  unsigned medianColl();
  unsigned collide(pHOT& cl, double Fn, double tau);
  void *collide_thread(void *arg);
  void collide_thread_fork(pHOT *tree, double Fn, double tau);

  double threadTime() { forkTime.reset(); return forkSoFar()*1.0e-6; }
  double joinedTime() { snglTime.reset(); return snglSoFar()*1.0e-6; }
  double waitngTime() { waitTime.reset(); return waitSoFar()*1.0e-6; }
  double joinngTime() { joinTime.reset(); return joinSoFar()*1.0e-6; }

  void colldeTime() { 
    for (int n=0; n<nthrds; n++) {
      cout << "   " << setw(3) << n << setw(18)
	   << collSoFar[n]()*1.0e-6 << setw(10) << collCnt[n] << endl;
      collTime[n].reset(); 
      collCnt[n] = 0;
    }
  }

  void tsdiag(ostream& out) {
    out << "-----Time step diagnostics---------------" << endl;
    for (unsigned k=0; k<numdiag; k++) {
      out << setw(15) << pow(10.0, -3.5+k)
	   << setw(15) << tdiag[k] << endl;
      tdiag[k] = 0;
    }
    out << "-----------------------------------------" << endl;
  }

};

#endif
