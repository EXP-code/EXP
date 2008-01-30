#ifndef _Collide_H
#define _Collide_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>

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

public:

  typedef pair<double, vector<double> > Precord;

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
  Normal *norm;

  vector<Precord> prec;

protected:

  int nthrds;

  std::vector<unsigned> numcnt;
  std::vector<unsigned> colcnt;
  std::vector<double>   tsrat, tdens, tvolc, ttemp, tdelt;

  static const unsigned Nphase = 4, Nmfp = 5;
  std::vector<Precord> tphase, tmfpst;
  std::vector< vector<Precord> > tphaseT, tmfpstT;

  std::vector< vector<unsigned> > numcntT;
  std::vector< vector<unsigned> > colcntT;
  std::vector<unsigned> error1T;
  std::vector<unsigned> col1T, epsm1T, Nepsm1T;
  std::vector<double> KEtotT, KElostT;
  std::vector< vector<double> >  mfpratT, tsratT, tdensT, tvolcT, ttempT, tdeltT;

  const static unsigned numdiag = 9;
  std::vector<unsigned> tdiag, tdiag1, tdiag0;
  std::vector< vector<unsigned> > tdiagT;

  std::vector< vector<pCell*> > cellist;

  double diam0, diam;
  unsigned coltot;
  unsigned epsmtot;
  unsigned epsmcells;
  unsigned errtot;
  unsigned numtot;

  int use_temp, use_dens, use_delt;

  virtual void initialize_cell(pCell* cell, double vrelmax,
			       double dt, double number, int id=0) { }

  virtual int inelastic(pHOT* tree, Particle* p1, Particle* p2,
			double *cr, int id=0) { return 0; }

  virtual void EPSM(pHOT* tree, pCell* cell, int id);
  vector<double> coolrate;

  void debug_list(pHOT& tree);

public:

  static int CNUM;		// Target CNUM collisions if nonzero
				// default: 0

  static bool CBA;		// Use the Consistent Boltzmann Algorithm
				// Alexander, Garcia, Alder 1995, PRL, 74, 5212

  static bool SORTED;		// Print out sorted cell parameters

  static bool PHASE;		// Print out T-rho plane for cells will mass weighting
  
				// Use EPSM when mean free path/cell size ratio 
				// is below EPSMratio (default = -1.0, never)
  static double EPSMratio;
				// Temperature floor in resampling dispersion for EPSM
  static double TFLOOR;


  Collide(double diameter, int nth=1);
  ~Collide();

  virtual void 
  set_temp_dens(int Tpos, int Dpos) { use_temp=Tpos; use_dens=Dpos; }
  virtual void set_timestep(int DTpos) { }

  unsigned total()  { unsigned t=coltot; coltot=0; return t; }
  unsigned EPSMtotal()  { unsigned t=epsmtot; epsmtot=0; return t; }
  unsigned EPSMcells()  { unsigned t=epsmcells; epsmcells=0; return t; }
  unsigned errors() { unsigned t=errtot; errtot=0; return t; }
  unsigned cellN()  { return numtot; }
  unsigned medianNumber();
  unsigned medianColl();
  unsigned collide(pHOT& cl, double Fn, double tau);
  void collQuantile(vector<double>& quantiles, vector<double>& coll_);
  void mfpsizeQuantile(vector<double>& quantiles, vector<double>& mfp_, vector<double>& ts_) ;

  virtual double Etotal() { cerr << "Oops [1]" << endl; return 0.0; }
  virtual double Elost()  { cerr << "Oops [2]" << endl; return 0.0; }

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
