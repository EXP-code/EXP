// This may look like C code, but it is really -*- C++ -*-

#ifndef _SphericalSL_h

#define _SphericalSL_h 1

#include <stdio.h>
#include <vector.h>
#include <logic.h>
#include <SLGridMP2.h>
#include <any_thread2.h>

/*
static struct CoefHeader {
  double time;
  int mmax;
  int nord;
  int nmax;
} coefheader;

static struct CoefHeader2 {
  double time;
  int mmax;
  int nmax;
} coefheader2;
*/

class SphericalSL : public any_thread
{
private:
  int NMAX, LMAX, n;

  double *MPIin, *MPIout;
  bool MPIset;

  SLGridSph *ortho;

  Vector *cosm, *sinm;
  Matrix expcoef, expcoef1, *expcoef0, *cc, *cc1;
  Matrix factorial, *potd, *dpot, *legs, *dlegs;
  Matrix normM, krnl, dend;

  bool coefs_made;

  bool initialized;
  void initialize(void);

  pthread_mutex_t cc_lock;
  int *use, compute;
  bool do_accel;
  void thread_call(void*);
  void determine_coefficients_SLsph_thread(void * arg);
  void determine_acceleration_and_potential_SLsph_thread(void * arg);

  void get_pot_coefs_SLsph_safe(int l, Vector& coef, double *p, double *dp,
				Matrix& potd1, Matrix& dpot1);


				// Error function

  void bomb(const char *s) {
    cerr << "SphericalSL ERROR: " << s << '\n';
    exit(-1);
  }


public:
				// Global parameters

  static double RMIN;
  static double RMAX;
  static int NUMR;

				// Constructors
  SphericalSL(void);
  SphericalSL(int lmax, int nmax);
  ~SphericalSL(void);

  void reset(int lmax, int nmax);

				// Parameter access
  int get_maxNR(void) {return NMAX;}
  int get_maxNL(void) {return LMAX;}
  SLGridSph* SL(void) {return ortho;}

				// Main member functions

  void determine_coefficients_SLsph(void);
  void get_acceleration_and_potential_SLsph(void);
  void determine_acceleration_and_potential_SLsph(void);

  void determine_fields_at_point_SLsph(
				       double r, double theta, double phi,
				       double *tdens, double *tpotl, 
				       double *tpotr, double *tpott, 
				       double *tpotp);


  void get_potl_dens_SLsph(int l, int n, double r) {
    ortho->get_dens(dend, r);
    ortho->get_pot(potd[0], r);
    ortho->get_force(dpot[0], r);
  }

  void get_dens_coefs_SLsph(int l, Vector& coef, double *p);
  void get_pot_coefs_SLsph(int l, Vector& coef, double *p, double *dp);

  void dump_coefs_SLsph(FILE *fout);
};


#endif
