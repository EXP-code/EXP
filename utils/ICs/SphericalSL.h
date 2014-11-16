// This may look like C code, but it is really -*- C++ -*-

#ifndef _SphericalSL_H
#define _SphericalSL_H

#include <vector>

#include <Vector.h>
#include <logic.h>
#include <Particle.H>
#include <SLGridMP2.h>

class SphericalSL
{
private:
  int NMAX, LMAX, n;

  double *MPIin, *MPIout;
  bool MPIset;

  SLGridSph *ortho;

  Vector cosm, sinm;
  Matrix expcoef, expcoef1, *cc, *cc1;
  Matrix factorial, potd, dpot, legs, dlegs;
  Matrix normM, krnl, dend;

  bool coefs_made;

  bool initialized;
  void initialize(void);

  int use, compute, used;


  void compute_coefficients(vector<Particle> &part);

  void parallel_gather_coefficients(Matrix& expcoef, Matrix& expcoef1,
				    Matrix*& cc, Matrix*& cc1,
				    int lmax);

  void parallel_distribute_coefficients(Matrix& expcoef, int lmax);

  void pca_hall(int compute, Matrix& expcoef, Matrix*& cc, Matrix& normM);

  void get_pot_coefs_safe(int l, Vector& coef, double *p, double *dp,
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

  static int selector;		// For PCA methods
  static int tk_type;
  static double tksmooth;
  static double tkcum;



				// Constructors
  SphericalSL(void);
  SphericalSL(int lmax, int nmax, double SCALE=1.0);
  ~SphericalSL(void);

  void reset(int lmax, int nmax, double SCALE);

				// Parameter access
  int get_maxNR(void) {return NMAX;}
  int get_maxNL(void) {return LMAX;}
  SLGridSph* SL(void) {return ortho;}

				// Main member functions

  void accumulate(vector<Particle>& particles);
  void determine_fields_at_point(
				 double r, double theta, double phi,
				 double *tdens, double *tpotl, 
				 double *tpotr, double *tpott, 
				 double *tpotp);
  

  void get_potl_dens(int l, int n, double r) {
    ortho->get_dens(dend, r);
    ortho->get_pot(potd, r);
    ortho->get_force(dpot, r);
  }

  void get_dens_coefs(int l, Vector& coef, double *p);
  void get_pot_coefs(int l, Vector& coef, double *p, double *dp);

  void dump_coefs(ofstream& out);
  void dump_basis(string& dumpname);

  void set_compute() { compute = 1;}
  void unset_compute() { compute = 0;}

};


#endif
