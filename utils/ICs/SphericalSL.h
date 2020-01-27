// This may look like C code, but it is really -*- C++ -*-

#ifndef _SphericalSL_H
#define _SphericalSL_H

#include <vector>

#include <Particle.H>
#include <Vector.h>
#include <logic.h>
#include <Particle.H>
#include <SLGridMP2.h>

class SphericalSL
{
private:
  int NMAX, LMAX, n, nthrds;

  double *MPIin, *MPIout;
  bool MPIset;

  SLGridSph *ortho;

  Matrix expcoef, *cc;

  std::vector<Matrix> expcoef1;
  std::vector<Matrix*> cc1;

  Matrix factorial;
  Matrix normM, krnl;
  std::vector<Matrix> potd, legs, dpot, dlegs, dend;
  std::vector<Vector> cosm, sinm;

  bool coefs_made;

  bool initialized;
  void initialize(void);

  std::vector<int> use;
  int compute, used;


  void compute_coefficients_single(vector<Particle> &part);
  void compute_coefficients_thread(vector<Particle> &part);
  void compute_coefficients_thread_call(int id, vector<Particle>* p);

  void parallel_gather_coefficients(Matrix& expcoef, Matrix& expcoef1,
				    Matrix*& cc, Matrix*& cc1,
				    int lmax);

  void parallel_distribute_coefficients(Matrix& expcoef, int lmax);

  void pca_hall(int compute, Matrix& expcoef, Matrix*& cc, Matrix& normM);

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
  SphericalSL(int Nth, int lmax, int nmax, int CMAP, double SCALE=1.0);
  ~SphericalSL(void);

  void reset(int Nth, int lmax, int nmax, int CMAP, double SCALE);

				// Parameter access
  int get_maxNR(void) {return NMAX;}
  int get_maxNL(void) {return LMAX;}
  SLGridSph* SL(void) {return ortho;}

				// Main member functions

  void accumulate(vector<Particle>& particles);

  
  void determine_fields_at_point(double r, double theta, double phi,
				 double *tdens, double *tpotl, 
				 double *tpotr, double *tpott, 
				 double *tpotp, int id=0);
  

  void get_potl_dens(int l, int n, double r, int id=0)
  {
    ortho->get_dens (dend[id], r);
    ortho->get_pot  (potd[id], r);
    ortho->get_force(dpot[id], r);
  }

  void get_dens_coefs(int l, Vector& coef, double *p, int id=0);
  void get_pot_coefs (int l, Vector& coef, double *p, double *dp, int id=0);

  void dump_coefs(ofstream& out);
  void dump_basis(string& dumpname);

  void set_compute()   { compute = 1;}
  void unset_compute() { compute = 0;}

};

typedef boost::shared_ptr<SphericalSL> SphericalSLptr;

#endif
