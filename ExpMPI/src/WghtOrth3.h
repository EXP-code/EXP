// This may look like C code, but it is really -*- C++ -*-

#ifndef _WghtOrth_h
#define _WghtOrth_h 1

#include <SLGridMP2.h>

#include <stdio.h>
#include <logic.h>

class WghtOrtho
{
public:

				// Error function

  void bomb(const char *s) {
    cerr << "WghtOrtho ERROR: " << s << '\n';
    exit(-1);
  }

};


class CylindricalSL : public WghtOrtho
{
private:
  double dk, dZ;
  int NMAX, NFFT;
  int MMAX, n;

  double *MPIin, *MPIout;
  bool MPIset;

  SLGridCyl *ortho;

  Vector zn;
  Matrix* accum_cos;
  Matrix* accum_sin;
  Matrix* table;
  Matrix* tablef;

  bool coefs_made;

  bool initialized;
  void initialize(void);

public:

				// Global parameters

  static double RMIN;
  static double RMAX;
  static int NUMR;
  static double ZMAX;
  static int NINT;


				// Constructors
  CylindricalSL(void);
  CylindricalSL(int mu, int nfft, int m);
  ~CylindricalSL(void);

  void reset(int mu, int nfft, int m);

				// Parameter access
  int get_maxNR(void) {return NMAX;}
  int get_maxNZ(void) {return n;}
  int get_maxM(void) {return MMAX;}
  SLGridCyl* SL(void) {return ortho;}

				// Main member functions
  double get_pot(int m, int k, int n, double r, double z);
  double get_dens(int m, int k, int n, double r, double z);

  Matrix inner_product_dens(int m, double (*fct)(double, double));
  Matrix inner_product_pot(int m, double (*fct)(double, double));

  double pot_eval(Matrix& coef, int m, double r, double z);
  double dens_eval(Matrix& coef, int m, double r, double z);
  double r_force_eval(Matrix& coef, int m, double r, double z);
  double z_force_eval(Matrix& coef, int m, double r, double z);
  void force_eval(Matrix& coef, int m, double r, double z, double& fr, double& fz);

  void setup_accumulation(void);
  void accumulate(double r, double z, double phi, double mass);
  void make_coefficients(void);
  void accumulated_eval(double r, double z, double phi,
			double& p, double& fr, double& fz, double& fp);
  double accumulated_dens_eval(double r, double z, double phi);
  double accumulated_vertical_mass_eval(double r, double z, double phi);
  double accumulated_potential_energy(void);
  
  void dump_coefs(FILE *fout);
  void dump_coefs_binary(FILE *fout, double time);

  void read_coefs(istream *in);
  void read_coefs_binary(istream *in, double* time);
};

#endif
