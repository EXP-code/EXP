// This may look like C code, but it is really -*- C++ -*-

#ifndef _massmodel_h
#ifdef __GNUG__
#pragma interface
#endif
#define _massmodel_h 1

const char rcsid_massmodel[] = "$Id$";

#include <ACG.h>
#include <Uniform.h>

#include <logic.h>
#include <Vector.h>

class QPDistF;

struct RUN {
  Vector x;
  Vector y;
  Vector y2;
  int num;
};

struct FDIST {
  Vector Q;
  Vector fQ;
  Vector ffQ;
  Vector fQ2;
  Vector ffQ2;
  double ra2;
  double off;
  int num;
};

class MassModel
{
public:

  int dim;
  String ModelID;
  Logic defined;

  virtual double get_mass(const double, const double, const double) = 0;
  virtual double get_density(const double, const double, const double) = 0;
  virtual double get_pot(const double, const double, const double) = 0;

  const int dof() { return dim; }

				// Error function

  void bomb(char *s) {
    cerr << "ERROR from " << ModelID << ": " << s << '\n';
    exit(-1);
  }

};
  
class AxiSymModel : public MassModel
{
private:
                                // Stuff for gen_point
  ACG *gen;
  Uniform *Unit;
  Logic gen_firstime;
  Vector gen_rloc, gen_mass, gen_fmax;
  int gen_seed, gen_N, gen_itmax;

  Vector gen_point_2d(int& ierr);
  Vector gen_point_3d(int& ierr);

public:
  Logic dist_defined;

  AxiSymModel(void) { 
    ModelID = "AxiSymModel";
    gen_firstime = TRUE;
    gen_seed = 11;
    gen_N = 400;
    gen_itmax=4000;
  }

  void set_seed(int s) { gen_seed = s;}
  void set_itmax(int s) { gen_itmax = s;}

  virtual double get_mass(const double) = 0;
  virtual double get_density(const double) = 0;
  virtual double get_pot(const double) = 0;
  virtual double get_dpot(const double) = 0;
  virtual double get_dpot2(const double) = 0;
  virtual void get_pot_dpot(const double, double&, double&) = 0;

  double get_mass(const double x1, const double x2, const double x3)
    { return get_mass(sqrt(x1*x1 + x2*x2 + x3*x3)); }

  double get_density(const double x1, const double x2, const double x3)
    { return get_density(sqrt(x1*x1 + x2*x2 + x3*x3)); }

  double get_pot(const double x1, const double x2, const double x3)
    { return get_pot(sqrt(x1*x1 + x2*x2 + x3*x3)); }

  // Addiional member functions

  virtual double get_min_radius(void) = 0;
  virtual double get_max_radius(void) = 0;
  virtual double distf(double, double) = 0;
  virtual double dfde(double, double) = 0;
  virtual double dfdl(double, double) = 0;
  virtual double d2fde2(double, double) = 0;

  Vector gen_point(int& ierr) {
    if (dof()==2)
      return gen_point_2d(ierr);
    else if (dof()==3)
      return gen_point_3d(ierr);
    else
      bomb( "dof must be 2 or 3" );
  }

};


class EmbeddedDiskModel : public AxiSymModel
{
private:
  AxiSymModel **t;
  double *m_scale;
  double *r_scale;
  int number;

  double rmin, rmax;

  QPDistF *df;

public:
  EmbeddedDiskModel(AxiSymModel **T, double *M_scale, double *R_scale, 
		      int NUMBER);

  double get_mass(const double);
  double get_density(const double);
  double get_pot(const double);
  double get_dpot(const double);
  double get_dpot2(const double);
  void get_pot_dpot(const double, double&, double&);

  double get_mass(const double x1, const double x2, const double x3)
    { return get_mass(sqrt(x1*x1 + x2*x2 + x3*x3)); }

  double get_density(const double x1, const double x2, const double x3)
    { return get_density(sqrt(x1*x1 + x2*x2 + x3*x3)); }

  double get_pot(const double x1, const double x2, const double x3)
    { return get_pot(sqrt(x1*x1 + x2*x2 + x3*x3)); }

  // Additional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }

  void setup_df(int egrid=10, int kgrid=5, int mgrid=20,
		double lambda=0.0, double alpha=-6.0, double beta=1.0,
		double gama=1.0, double sigma=2.0, double rmmax=-1, 
		double roff=0.05, double eoff=0.5, double koff=0.5, 
		int nint=20, int numt=20);

  // Read in from state file

  void setup_df(String& file);

  void verbose_df(void);
  double distf(double E, double L);
  double dfde(double E, double L);
  double dfdl(double E, double L);
  double d2fde2(double E, double L);
  void save_df(String& file);

};


class SphericalModelTable : public AxiSymModel
{
private:
  struct RUN mass;
  struct RUN density;
  struct RUN pot;
  struct FDIST df;
  int num;
  int numdf;
  int num_params;
  double *params;
  double diverge_rfac;
  int diverge;
  int external;

public:

  SphericalModelTable(String filename, 
		 int DIVERGE = 0, double DIVERGE_RFAC = 1.0, int EXTERNAL = 0);

  SphericalModelTable::SphericalModelTable(int num, 
		 double *r, double *d, double *m, double *p,
	         int DIVERGE = 0, double DIVERGE_RFAC = 0, int EXTERNAL = 0,
		 String ID = "" );

//  ~SphericalModelTable();

  // Required member functions

  double get_mass(const double);
  double get_density(const double);
  double get_pot(const double);
  double get_dpot(const double);
  double get_dpot2(const double);
  void get_pot_dpot(const double, double&, double&);
  
  // Additional member functions

  const int get_num_param(void) { return num_params; }
  const double get_param(int i) { return params[i-1]; }

  double get_min_radius(void) { return mass.x[1]; }
  double get_max_radius(void) { return mass.x[mass.num]; }
  int grid_size(void) { return num; }

  void setup_df(int NUM, double RA=1.0e20);
  const double get_ra2(void) { return df.ra2; }

  double distf(double E, double L);
  double dfde(double E, double L);
  double dfdl(double E, double L);
  double d2fde2(double E, double L);
};

#include <QPDistF.h>

#endif
