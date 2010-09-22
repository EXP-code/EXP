// This may look like C code, but it is really -*- C++ -*-

#ifndef _massmodel_h
#define _massmodel_h 1

const char rcsid_massmodel[] = "$Id$";

#include <vector>

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>

#include <Vector.h>
#include <orbit.h>

class QPDistF;

class RUN 
{
public:
  Vector x;
  Vector y;
  Vector y2;
  int num;
};

class FDIST 
{
public:
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
  string ModelID;
  bool defined;

  virtual ~MassModel() {}

  virtual double get_mass(const double, const double, const double) = 0;
  virtual double get_density(const double, const double, const double) = 0;
  virtual double get_pot(const double, const double, const double) = 0;

  int dof() { return dim; }

				// Error function

  void bomb(const char *s) {
    cerr << "ERROR from " << ModelID << ": " << s << '\n';
    exit(-1);
  }

};
  

class AxiSymModel : public MassModel
{
protected:
  // Stuff for gen_point
  ACG *gen;
  Uniform *Unit;
  Normal *Gauss;
  bool gen_firstime;
  bool gen_firstime_E;
  bool gen_firstime_jeans;
  Vector gen_rloc, gen_mass, gen_fmax;
  SphericalOrbit gen_orb;
  double gen_fomax;
  
  Vector gen_point_2d(int& ierr);
  Vector gen_point_2d(double r, int& ierr);
  Vector gen_point_3d(int& ierr);
  Vector gen_point_3d(double Emin, double Emax, double Kmin, double Kmax, int& ierr);
  Vector gen_point_jeans_3d(int& ierr);
  
  double Emin_grid, Emax_grid, dEgrid, dKgrid;
  vector<double> Egrid, Kgrid, EgridMass;
  vector<double> Jmax;
  
  class WRgrid
  {
  public:
    vector<double> w1;
    vector<double> r;
  };
  
  typedef vector<WRgrid> wrvector;
  vector<wrvector> Rgrid;
  
public:
  // Stuff for gen_point
  static bool gen_EJ;
  static int numr, numj;
  static int gen_N;
  static int gen_E;
  static int gen_K;
  static int gen_itmax;
  static int gen_logr;
  static double gen_rmin;
  static double gen_kmin;
  static double gen_tolE, gen_tolK;
  static unsigned int gen_seed;
  
  bool dist_defined;
  
  AxiSymModel(void) { 
    ModelID = "AxiSymModel";
    gen_firstime = true;
    gen_firstime_E = true;
    gen_firstime_jeans = true;
  }
  
  virtual ~AxiSymModel() {}
  
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
  
  virtual Vector gen_point(int& ierr) {
    if (dof()==2)
      return gen_point_2d(ierr);
    else if (dof()==3)
      return gen_point_3d(ierr);
    else
      bomb( "dof must be 2 or 3" );
    
    return Vector();
  }
  
  virtual Vector gen_point_jeans(int& ierr) {
    if (dof()==2)
      bomb( "AxiSymModel: gen_point_jeans_2d(ierr) not implemented!" );
    else if (dof()==3)
      return gen_point_jeans_3d(ierr);
    else
      bomb( "dof must be 2 or 3" );
    
    return Vector();
  }
  
  virtual Vector gen_point(double r, int& ierr) {
    if (dof()==2)
      return gen_point_2d(r, ierr);
    else if (dof()==3)
      bomb( "AxiSymModel: gen_point(r, ierr) not implemented!" );
    else
      bomb( "AxiSymModel: dof must be 2 or 3" );
    
    return Vector();
  }
  
  virtual Vector gen_point(double Emin, double Emax, double Kmin, double Kmax, int& ierr) {
    if (dof()==2)
      bomb( "AxiSymModel: gen_point(r, ierr) not implemented!" );
    else if (dof()==3)
      return gen_point_3d(Emin, Emax, Kmin, Kmax, ierr);
    else
      bomb( "AxiSymModel: dof must be 2 or 3" );
    
    return Vector();
  }
  
  virtual void gen_velocity(double *pos, double *vel, int& ierr);
  
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

  virtual ~EmbeddedDiskModel();

  virtual double get_mass(const double);
  virtual double get_density(const double);
  virtual double get_pot(const double);
  virtual double get_dpot(const double);
  virtual double get_dpot2(const double);
  virtual void get_pot_dpot(const double, double&, double&);

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
		double kmin=0.0, double kmax=1.0,
		int nint=20, int numt=20);

  // Read in from state file

  void setup_df(string& file);

  void verbose_df(void);
  double distf(double E, double L);
  double dfde(double E, double L);
  double dfdl(double E, double L);
  double d2fde2(double E, double L);
  void save_df(string& file);

};


class SphericalModelTable : public AxiSymModel
{
private:
  RUN mass;
  RUN density;
  RUN pot;
  struct FDIST df;
  int num;
  int numdf;
  int num_params;
  double *params;
  double diverge_rfac;
  int diverge;
  int external;

public:

  static int count;		// Count instantiations
  static int even;		// Assume even spacing (default: yes)
  static int logscale;		// Log scale in df computation (default: yes)
  static int linear;		// Linear interpolation in model (default: no)

  SphericalModelTable(string filename, 
		 int DIVERGE = 0, double DIVERGE_RFAC = 1.0, int EXTERNAL = 0);

  SphericalModelTable(int num, 
		 double *r, double *d, double *m, double *p,
	         int DIVERGE = 0, double DIVERGE_RFAC = 0, int EXTERNAL = 0,
					   string ID = "" );

  virtual ~SphericalModelTable();

  // Required member functions

  virtual double get_mass(const double);
  virtual double get_density(const double);
  virtual double get_pot(const double);
  virtual double get_dpot(const double);
  virtual double get_dpot2(const double);
  virtual void get_pot_dpot(const double, double&, double&);
  
  // Additional member functions

  int get_num_param(void) { return num_params; }
  double get_param(int i) { return params[i-1]; }

  double get_min_radius(void) { return mass.x[1]; }
  double get_max_radius(void) { return mass.x[mass.num]; }
  int grid_size(void) { return num; }

  void setup_df(int NUM, double RA=1.0e20);
  void print_model(char const *name);
  void print_df(char const *name);
  double get_ra2(void) { return df.ra2; }

  double distf(double E, double L);
  double dfde(double E, double L);
  double dfdl(double E, double L);
  double d2fde2(double E, double L);
};


class SphericalModelMulti : public AxiSymModel
{
protected:
  AxiSymModel* real;
  AxiSymModel* fake;

  SphericalOrbit orb, gen_orb;

  double rmin_gen, rmax_gen;

public:

  SphericalModelMulti(AxiSymModel* Real, AxiSymModel* Fake);

  // Required member functions

  double get_mass(const double r) { return real->get_mass(r); }
  double get_density(const double r) { return real->get_density(r); }
  double get_pot(const double r) { return real->get_pot(r); }
  double get_dpot(const double r)  { return real->get_dpot(r); }
  double get_dpot2(const double r)  { return real->get_dpot2(r); }
  void get_pot_dpot(const double r, double& p, double& dp) 
    { real->get_pot_dpot(r, p, dp); }
  
  // Additional member functions

  double get_min_radius(void) { return real->get_min_radius(); }
  double get_max_radius(void) { return real->get_max_radius(); }

  double distf(double E, double L)  { return real->distf(E, L); }
  double dfde(double E, double L)   { return real->dfde(E, L); }
  double dfdl(double E, double L)   { return real->dfdl(E, L); }
  double d2fde2(double E, double L) { return real->d2fde2(E, L); }

  // Overloaded to provide mass distribution from Real and Number distribution from Fake
  Vector gen_point(int& ierr);
  Vector gen_point(double r, int& ierr);
  Vector gen_point(double Emin, double Emax, double Kmin, double Kmax, int& ierr);


  // Set new minimum and maximum for realization
  void set_min_radius(const double& r) { rmin_gen = r; }
  void set_max_radius(const double& r) { rmax_gen = r; }
};


#include <QPDistF.h>

#endif
