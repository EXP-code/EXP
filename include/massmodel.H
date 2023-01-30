#ifndef _massmodel_H
#define _massmodel_H

#include <vector>
#include <memory>
#include <random>

#include <Eigen/Eigen>

#include <interp.H>
#include <orbit.H>

class QPDistF;

//! A one-dimensional profile
class RUN 
{
public:
  Eigen::VectorXd x;
  Eigen::VectorXd y;
  Eigen::VectorXd y2;
  int num;

  //! Simple constructor.  Use value of num to check for initialization.
  RUN() : num(0) {}

  void resize(int N)
  {
    num = N;
    x.resize(N);
    y.resize(N);
    y2.resize(N);
  }
};

//! A one-dimensional Merritt-Osipkov distribution
class FDIST 
{
public:
  Eigen::VectorXd Q;
  Eigen::VectorXd fQ;
  Eigen::VectorXd ffQ;
  Eigen::VectorXd fQ2;
  Eigen::VectorXd ffQ2;
  double ra2;
  double off;
  int num;

  //! Simple constructor.  Use value of num to check for initialization.
  FDIST() : num(0) {}

};

//! Chebyshev version
class FDISTC
{
public:

  Eigen::VectorXd  Q;
  Eigen::VectorXd  fQ;
  Eigen::VectorXd  ffQ;
  Cheby1d FF, GG;
  double ra2;
  double off;
  int num;

  //! Simple constructor.  Use value of num to check for initialization.
  FDISTC() : num(0) {}

};

//! A three-dimensional density-potential model
class MassModel

{
public:

  //! The number of dimensions (degrees of freedom)
  int dim;
  //! Model identifier
  string ModelID;
  //! True if model is assigned
  bool defined;

  //! Destructor
  virtual ~MassModel() {}

  /** Get the mass inside of position (x, y, z) (this is ambiguous in
      general and will depend on the model) */
  virtual double get_mass(const double, const double, const double) = 0;
  //! The density at (x, y, z)
  virtual double get_density(const double, const double, const double) = 0;
  //! The potential at (x, y, z)
  virtual double get_pot(const double, const double, const double) = 0;

  //! Return the number of dimensions (degrees of freedom)
  int dof() { return dim; }

  //! Exception handling using std::exception
  void bomb(const char *s) {
    std::ostringstream sout;
    sout << "MassModel error from " << ModelID << ": " << s;
    throw std::runtime_error(sout.str());
  }

};
  

//! A specification of a MassModel to the axisymmetric case
class AxiSymModel : public MassModel,
		    public std::enable_shared_from_this<AxiSymModel>
{
protected:
  //@{
  //! Stuff for gen_point

  std::uniform_real_distribution<> Unit;
  std::normal_distribution<> Gauss;

  bool gen_firstime;
  bool gen_firstime_E;
  bool gen_firstime_jeans;
  Eigen::VectorXd gen_rloc, gen_mass, gen_fmax;
  SphericalOrbit gen_orb;
  double gen_fomax, gen_ecut;
  //@}

  Eigen::VectorXd gen_point_2d(int& ierr);
  Eigen::VectorXd gen_point_2d(double r, int& ierr);
  Eigen::VectorXd gen_point_3d(int& ierr);
  Eigen::VectorXd gen_point_3d(double Emin, double Emax, double Kmin, double Kmax, int& ierr);
  Eigen::VectorXd gen_point_jeans_3d(int& ierr);
  
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
  //@{
  //! Stuff for gen_point
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
  //@}
  
  //! The distribution is defined and computed
  bool dist_defined;
  
  //! Constructor
  AxiSymModel(void) { 
    ModelID = "AxiSymModel";
    gen_firstime = true;
    gen_firstime_E = true;
    gen_firstime_jeans = true;
    gen_ecut = 1.0;
  }
  
  //! Destructor
  virtual ~AxiSymModel() {}
  
  //! Maximum number of iterations for obtaining a valid variate
  void set_itmax(int s) { gen_itmax = s;}
  
  //@{
  //! Access the profile
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
  //@}
  
  //@{
  // Addiional member access functions
  virtual double get_min_radius(void) = 0;
  virtual double get_max_radius(void) = 0;
  virtual double distf(double, double) = 0;
  virtual double dfde(double, double) = 0;
  virtual double dfdl(double, double) = 0;
  virtual double d2fde2(double, double) = 0;
  //@}
  
  //! Set cutoff on multimass realization grid
  void set_Ecut(double cut) { gen_ecut = cut; }

  //! Generate a phase-space point
  virtual Eigen::VectorXd gen_point(int& ierr) {
    if (dof()==2)
      return gen_point_2d(ierr);
    else if (dof()==3)
      return gen_point_3d(ierr);
    else
      bomb( "dof must be 2 or 3" );
    
    return Eigen::VectorXd();
  }
  
  //! Generate a phase-space point using Jeans' equations
  virtual Eigen::VectorXd gen_point_jeans(int& ierr) {
    if (dof()==2)
      bomb( "AxiSymModel: gen_point_jeans_2d(ierr) not implemented!" );
    else if (dof()==3)
      return gen_point_jeans_3d(ierr);
    else
      bomb( "dof must be 2 or 3" );
    
    return Eigen::VectorXd();
  }
  
  //! Generate a phase-space point at a particular radius
  virtual Eigen::VectorXd gen_point(double r, int& ierr) {
    if (dof()==2)
      return gen_point_2d(r, ierr);
    else if (dof()==3)
      bomb( "AxiSymModel: gen_point(r, ierr) not implemented!" );
    else
      bomb( "AxiSymModel: dof must be 2 or 3" );
    
    return Eigen::VectorXd();
  }
  
  //! Generate a phase-space point at a particular energy and angular momentum
  virtual Eigen::VectorXd gen_point(double Emin, double Emax, double Kmin, double Kmax, int& ierr) {
    if (dof()==2)
      bomb( "AxiSymModel: gen_point(r, ierr) not implemented!" );
    else if (dof()==3)
      return gen_point_3d(Emin, Emax, Kmin, Kmax, ierr);
    else
      bomb( "AxiSymModel: dof must be 2 or 3" );
    
    return Eigen::VectorXd();
  }
  
  //! Generate a the velocity variate from a position
  virtual void gen_velocity(double *pos, double *vel, int& ierr);
  
};

class EmbeddedDiskModel : public AxiSymModel
{
private:
  std::vector<std::shared_ptr<AxiSymModel>> t;
  std::vector<double> m_scale;
  std::vector<double> r_scale;
  int number;

  double rmin, rmax;

  std::shared_ptr<QPDistF> df;

public:
  EmbeddedDiskModel(std::vector<std::shared_ptr<AxiSymModel>>& T,
		    std::vector<double>& M_scale, std::vector<double>& R_scale, 
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

  void   verbose_df(void);
  double distf(double E, double L);
  double dfde(double E, double L);
  double dfdl(double E, double L);
  double d2fde2(double E, double L);
  void   save_df(string& file);

};


/** Describe a spherical model from a four-column table of radius,
    density, enclosed mass, and gravitational potential

    The file may have any number of leading comment lines, any line
    that starts with a "!" or a "#".  The first non-comment record is
    assumed to be an integer containing the number of lines.  The
    subseuqnet lines are assumed to be records beginning the radius
    and followed by the triple of (density, enclosed mass,
    gravitational potential).  
*/
class SphericalModelTable : public AxiSymModel
{
private:
  RUN mass;
  RUN density;
  RUN pot;
  FDIST df;
  FDISTC dfc;
  int num;
  int numdf;
  int num_params;
  std::vector<double> params;
  double diverge_rfac;
  double foffset;
  int diverge;
  int external;

  void debug_fdist();

public:

  //! Count the number of instantiations (for debugging)
  static int count;		

  //! Assume even spacing in the mass model table (default: yes)
  static int even;		

  //! Chebyshev order (default: 0, no Chebyshev)
  static int chebyN;

  //! Linear interpolation in model (default: no)
  static int linear;		

  //! Constructor (from a table)
  SphericalModelTable(string filename, 
		      int DIVERGE = 0, double DIVERGE_RFAC = 1.0, int EXTERNAL = 0);

  //! Constructor (from array pointers)
  SphericalModelTable(int num, double *r, double *d, double *m, double *p,
		      int DIVERGE = 0, double DIVERGE_RFAC = 0,
		      int EXTERNAL = 0, std::string ID = "" );

  //! Constructor (from vectors)
  SphericalModelTable(std::vector<double>& r, std::vector<double>& d,
		      std::vector<double>& m, std::vector<double>& p,
		      int DIVERGE = 0, double DIVERGE_RFAC = 0,
		      int EXTERNAL = 0, std::string ID = "" );

  //! Destructor
  virtual ~SphericalModelTable();

  //@{
  //! Evaluate the profile at a particular radius
  virtual double get_mass(const double) override;
  virtual double get_density(const double) override;
  virtual double get_pot(const double) override;
  virtual double get_dpot(const double) override;
  virtual double get_dpot2(const double) override;
  virtual void   get_pot_dpot(const double, double&, double&) override;
  //@}
  
  //@{
  //! Get info about the model
  int get_num_param(void) { return num_params; }
  double get_param(int i) { return params[i]; }

  double get_min_radius(void) override { return mass.x[0]; }
  double get_max_radius(void) override { return mass.x[mass.num-1]; }
  int grid_size(void) { return num; }
  void print_model(const std::string& name);
  void print_model_eval(const std::string& name, int number);
  //@}

  //@{
  //! Set up and get info about the phase-space distribution
  void setup_df(int NUM, double RA=1.0e20);
  void print_df(const std::string& name);
  double get_ra2(void) { return df.ra2; }
  //@}

  //@{
  //! Evaluate the distribution function
  double distf(double E, double L) override;
  double dfde(double E, double L) override;
  double dfdl(double E, double L) override;
  double d2fde2(double E, double L) override;
  //@}
};


/** 
    This uses two SphericalModelTables to define a mass profile and a
    number profile for generating a distribution with variable mass
    per particle.
*/
class SphericalModelMulti : public AxiSymModel
{
protected:
  std::shared_ptr<AxiSymModel> real, fake;

  SphericalOrbit orb, gen_orb;

  double rmin_gen, rmax_gen;
  bool noneg;

public:

  //! Constructor
  SphericalModelMulti(std::shared_ptr<AxiSymModel> Real,
		      std::shared_ptr<AxiSymModel> Fake);

  //@{
  //! Evaluate the profile at a radial point
  virtual double get_mass(const double r) override
  { return real->get_mass(r); }

  virtual double get_density(const double r) override
  { return real->get_density(r); }
  
  virtual double get_pot(const double r) override
  { return real->get_pot(r); }

  virtual double get_dpot(const double r)  override
  { return real->get_dpot(r); }

  virtual double get_dpot2(const double r)  override
  { return real->get_dpot2(r); }

  virtual void get_pot_dpot(const double r, double& p, double& dp) override
  { real->get_pot_dpot(r, p, dp); }
  //@}
  
  //@{
  //! Additional member functions
  double get_min_radius(void) override { return real->get_min_radius(); }
  double get_max_radius(void) override { return real->get_max_radius(); }
  //@}

  //@{
  //! Evaluate distribution function and its partial derivatives
  double distf(double E, double L) override { return real->distf(E, L);  }
  double dfde (double E, double L) override { return real->dfde(E, L);   }
  double dfdl (double E, double L) override { return real->dfdl(E, L);   }
  double d2fde2(double E, double L) override { return real->d2fde2(E, L); }
  //@}

  //@{
  //! Overloaded to provide mass distribution from Real and Number distribution from Fake
  Eigen::VectorXd gen_point(int& ierr) override;
  Eigen::VectorXd gen_point(double r, int& ierr) override;
  Eigen::VectorXd gen_point(double Emin, double Emax, double Kmin, double Kmax, int& ierr) override;
  //@}


  //@{
  //! Set new minimum and maximum for realization
  void set_min_radius(const double& r) { rmin_gen = r; }
  void set_max_radius(const double& r) { rmax_gen = r; }
  //@}

  //! Turn off negative mass suppression by NOT requeueing particle
  //! selection if the DF ratio is negative
  void allowNegativeMass() { noneg = false; }

};

typedef std::shared_ptr<AxiSymModel>         AxiSymModPtr;
typedef std::shared_ptr<EmbeddedDiskModel>   EmbDiskModPtr;
typedef std::shared_ptr<SphericalModelMulti> SphModMultPtr;
typedef std::shared_ptr<SphericalModelTable> SphModTblPtr;

#include <QPDistF.H>

#endif
