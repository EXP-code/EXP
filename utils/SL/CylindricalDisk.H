#include <vector>
#include <memory>
#include <cmath>

using namespace std;

#include "massmodel.H"
#include "SLSphere.H"

/** 
    Functionoid defining the radial and vertical profiles for the disk
*/
class DiskDensityProfile
{
public:
  virtual double Density(vector<double>& x)     = 0;
  virtual void SetParameters(vector<double>& p) = 0;
  virtual ~DiskDensityProfile() = 0;
};

inline DiskDensityProfile::~DiskDensityProfile() {}

/**
   Class to compute multipole profiles for a cylindrical disk
   with an exponential surface density and sech^2 vertical density

   The default profile is exponential with a sech^2 vertical
   profile.  The three input parameters are:

   0: disk mass
   1: scale length
   2: scale height
*/
class CylindricalDisk : public DiskDensityProfile
{
private:
  double M;
  double A;
  double H;
  double normS, norm;

protected:

  std::shared_ptr<SphericalModelTable> model;
  std::shared_ptr<SLSphere> ortho;

  std::vector<Eigen::VectorXd> coefs;
  double rmin, rmax;
  int lmax, nmax, numr, numt, numg;

  void make_grids();
  bool get_interp(double x, double z);

  double dR, rgmin, rgmax;
  vector< vector< vector<double> > > frc;
  vector<int> indx;
  vector<double> aint;
  bool grid, linear;
  bool initialized;

public:

  /// Constructor
  CylindricalDisk() : M(0.1), A(0.01), H(0.001), 
		      normS(M/(2.0*M_PI*A*A)), norm(normS/(2.0*H)),
		      grid(false), initialized(false) {}

  /// Destructor
  virtual ~CylindricalDisk();

  /// Compute the force grid
  void Initialize(double rmin, double rmax, bool logr,
		  int nmax, int lmax, int numr, int numt, int numg,
		  vector<double> param);

  /// Evaluation of density profile
  double density_eval(double x, double y, double z);

  /// Evaluation of potential profile
  double potential_eval(double x, double y, double z);

  /// Evaluation of force on table
  void force_eval(double x, double y, double z, double &fr, double &fz);

  /// Evaluation of potential and density table
  void pot_dens_eval(double x, double y, double z, double &potl, double &dens);

  /// Dump coefficients
  void dump_coefficients(string filename);

  /// The disk volume density
  virtual double Density(vector<double>& x)
  {
    double z = 1.0/cosh(x[2]/H);
    double R = sqrt(x[0]*x[0]+x[1]*x[1]);
    return norm * exp(-R/A)*z*z;
  }

  /// The disk surface density
  virtual double Density(double r)
  {
    return normS * exp(-r/A);
  }

  /// For setting the disk parameters
  virtual void SetParameters(vector<double>& p) 
  {
    if (p.size()<3) 
      throw "CylindricalDisk: parameter vector size must be >=3";

    M     = p[0];
    A     = p[1];
    H     = p[2];
    normS = M/(2.0*M_PI*A*A);
    norm  = normS/(2.0*H);
  }

};


/**
   Class to compute multipole profiles for a cylindrical disk
   with an Kuzmin surface density and sech^2 vertical density
*/
class KuzminDisk : public CylindricalDisk
{
private:
  double V, A, H, normS, norm;

public:

  KuzminDisk() : V(1.0), A(0.01), H(0.001), 
		 normS(V*V/(2.0*M_PI*A)), norm(normS/(2.0*H)) {}

  virtual double Density(vector<double>& x)
  {
    double z = 1.0/cosh(x[2]/H);
    double R = sqrt(x[0]*x[0]+x[1]*x[1])/A;
    return norm*pow(1.0 + R*R, -1.5)*z*z;
  }

  virtual double Density(double r)
  {
    double R = r/A;
    return normS*pow(1.0 + R*R, -1.5);
  }

  virtual void SetParameters(vector<double>& p) 
  {
    if (p.size()<3) 
      throw "KuzminDisk: parameter vector size must be >=3";

    V     = p[0];
    A     = p[1];
    H     = p[2];
    normS = V*V/(2.0*M_PI*A);
    norm  = normS/(2.0*H);
  }

};
