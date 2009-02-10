#include <vector>
#include <cmath>

using namespace std;

#include <Vector.h>
#include <massmodel.h>

#include <SLSphere.H>

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
  double norm;

  SphericalModelTable *model;
  SLSphere *ortho;

  vector<Vector> coefs;
  double rmin, rmax;
  int lmax, nmax, numr, numt, numg;

  void make_grids();
  bool get_interp(double x, double z);

  double dR, rgmin, rgmax;
  vector< vector< vector<double> > > frc;
  vector<int> indx;
  vector<double> aint;
  bool grid, linear;

public:

  /// Constructor
  CylindricalDisk(double rmin, double rmax, 
		  int nmax, int lmax, int numr, int numt, int numg,
		  vector<double> param);

  /// Destructor
  ~CylindricalDisk();

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

  /// The disk density
  virtual double Density(vector<double>& x)
  {
    
    double z = 1.0/cosh(x[2]/H);
    double R = sqrt(x[0]*x[0]+x[1]*x[1]);
    return exp(-R/A)*z*z/norm;
  }

  /// For setting the disk parameters
  virtual void SetParameters(vector<double>& p) 
  {
    if (p.size()<3) 
      throw "DiskDensityProfile: parameter vector size must be >=3";
    M = p[0];
    A = p[1];
    H = p[2];
    norm = 4.0*M_PI*A*A*H/M;
  }

};


/**
   Class to compute multipole profiles for a cylindrical disk
   with an Mestel surface density and sech^2 vertical density
*/
class MestelDisk : public CylindricalDisk
{
private:
  double V, A, H, norm;

public:

  MestelDisk(double rmin, double rmax, 
	     int nmax, int lmax, int numr, int numt, int numg,
	     vector<double> param) :
    V(1.0), A(0.01), H(0.001), 
    CylindricalDisk(rmin, rmax, nmax, lmax, numr, numt, numg, param)
  {
    norm = V*V/(4.0*M_PI*H);
  }

  virtual double Density(vector<double>& x)
  {
    double z = 1.0/cosh(x[2]/H);
    double R = sqrt(x[0]*x[0]+x[1]*x[1])/A;
    return norm*pow(1.0 + R*R, -1.5)*z*z;
  }

  virtual void SetParameters(vector<double>& p) 
  {
    if (p.size()<3) 
      throw "MestelDisk: parameter vector size must be >=3";
    V = p[0];
    A = p[1];
    H = p[2];
    norm = V*V/(4.0*M_PI*H);
  }

};
