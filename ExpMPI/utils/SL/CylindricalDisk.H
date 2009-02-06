#include <vector>
#include <cmath>

using namespace std;

#include <Vector.h>
#include <massmodel.h>

#include <SLSphere.H>

/** Functionoid defining the radial and vertical profiles for the
    disk.  The default profile is exponential with a sech^2 vertical
    profile.  The three input parameters are:

    0: disk mass
    1: scale length
    2: scale height
*/
class DiskDensityProfile
{
private:
  double M;
  double A;
  double H;
  double norm;

public:
  DiskDensityProfile() : M(0.1), A(0.01), H(0.001) 
  {
    norm = 4.0*M_PI*A*A*H;
  }
  virtual double Density(vector<double>& x)
  {
    
    double z = 1.0/cosh(x[2]/H);
    double R = sqrt(x[0]*x[0]+x[1]*x[1]);
    return exp(-R/A)*z*z/norm;
  }
  virtual void SetParameters(vector<double>& p) 
  {
    if (p.size()<3) 
      throw "DiskDensityProfile: parameter vector size must be >=3";
    M = p[0];
    A = p[1];
    H = p[2];
    norm = 4.0*M_PI*A*A*H;
  }
  virtual ~DiskDensityProfile() = 0;
};

inline DiskDensityProfile::~DiskDensityProfile() {}

/**
   Class to compute multipole profiles for a cylindrical disk
*/
class CylindricalDisk : public DiskDensityProfile
{
private:
  SphericalModelTable *model;
  SLSphere *ortho;

  vector<Vector> coefs;
  double dR;
  double rmin, rmax;
  int lmax, nmax, numr;

public:

  // Constructor
  CylindricalDisk(double rmin, double rmax, 
		  int nmax, int lmax, int numr, int numt);

  // Destructor
  ~CylindricalDisk();

  // Evaluation of density profile
  double density_eval(double x, double y, double z);

  // Evaluation of potential profile
  double potential_eval(double x, double y, double z);

  // Dump coefficients
  void dump_coefficients(string filename);
};


