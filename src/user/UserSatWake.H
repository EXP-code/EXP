#ifndef _UserSatWake_H
#define _UserSatWake_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Halo model based on SphericalModelTable

    @param LMIN:		Minimum spherical harmonic index
    @param LMAX:		Maximum spherical harmonic index
    @param MMIN			Minimum azimuthal spherical harmonic index
    @param MMAX:		Maximum azimuthal spherical harmonic index
    @param lmax:		Maximum L1 index
    @param nmax:		Maximum radial index
    @param numr:		Maximum radial grid for SphereSL
    @param rmin:		Minimum radius for SphereSL
    @param rmax:		Maximum radius for SphereSL
    @param nfreqs:		Number of points in frequency grid for inverse Laplace transform
    @param nint:		Number of points in time grid for Laplace transform segments
    @param scale:		Coordinate scale for SphereSL
    @param HALO_TRUNC:		Truncate halo expansion
    @param nptsE:		Power of 2 quadrature index for Energy integral in response matrices
    @param nptsK:		Power of 2 quadrature index for Kappa integral in response matrices
    @param CAUCHY:		Use Cauchy's theorem for integral evaluation
    @param RATINT:		Use rational function integrator
    @param HALO_TYPE:		Halo model type
    @param PTGRID:		value for PTGRID parameter in response calculation
    @param NRECS:		quadrature points in frequency evaluation
    @param DIVERGE:		Use power law extrapolation
    @param DIVEXPON:		Power law exponent for extrapolation
    @param OLD:			Old response type calculation
    @param VERBOSE:		Verbose diagnostic output (default: 0)
    @param SITYPE:		Non-singular integration method (default: 2=rombint)
    @param RMODMAX:		Maximum radius for halo model
    @param ROFFSET:		Radial offset for EmbededDisk
    @param MODRADIUS:		Model radius for TwoPowerSphere
    @param MODMASS:		Model mass for TwoPowerSphere
    @param AEXPON:		Core power law for TwoPowerSphere
    @param BEXPON:		Envelope power law for TwoPowerSphere
    @param DELTA:		Delta value for halo response computation
    @param HALO_MODEL:		Halo model index (0=file)
    @param OMPI:		Constant complex frequency offset
    @param NUMDF:		Number of points for Eddington inversion
    @param NICE:		Set process priority
    @param RA:			Anistropic radius for Oskipkov-Merritt inversion
    @param INCLINE:		Euler angles in degrees (2nd rotation about Y axis)
    @param PSI:			Euler angles in degrees (3rd rotation about Z axis)
    @param PHIP:		Euler angles in degrees (1st rotation about Z axis)
    @param TYPE:		Expansion basis type (0=bessel, 3=Sturm)
    @param MODEL:		Halo model type (0=input file)
    @param NUMT:		Number of times for wake computatoin (e.g. for movies)
    @param DIVERGE:		Use power law extrapolation
    @param DIVEXPON:		Power law exponent for extrapolation
    @param Rfac:		Radius for evaluation of characteristic frequency
    @param Mfac:		Multiplicative factor for maximum frequency
    @param E:			Orbital energy
    @param Rperi:		Pericenter value
    @param Rsoft:		Satellite softening radius
    @param Tmax:		Maximum time
    @param delT:		Orbital time step
    @param Toffset:		Offset from pericenter (which is T=0 by default)
    @param MASS:		Satellite mass
    @param logL:		Log Lambda value for satellite orbit decay
    @param INFILE:		Model file
    @param CACHEDIR:		Directory containing cache files
    @param UseCache:		Attempt to read and write cache
    @param Circ:		Interpret Rperi as the circular orbit radius
*/
class UserSatWake : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  //============================================================
  // Parameters
  //============================================================

  int Number, Pmax, LMIN, LMAX, MMIN, MMAX, lmax, nmax, HALO_TRUNC;
  int nptsE, nptsK, CAUCHY, RATINT, HALO_TYPE, PTGRID, NRECS, nfreqs;
  int OLD, VERBOSE, HALO_MODEL, NUMDF, numr, nint;
  int TYPE, MODEL, NUMT, DIVERGE, NUMXY;
  double RMODMAX, RA, INCLINE, PSI, PHIP, E, Rperi, Rsoft, Tmax, delT;
  double DIVEXPON, DELTA, OMPI, satmass, Toffset, Rfac, Mfac, logL;
  double rmin, rmax, scale, XYMAX;
  string INFILE, CACHEDIR, ctr_name;
  IntegrationType SITYPE;
  bool UseCache, RespChk, Circ;

  //============================================================

  vector<Eigen::VectorXcd> curcoefs;
  vector< vector<Eigen::VectorXcd> > rcoefs;
  vector<double> Times;
  vector<int> Lhalo, Mhalo;

  AxiSymModPtr halo_model;
  AxiSymBioPtr u;
  Component    *c0;
  
  void userinfo();
  void initialize_coefficients();
  void compute_coefficients();
  void check_response();
  void gnuplot_out(std::vector<Eigen::VectorXcd>& coefs, 
		   RespMat::gravity grav, RespMat::response type,
		   bool relative, string& file);
  std::complex<double> I;
  double rtol;
  int Nhalo;
  string cachefile;

public:
  

  //! Constructor
  UserSatWake(const YAML::Node& conf);

  //! Destructor
  ~UserSatWake();
  
};

#endif
