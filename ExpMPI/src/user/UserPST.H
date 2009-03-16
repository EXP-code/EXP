#ifndef _UserPST_H
#define _UserPST_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <CylindricalDisk.H>
#include <EllipsoidForce.H>

/** Piner-Stone-Teuben model for a fixed rotating bar potential for gas forcing
    Following their paper (1995, ApJ 449,508), their model has a Kuzmin
    disk,
    \f[
    \Sigma(r) = \left(\frac{v_o^2}{2\pi G r_o}\right)
    \left(1 + \frac{r^2}{r_o^2}\right)^{-3/2},
    \f]
    where the constants \f$v_o=200\hbox{km/s}\f$ and \f$r_o=14.1\hbox{kpc}\f$,
    a bulge with the density
    \f[
    \rho(r) = \rho_o\left(1 + \frac{r^2}{r_o^2}\right)^{-3/2},
    \f]
    and a Ferrer's ellipsoid with \f$a=5\hbox{kpc}\f$ and \f$b=c\f$.

    @param rmin is the minimum model radius
    @param rmax is the maximum model radius
    @param numr is the number of radial points for the bulge model
    @param blog use log spacing for bulge model if true
    @param dlog use log spacing for disk model if true
    @param arat is the ratio of the bar model major to minor axes
    @param Qm is the bar quadrupole amplitude in (solar masses)*(kpc)^2
    @param rL is the corotation radius in kpc
    @param rhoC is the central density in (solar masses)/(kpc)^3
    @param nu is Ferrer's profile index
    @param Lmax is the maximum harmonic index for disk pot expansion
    @param Nmax is the maximum radial index for expansion
    @param numR is the number of radial table elements for the expansion
    @param numt is the number Integration knots for coefficient integration in the theta direction
    @param numg
    @param DeltaT is the spread of the turn on
    @param Ton is the time at the center of the <code>erf</code> turn on

    The slow bar turn on uses the following function:
    \f[
    \eta(t) = \frac12\left[
    1+\hbox{erf}\left({t_{now}-T_{on}\over\delta T_{on}}\right)\right].
    \f]
    The total mass of the model is held fixed during turn on by
    letting the bar mass be \f$M_{bar}\eta(t)\f$ and increasing the
    mass of the disk+bulge by \f$M_{bar}[1-\eta(t)]\f$.
*/
class UserPST : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();


  double rmin;		// Minimum model radius
  double rmax;		// Maximum model radius

  int    numr;		// Number of radial points for bulge model

  bool   blog;		// Use log spacing for bulge model
  bool   dlog;		// Use log spacing for disk model

  double arat ;		// Maximum to minimum axis ratio for bar
  double Qm   ;		// Quadrupole strength
  double rL   ;		// Corotation radius
  double rhoC ;		// Central density
  double nu   ;		// Ferrer's index

  double Mbar ;		// Implied bar mass
  double a    ;		// Implied semimajor axis
  double b    ;		// Implied semiminor axis
  double rho0 ;		// Implied central bar density
  double rhoB ;		// Implied central bulge density

  int Lmax;		// Maximum harmonic index for disk pot expansion
  int Nmax;		// Maximum radial index for expansion
  int numR;		// Number of radial table elements
  int numt;		// Integration knots for theta
  int numg;		// 


  double Ton;		// Turn on start time
  double DeltaT;	// Turn on width

  double Mscale;
  double omega;
  string filename;

  void userinfo();
  double get_fid_bulge_dens();

  SphericalModelTable *bulge;
  EllipsoidForce      *bar;
  KuzminDisk          disk;

public:

  //! Constructor
  UserPST(string &line);

  //! Destructor
  ~UserPST();

};

#endif
