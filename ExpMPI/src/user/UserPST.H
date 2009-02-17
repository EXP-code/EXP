#ifndef _UserPST_H
#define _UserPST_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <CylindricalDisk.H>
#include <EllipsoidForce.H>

/** Piner-Stone-Teuben model for a fixed rotating bar potential for gas forcing

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


    The logic for pattern speed determination is as follows:
    <ul>
    <li> By default the pattern speed is determined from the bar length
    <code>length</code>, <code>Fcorot</code> using the background potential.
    If <code>omega</code>>0, then this value is used rather than the default
    one.
    <li> If <code>fixed</code> is <code>false</code>, then the pattern speed
    changes are determined self-consistently from the torque
    <li> if <code>self</code> is <code>true</code> and <code>dtom</code>>0,
    then the pattern speed is determined as follows:

    \f[
    \Omega = \Omega_0\left[1 + {\Delta\Omega\over 2}\left\{
    1+\hbox{erf}\left({t_{now}-T_0\over\delta t_{om}}\right)\right\}\right]
    \f]

    <li> If <code>self</code> is <code>false</code>, then the pattern speed
    is determined as follows:

    \f[
    \Omega = \Omega_0\left[1 + \Delta\Omega(t_{now}-T_0)\right]
    \f]

    </ul>
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

  double arat;		// Maximum to minimum axis ratio for bar
  double Qm  ;		// Quadrupole strength
  double rL  ;		// Corotation radius
  double rhoC;		// Central density
  double nu  ;		// Ferrer's index

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
