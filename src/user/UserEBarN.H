#ifndef _UserEBarN_H
#define _UserEBarN_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <ZBrent.H>

/** Bar with density stratified on confocal ellipsoidal shells.

    @param length is the bar length
    @param bratio is the ratio of the submajor axis to the major axis
    @param cratio is the ratio of the minor axis to the submajor axis
    @param amp is the bar amplitude in fraction of the monopole mass
    @param barmass is the total mass of bar (monopole part)
    @param Ton is the time at the center of the <code>erf</code> turn on for the quadrupole
    @param Toff is the time at the center of the <code>erf</code> turn off for the quadrupole
    @param DeltaT is the spread of the turn on
    @param TmonoOn is the time at the center of the <code>erf</code> turn on
    for the monopole
    @param TmonoOff is the time at the center of the <code>erf</code> turn off
    for the quadrupole
    @param DeltaMonoT is the spread of the turn on for the quadrupole
    @param dtom is the width of forced bar slow down (>0 for erf)
    @param Fcorot is fraction of the length relative to the corotation radius
    @param fixed if nonzero, pattern speed is computed at every step from the current potential, if zero, pattern speed is evaluated once to start (default: false).
    @param self if the opposite of <code>fixed</code> included here for consistency with other routines
    @param monopole force is self-consistently computed from Newton's third law (default: true)
    @param follow computes monopole center self-consitently (default: true)
    @param onoff to apply turn-on, turn-off to monopole (default: false)
    @param monofrac is the fraction of the monopole to turn off (default: 1.0)
    @param quadfrac is the fraction of the quadrupole to turn off (default: 1.0)
    @param alpha is the power hardness (softness) of the bar quadrupole (default: 5)
    @param angmomfac factor artificially changing the self-consistent total bar angular momentum (default: 1.0)
    @param ctrname defines the component that defines the center (defeats inertial centering)
    @param angmname defines the component that defines the L_z reservoir (null name: L_z[init]=0)
    @param tblname sets and defines the name of a file containing a list of <code>time, b1, b5</code> values (off by default)
    @param bartype sets the profile type (see BarType)
    @param model parameter is the exponent of a power law density profile for if <code>bartype=powerlaw</code>, the Ferrers index if <code>bartype=ferrers</code> or the exponential scale length if <code>bartype=expon</code>.
    @param rmin is the minimum radius for the quadrupole table
    @param rmax is the maximum radius for the quadrupole table
    @param numt is the number of entries for the table

    The potential is evaluated using the tecniques described in
    Chandrasekhar (1969) EFE as follows.  The density is assumed to be
    stratified on confocal ellipsoids:
    \f[
    \rho({\bf r}) = \rho(m^2)
    \f]
    where
    \f[
    m^2 = \sum_{i=1}^3 {x_i^2\over a_i^2}.
    \f]
    Note that \f$m\f$ is analogous to the radius and is proptional to the true 
    radius for a spherical distribution (\f$a_1=a_2=a_3\f$).
    Inside the ellipsoid, the gravitational potential is given by:
    \f[
    \Phi_{in}({\bf r}) = -G\pi a_1a_2a_3 \int_0^\infty {du\over\Delta} 
    \int^\infty_{m^2(u)} dm^2 \rho(m^2)
    \f]
    where 
    \f[
    m^2(u) = \sum_{i=1}^3 {x_i^2\over a_i^2 + u}
    \f]
    and
    \f[
    \Delta^2 = \prod_{i=1}^3(a_i^2 + u).
    \f]
    Outside the ellipsoid, the gravitational potential is given by:
    \f[
    \Phi_{out}({\bf r}) = -G\pi a_1a_2a_3 \int_\lambda^\infty {du\over\Delta} 
    \int^\infty_{m^2(u)} dm^2 \rho(m^2)
    \f]
    where 
    \f[
    1 = \sum_{i=1}^3 {x_i^2\over a_i^2 + \lambda}.
    \f]

    Two density profiles are included.  A power-law profile of the form:
    \f[
    \rho = \rho_o m^{2\alpha}
    \f]
    where \f$\alpha>-1\f$ and an exponential profile
    \f[
    \rho = \rho_o {e^{-a_1 m/ h} \over m}
    \f]
    where \f$h\f$ is the disk scale length.  The latter profile should 
    approximate an exponential disk in projection when \f$a_1, a_2 \gg a_3\f$.


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
class UserEBarN : public ExternalForce
{
private:
  
  enum BarType {powerlaw, ferrers, expon} bartype;

  string ctr_name, angm_name;
  Component *c0, *c1;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  std::vector<double> teval;
  std::vector< std::vector<double> > tacc;
  double bps[3], vel[3], acc[3], acc1[3];

  double length, bratio, cratio, amplitude, barmass, Ton, Toff, DeltaT, Fcorot;
  double TmonoOn, TmonoOff, DeltaMonoT, modelp, monoamp;
  double monopole_frac, quadrupole_frac, alpha, model_param, rsmin, rsmax;
  bool fixed, monopole, monopole_follow, monopole_onoff, table;
  int numt, N;
  string filename;

  static constexpr double numfac = 3.86274202023190e-01;
  bool firstime;
  double posang, lastomega, lasttime, angmomfac;
  double Lz, Lz0, Lzbar, omega, omega0, Iz;
  double T0, DOmega, dtom;
  string name;

  double ldr, lrmin, lrmax;
  vector<double> rr, mm, uu, pp;

  void userinfo();
  
  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

  double rho0;
  std::shared_ptr<LegeQuad> gq, gt;
  double Density(vector<double> x);
  double Potential(vector<double> x);
  double RhoBar(double r);
  double U22(double r);
  void Inertia(vector<double>& I);
  bool quadpot(double r, double& p, double& f, double& fr, double& m);

public:

  //! Constructor
  UserEBarN(const YAML::Node& conf);

  //! Destructor
  ~UserEBarN();

};

#endif
