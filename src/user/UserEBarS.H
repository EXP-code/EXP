#ifndef _UserEBarS_H
#define _UserEBarS_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <EllipForce.H>

/** Bar rotating at specified fraction of corotation with a monopole
    This version computes the torque from the phase space using
    Newton's 3rd law

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
    @param tom0 is the midpoint of forced bar slow down
    @param dtom is the width of forced bar slow down (>0 for erf)
    @param Fcorot is fraction of the length relative to the corotation radius
    @param fixed if nonzero, pattern speed is computed at every step from the current potential, if zero, pattern speed is evaluated once to start (default: false).
    @param self if the opposite of <code>fixed</code> included here for consistency with other routines
    @param monopole force is self-consistently computed from Newton's third law (default: true)
    @param follow computes monopole center self-consitently (default: true)
    @param onoff to apply turn-on, turn-off to monopole (default: false)
    @param monofrac is the fraction of the monopole to turn off (default: 1.0)
    @param quadfrac is the fraction of the quadrupole to turn off (default: 1.0)
    @param soft use the shallow form of the potential (default: false)
    @param angmomfac factor artificially changing the self-consistent total bar angular momentum (default: 1.0)
    @param ctrname defines the component that defines the center (defeats inertial centering)
    @param tblname sets and defines the name of a file containing a list of <code>time, b1, b5</code> values (off by default)

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
class UserEBarS : public ExternalForce
{
private:
  
  string ctr_name, table_name;
  Component *c0;
  std::shared_ptr<EllipForce> ellip;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double tevalS;
  std::vector<double> teval;
  std::vector< std::vector<double> > tacc;
  double bps[3], vel[3], acc[3], acc1[3];
  std::vector<double> TzM, TzN;
  int okM;
  int cid;

  double length, bratio, cratio, amplitude, barmass, Ton, Toff, DeltaT, Fcorot;
  double TmonoOn, TmonoOff, DeltaMonoT;
  double monopole_frac, quadrupole_frac;
  bool fixed, monopole, monopole_follow, monopole_onoff, soft, table;
  string filename;

  static constexpr double numfac = 3.86274202023190e-01;
  bool firstime;
  double posang, lastomega, lasttime, angmomfac;
  double Lz, Tz, omega, omega0, Iz, afac, b5;
  double T0, DOmega, tom0, dtom;
  vector<double> timeq, ampq, b5q;
  int qlast;
  string name;

  void userinfo();

public:

  //! Constructor
  UserEBarS(const YAML::Node& conf);

  //! Destructor
  ~UserEBarS();

};

#endif
