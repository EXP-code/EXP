#ifndef _UserEBarP_H
#define _UserEBarP_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <EllipForce.H>

/** Bar rotating at specified fraction of corotation with a monopole

    @param length is the bar length
    @param bratio is the ratio of the submajor axis to the major axis
    @param cratio is the ratio of the minor axis to the submajor axis
    @param amp is the bar amplitude in fraction of the monopole mass
    @param barmass is the total mass of bar (monopole part)
    @param angmomfac factor artificially changing the self-consistent total bar angular momentum (default: 1.0)
    @param Ton is the time at the center of the <code>erf</code> turn on for the quadrupole
    @param Toff is the time at the center of the <code>erf</code> turn off for the quadrupole
    @param DeltaT is the spread of the turn on
    @param TmonoOn is the time at the center of the <code>erf</code> turn on
    for the monopole
    @param TmonoOff is the time at the center of the <code>erf</code> turn off
    for the quadrupole
    @param DeltaMonoT is the spread of the turn on for the quadrupole
    @param monopole force is self-consistently computed from Newton's third law (default: true)
    @param onoff to apply turn-on, turn-off to monopole (default: false)
    @param monofrac is the fraction of the monopole to turn off (default: 1.0)
    @param quadfrac is the fraction of the quadrupole to turn off (default: 1.0)
    @param soft use the shallow form of the potential (default: false)
    @param ctrname defines the component that defines the center (defeats inertial centering)
    @param angmname defines the component that defines the L_z reservoir (null name: L_z[init]=0)
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
class UserEBarP : public ExternalForce
{
private:
  
  string ctr_name, angm_name, table_name;
  Component *c0, *c1;
  EllipForce *ellip;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  vector<double> teval;
  double bps[3], vel[3], acc[3], acc1[3], **tacc;

  double length, bratio, cratio, amplitude, barmass, Ton, Toff, DeltaT, Fcorot;
  double TmonoOn, TmonoOff, DeltaMonoT;
  double monopole_frac, quadrupole_frac;
  bool monopole, monopole_onoff, soft, table;
  string filename, fileomega;

  static constexpr double numfac = 3.86274202023190e-01;
  bool firstime;
  double posang, lastomega, lasttime, angmomfac;
  double Lz, Lz0, Lzbar, omega, omega0, Iz, afac, b5;
  double T0, DOmega, tom0, dtom;
  vector<double> timeq, ampq, b5q;
  vector<double> Time, Omega;
  int qlast;
  string name;

  double get_omega(double);
  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserEBarP(const YAML::Node& conf);

  //! Destructor
  ~UserEBarP();

};

#endif
