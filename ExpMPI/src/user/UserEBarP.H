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
    @param Ton is the time at the center of the <code>erf</code> turn on
    @param Toff is the time at the center of the <code>erf</code> turn off
    @param DeltaT is the spread of the turn on
    @param Fcorot is fraction of the length relative to the corotation radius
    @param fixed if nonzero, pattern speed is computed at every step from the current potential, if zero, pattern speed is evaluated once to start (default: false).
    @param monopole force is self-consistently computed from Newton's third law (default: true)
    @param soft use the shallow form of the potential (default: false)
    @param ctrname defines the component that defines the center (defeats inertial centering)
    @param angmname defines the component that defines the L_z reservoir (null name: L_z[init]=0)
    @param tblname sets and defines the name of a file containing a list of <code>time, b1, b5</code> values (off by default)
*/
class UserEBarP : public ExternalForce
{
private:
  
  string ctr_name, angm_name, table_name;
  Component *c0, *c1;
  EllipForce *ellip;

  void determine_acceleration_and_potential(void);
  void determine_acceleration_and_potential(vector<Particle>*);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double teval;
  double bps[3], vel[3], acc[3], acc1[3], **tacc;

  double length, bratio, cratio, amplitude, barmass, Ton, Toff, DeltaT, Fcorot;
  bool fixed, omegatab, monopole, soft, table, noself;
  string filename, fileomega;

  static const double numfac = 3.86274202023190e-01;
  bool firstime;
  double posang, lastomega, lasttime;
  double Lz, Lz0, Lzbar, omega, Iz, afac, b5;
  vector<double> timeq, ampq, b5q;
  vector<double> Time, Omega;
  int qlast;
  string name;

  double get_omega(double);
  void userinfo();

public:

  //! Constructor
  UserEBarP(string &line);

  //! Destructor
  ~UserEBarP();

};

#endif
