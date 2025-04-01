#ifndef _UserBar_H
#define _UserBar_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Bar rotating at specified fraction of corotation

    @param length is the bar length
    @param bratio is the ratio of the submajor axis to the major axis
    @param cratio is the ratio of the minor axis to the submajor axis
    @param amp is the bar amplitude in fraction of the monopole mass
    @param Ton is the time at the center of the <code>erf</code> turn on
    @param Toff is the time at the center of the <code>erf</code> turn off
    @param DeltaT is the spread of the turn on
    @param Fcorot is fraction of the length relative to the corotation radius
    @param fixed if nonzero, pattern speed is computed at every step from the current potential, if zero, pattern speed is evaluated once to start (default: false).
    @param soft use the shallow form of the potential (default: false)
    @param ctr_name defines the component that defines the bar center
    @param angm_name defines the component that defines the L_z reservoir
*/
class UserBar : public ExternalForce
{
private:
  
  std::string angm_name, ctr_name;
  Component *c0, *c1;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

#if HAVE_LIBCUDA==1
  //! Cuda implementation
  void determine_acceleration_and_potential_cuda();
#endif

  double length, bratio, cratio, amplitude, Ton, Toff, DeltaT, Fcorot;
  bool fixed, soft;
  string filename;

  static constexpr double numfac = 3.86274202023190e-01;
  bool firstime;
  double posang, lastomega, lasttime;
  double Lz, Lz0, Lzbar, omega, Iz, afac, b5;
  string name;

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserBar(const YAML::Node& conf);

  //! Destructor
  ~UserBar();

};

#endif
