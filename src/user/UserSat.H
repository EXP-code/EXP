#ifndef _UserSat_H
#define _UserSat_H

#include <memory>
#include <cmath>

#include "expand.H"
#include <localmpi.H>
#include <UnboundOrbit.H>
#include <LinearOrbit.H>
#include <SatelliteOrbit.H>

#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <FindOrb.H>

class UserSat : public ExternalForce
{
private:
  
  string com_name, orbfile;
  Component *c0;

  YAML::Node config;

  std::shared_ptr<Trajectory> traj;

  //! Overload to switch between cuda and cpu threaded implementation
  void determine_acceleration_and_potential();

#if HAVE_LIBCUDA==1
  //! Cuda implementation
  void determine_acceration_and_potential_cuda();
#endif
  //! CPU threaded implementation
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double core, mass, ton, toff, delta, toffset;

  bool orbit;			//! Print out orbit trajectory
  bool shadow;			//! Add a reflection-symmetric satellite
  bool verbose;			//! Print diagnostic messages
  bool zbflag;			//! Zero body flag
  bool pinning;			//! Pin center to a component

  double omega, phase, r0, tlast;

  void userinfo();

  enum TrajType {circ, bound, unbound, linear} traj_type;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:
  
  //! For debugging . . .
  static int instances;

  //! Constructor
  UserSat(const YAML::Node& conf);

  //! Destuctor
  ~UserSat();

};

#endif

