#ifndef _UserDisk_H
#define _UserDisk_H

#include "Particle.H"
#include "AxisymmetricBasis.H"
#include "ExternalCollection.H"

/** Potential and force of a thin exponential disk.  The force
    is tabulated on a two-dimensional grid with positive Z and
    interpolated using bilinear interpolation.

    The default values for Nscale, Ngrid and Nint should suffice
    for most applications.

    @param a is the disk scale length
    @param mass is the total mass of the disk
    @param Ton is the time at the center of the <code>erf</code> turn on
    @param Toff is the time at the center of the <code>erf</code> turn off
    @param DeltaT is the spread of the turn on
    @param ctrname defines the component that defines the center (otherwise the center is fixed to the origin)
    @param Nscale is the physical scale of the grid in numbers of disk scale lengths
    @param Ngrid is the number of table points in each dimension
    @param Nint is the number of integration knots
*/
class UserDisk : public ExternalForce
{
private:
  
  string ctr_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double a, mass, Ton, Toff, DeltaT;

  double *Rtable, *Ptable, *Ztable;
  double Nscale, Rmax, Zmax, dR, dZ, dfac;
  int Ngrid, Nint;
  bool debug;

  string name;

  void genTable();
  void printTable();
  void getTable(double R, double Z,
		double &pot, double &fr, double &fz);

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserDisk(const YAML::Node& conf);

  //! Destructor
  ~UserDisk();

};

#endif
