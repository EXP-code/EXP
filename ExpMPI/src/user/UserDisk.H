#ifndef _UserDisk_H
#define _UserDisk_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Bar rotating at specified fraction of corotation with a monopole

    @param a is the disk scale length
    @param mass is the total mass of the disk
    @param Ton is the time at the center of the <code>erf</code> turn on
    @param Toff is the time at the center of the <code>erf</code> turn off
    @param DeltaT is the spread of the turn on
    @param ctrname defines the component that defines the center (otherwise the center is fixed to the origin)
    @param Rmax is the maximum cylindrical radius for the potential/force table
    @param Zmax is the maximum height for the potential/force table
    @param Rnum is the number of radial table points
    @param Znum is the number of vertical table points
    @param Nint is the number of integration knots
*/
class UserDisk : public ExternalForce
{
private:
  
  string ctr_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void determine_acceleration_and_potential(vector<Particle>*);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double a, mass, Ton, Toff, DeltaT;

  double *Rtable, *Ptable, *Ztable;
  double Rmax, Zmax, dR, dZ, dfac;
  int Rnum, Znum, Nint;
  bool debug;

  string name;

  void genTable();
  void printTable();
  void getTable(double R, double Z,
		double &pot, double &fr, double &fz);

  void userinfo();

public:

  //! Constructor
  UserDisk(string &line);

  //! Destructor
  ~UserDisk();

};

#endif
