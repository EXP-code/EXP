#ifndef _UserSNheat_H
#define _UserSNheat_H

#include <random>
#include <pthread.h>

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Set off SN at a specified rate and add thermal energy to the gas
    within a specified radius of a specified origin.  The supernova
    fire off stochastically with a rate computed from the specified
    mean interval between supernova.
    @param ctrname is the name of the component whose particles will be heated
    @param delay is the time to wait after starting before setting off the 
    first SN
    @param dT is the mean spacing of SN explosion in years
    @param SNe is the energy per supernova in ergs
    @param radius is the radius for selecting particles to heat
    @param X is the origin in the x-dimension for heating
    @param Y is the origin in the y-dimension for heating
    @param Z is the origin in the z-dimension for heating
    @param Lunit is the physical scale in cm per system units (default: 300 kpc)
    @param Tunit is the physical time scale in years per system unit (default: gravitational)
    @param Munit is the physical mass in units of mp per system units (default: 1e12 Msun)
*/
class UserSNheat : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  void write_trace();

  vector<double> origin;
  string comp_name, filename;
  double delay, dT, tlast, dE, radius;
  vector< set<int> > plist;
  bool firstime, verbose;
  int N, ncount, nSN;

  std::normal_distribution<> norm;
  std::uniform_real_distribution<> unit;

  Component *c0;

  int arrivalTime(double dt);
  void userinfo();

  vector<double> ke0, ke1, mm0, mom;
  vector< vector<double> > pp0;
  double ketot0, ketot1, mass0, factor;

  static double Lunit, Munit, Tunit, Vunit, Eunit;

  pthread_barrier_t barr;
  //
  // Barrier initialization
  //
  void threading_init() 
  {
    if (pthread_barrier_init(&barr, 0, nthrds)) {
      cout << "UserSNheat: process " << myid 
	   << " could not create a barrier" << endl;
      MPI_Abort(MPI_COMM_WORLD, 36);
    }
  }

  //
  // Synchronization barrier
  //
  void Pbarrier(int id, int err) 
  {
    int rc = pthread_barrier_wait(&barr);
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {
      cerr << "UserSNheat: process " << myid << ", thread " << id
	   << " can't find fiducial component <" << comp_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, err);
    }
  }

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserSNheat(const YAML::Node& conf);

  //! Destructor
  ~UserSNheat();

};

#endif
