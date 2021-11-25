#ifndef _UserResPotOrb_H
#define _UserResPotOrb_H

#include <massmodel.H>
#include <biorth.H>
#include <ResPotOrb.H>
#include <TwoBodyDiffuse.H>

/** Apply specific component of resonance potential and optionally,
    the spherical restoring force for a satellite orbit specified a
    priori.  At each step, the position (and the pattern speed) is
    computed by table lookup; the orbit decay is NOT self consistent.

    @param LMAX is the maximum harmonic order
    @param NMAX is the maximum radial order
    @param NUMR is the radial grid for the biorthogonal expansion
    @param ALPHA is the power law index for J grid point distribution
    @param L0 is the "l" index for the resonance
    @param M is the "m" index for the resonance
    @param L1 is the "l_1" index for the resonance
    @param L2 is the "l_2" index for the resonance
    @param ton is the time at the center of the <code>erf</code> turn on
    @param toff is the time at the center of the <code>erf</code> turn off
    @param delta is the spread of the turn on
    @param ITMAX is the maximum number of iterations for mapping solution
    @param DELE is the fractional offset in the energy grid
    @param DELK is the offset in the kappa grid
    @param DELB is the offset in the beta grid
    @param MASS of the satellite
    @param AMP is a mass prefactor
    @param NUME is the grid in energy space for the force evaluation
    @param NUMK is the grid in kappa space for the force evaluation
    @param pmass is the mass for two-body diffusion (<=0 is off)
    @param usetag is the index of the use flag (if defined)
    @param use_background is 1 (0) if spherical potential is applied (not)
    @param model is the name of the file containing the model
    @param data is the name of the file containing the orbit
    @param cgtrname is the name of component to whose center this force will be pinned

    The phase and pattern speed are determined from the data file.
    The data file should contain the phase space of the orbit as a function
    of time with each record having the following 7 fields:

    'time' 'x-pos' 'y-pos' 'z-pos' 'x-vel' 'y-vel' 'z-vel'

*/
class UserResPotOrb : public ExternalForce
{
private:
  
  string com_name;
  Component *c0;

  AxiSymModPtr halo_model;
  std::vector<std::shared_ptr<ResPotOrb>> respot;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double ton, toff, delta, toffset, phase0, omega0, omega, Iz;
  
  int LMAX, NMAX, NUMR, L0, M0, usetag;
  double ALPHA, scale, Klim;
  std::vector<int> L1, L2;
  int numRes;
  
  double MASS, MFRAC, LENGTH, AMP, COROT, A21, A32;

  int NUMX, NUME, RECS;
  int ITMAX;
  double DELE, DELK, DELB;

  bool first;
  double phase, tlast, omlast;

  bool self;
  double domega, tom0, dtom;

  bool debug;

  std::vector<int> btotn;
  std::vector< vector<int> > bcount;
  std::vector< vector<double> > difLz;
  std::vector<double> difLz0, difLNP;

  std::string ctr_name;
  std::string model_file;
  std::string data_file;
  std::string filename;
  std::string fileomega;

  std::vector<double> Time, Radius, Phase, Omega;
  double get_radius(double);
  double get_phase(double);
  double get_omega(double);

  double pmass;
  std::shared_ptr<TwoBodyDiffuse> diffuse;

  void userinfo();

public:

  //! Constructor
  UserResPotOrb(const YAML::Node& conf);

  //! Destructor
  ~UserResPotOrb();

};


#endif
