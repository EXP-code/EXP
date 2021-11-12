#ifndef _UserDSMC_H
#define _UserDSMC_H

#include <random>

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

#include <global.H>
#include <Timer.H>

/**
   Helper class describing a comparison for sorting particles in radius
*/
class RadiusSort
{
public:
  double radius;
  int curindx;
  int curnode;

  bool operator() (const RadiusSort& x, const RadiusSort& y)
  {
    return x.radius < y.radius;
  }

  bool operator() (const double& x, const RadiusSort& y)
  {
    return x < y.radius;
  }

  bool operator() (const RadiusSort& x, const double& y)
  {
    return x.radius < y;
  }
};


/**
   Helper class describing a comparison for sorting the most bound particles
*/
class BoundSort
{
public:

  typedef pair<double, int> Bound;

  bool operator() (const Bound& x, const Bound& y)
  {
    return x.first < y.first;
  }

  bool operator() (const double& x, const Bound& y)
  {
    return x < y.first;
  }

  bool operator() (const Bound& x, const double& y)
  {
    return x.first < y;
  }
};


/** DSMC particle routine

    @param Rcloud is the effective radius of interaction (e.g. diameter of the sphere)
    @param NR is the number of radial bins
    @param NZ is the number of vertidal bins
    @param NP is the number of azimuthal bins
    @param Rmax is the maximum radius of the binned cylinders
    @param Zmax is the maximum extent of the cylinder above and below the midplane
    @param efftol is the minimum efficient allowed before triggering a repartition
    @param vfrac is the fraction of the mean velocity used to estimate the peak relative particle velocity
    @param seed is the random number seed
    @param debug may be set to true (non-zero value) to turn on debugging output
*/
class UserDSMC : public ExternalForce
{
private:
  
  string comp_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double rcloud, rmin, rmax, efftol, zmax, vfrac;
  int NR, NZ, NP, seed;

  bool debug;

				// Internal variables
				// ------------------

  int dN;			// Number of radial bins per process
  double dPhi;			// Aximuthal bin spacing

				// For determining the bin boundaries
				// 
  vector<float> radiiT, rbound;
  vector<RadiusSort> radii;
  vector< pair<double, int> > rboundI, zboundI;
  vector< pair<double, int> >::iterator iI;

				// Process assignment for each radial bin
  vector<unsigned short> owners;
				// Particle numbers per process
  vector<unsigned int> counts, countsT, nbegpts;
				// Efficiency estimate for repartitioning
  vector<unsigned int> eff, effT;
  bool regrid_flag;

				// Redistribution lists
  vector<int> redist, redistT, redistT2, nbin;

				// Bin assignments per processor
  vector< vector<int> > binlist;

				// Variables for collision computation
  vector< vector<double> > vcom, vrel, delv;
  vector<double> vmax;

				// Pseudo-random-number generation
  std::uniform_real_distribution<> unif;

  void makeGrid();
  void makeSort();
  void userinfo();

#ifdef TIMER
  Timer *timer;
#endif

public:

  //! Constructor
  UserDSMC(string &line);

  //! Destructor
  ~UserDSMC();

};

#endif
