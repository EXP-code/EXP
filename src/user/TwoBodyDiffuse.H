#ifndef _TwoBodyDiffuse_H
#define _TwoBodyDiffuse_H

using namespace std;

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <random>

#include <gaussQ.H>
#include <massmodel.H>

/** Random walk of each orbit on energy surface

    @param seed for the pseudorandom number generator
    @param pmass is the mass of virtual two-body particle
    @param logL for diffusion coefficients
    @param clip the Gaussian variates at this number of "Sigma"
    @param numr is the number of radial entries in diff. coefficient table
    @param numv is the number of velocity  entries in diff. coefficient table
    @param nume is the number of energy integration knots for diff. coef.
    @param logr is set to 1 for radial logarithmic spacing in diff. coef. table
    @param modfile contains the background model
    @param diverge set to extrapolate model for a power law cusp
    @param diverge_rfac exponent for the power law
*/
class TwoBodyDiffuse
{

  typedef vector<double> dvector;

private:
  
  long int seed;
  double pmass, logL, clip;
  int numr;
  int numv;
  int nume;
  bool logr;
  string modfile;
  int diverge;
  double diverge_rfac;

  SphericalModelTable *model;

  void compute_diffuse();
  void get_coefs(double r, double v,
		 double& dvpara1, double& dvpara2, double& dvperp2);


  double delR, delV, Gamma, Rmin, Rmax, Emin, Emax;
  vector<double> R, Vmax;
  vector<dvector> dVpara1, dVpara2, dVperp2;

  LegeQuad *jq;

  std::uniform_real_distribution<> urand;
  std::normal_distribution<> nrand;

public:

  //! Constructor
  TwoBodyDiffuse(double PMASS=1.0e-6,
		 bool   LOGR=true,
		 double LogL=5.7,
		 double Clip=100.0,
		 int    Numr=200,
		 int    Numv=64,
		 int    Nume=128,
		 string ModFile="SLGridSph.model",
		 int Diverge=0,
		 double Diverge_rfac=-1.0);

  //! Destructor
  ~TwoBodyDiffuse();


  //! Get velocity kick
  void get_diffusion(double dtime, double* pos, double* vel, double* dvel);

  //! Print out coefficient grid in gnuplot format
  void dump_grid(ostream* out);

};

#endif
