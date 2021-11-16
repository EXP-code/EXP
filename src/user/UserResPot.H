#ifndef _UserResPoT_H
#define _UserResPoT_H

#include <massmodel.H>
#include <biorth.H>
#include <ResPot.H>

/** Apply specific component of resonance potential and
    optionally, the spherical restoring force

    @param LMAX is the maximum harmonic order
    @param NMAX is the maximum radial order
    @param NUMR is the radial grid for the biorthogonal expansion
    @param L0 is the "l" index for the resonance
    @param M is the "m" index for the resonance
    @param L1 is the "l_1" index for the resonance
    @param L2 is the "l_2" index for the resonance
    @param rmin is the minimum radius for the biorthogonal expansion
    @param rmax is the maximum radius for the biorthogonal expansion
    @param scale for coordinate mapping in expansion
    @param drfac is the relative radial offset for numerical derivatives
    @param ton is the time at the center of the <code>erf</code> turn on
    @param toff is the time at the center of the <code>erf</code> turn off
    @param delta is the spread of the turn on
    @param toffset is the origin shift for the turn on
    @param domega is the half width of any frequency change
    @param t0 is the central point of any frequency change interval
    @param delt is the characteristic width of the frequency change interval
    @param phase0 is the initial phase
    @param domega slow down rate if <code>self</code> is false
    @param tom0 is the midpoint of forced bar slow down
    @param dtom is the width of forced bar slow down (>0 for erf)
    @param ITMAX is the maximum number of iterations for mapping solution
    @param MASS of the bar (postive to use)
    @param MFRAC is the mass fraction of the bar inside LENGTH (postive)
    @param AMP is a mass prefactor (any sign)
    @param LENGTH of the length
    @param A21 is the ratio of the submajor axis to the major axis
    @param A32 is the ratio of the minor axis to the submajor axis
    @param NUME is the grid in energy space for the force evaluation
    @param NUMK is the grid in kappa space for the force evaluation
    @param usetag is the index of the use flag (if defined)
    @param use_background is 1 (0) if spherical potential is applied (not)
    @param model is the name of the file containing the model
    @param cgtrname is the name of component to whose center this force will be pinned

    The logic for pattern speed determination is as follows:
    <ul>
    <li> By default the pattern speed is determined from the bar length 
    <code>LENGTH</code>, <code>Fcorot</code> using the background potential
    <li> If <code>self</code> is <code>true</code>, then the pattern speed 
    changes are determined self-consistently from the torque
    <li> if <code>self</code> is <code>false</code> and <code>dtom</code>>0, 
    then the pattern speed is determined as follows:
    \f[
    \Omega = \Omega_0\left[1 + {\Delta\Omega\over2}
    \left\{ 1 + \hbox{erf}\left({t_{now}-t_{om0}\over \delta t_{om}}\right)
    \right\}\right]
    \f]
    <li> If <code>self</code> is <code>false</code> and <code>dtom</code><0, 
    then the pattern speed changes at the rate:
    \f[
    {d\Omega\over dt} = \Omega_0\Delta\Omega\times dT \times
    {1\over 2}\left[1+\hbox{erf}\left({t_{now}-t_{om0}\over\delta}\right)\right]
    \times
    {1\over2}\left[1+\hbox{erf}\left({t_{off}-t_{now}\over\delta}\right)\right] 
    \f]
    </ul>
*/
class UserResPot : public ExternalForce
{
private:
  
  string com_name;
  Component *c0;

  std::shared_ptr<AxiSymModel>  halo_model;
  std::shared_ptr<Perturbation> pert;
  std::shared_ptr<ResPot>       respot;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double ton, toff, delta, toffset, phase0, omega0, omega, Iz;
  
  int LMAX, NMAX, NUMR, L0, M0, L1, L2, usetag;
  double rmin, rmax, scale, drfac, Klim;
  
  double MASS, MFRAC, LENGTH, AMP, COROT, A21, A32;

  int NUMX, NUME, RECS;
  int ITMAX;

  bool first;
  double phase, tlast, omlast;
  double domega, tom0, dtom;

  bool usebar, useorb;

  bool debug;

  vector<int> btotn;
  vector< vector<int> > bcount;
  vector<double> difLz;
  double difLz0;

  string ctr_name;
  string model_file;
  string filename;
  string fileomega;

  vector<double> Time, Phase, Radius, Omega;
  double get_omega(double);

  void userinfo();

public:

  //! Constructor
  UserResPot(const YAML::Node& conf);

  //! Destructor
  ~UserResPot();

};


#endif
