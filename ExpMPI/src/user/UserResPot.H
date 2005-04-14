#ifndef _UserResPoT_H
#define _UserResPoT_H

#include <massmodel.h>
#include <biorth.h>
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
    @param ITMAX is the maximum number of iterations for mapping solution
    @param MASS of the bar (postive)
    @param AMP is a mass prefactor (any sign)
    @param LENGTH of the length
    @param A21 is the ratio of the submajor axis to the major axis
    @param A32 is the ratio of the minor axis to the submajor axis
    @param NUME is the grid in energy space for the force evaluation
    @param NUMK is the grid in kappa space for the force evaluation
    @param tag is the index of the use flag (if defined)
    @param use_background is 1 (0) if spherical potential is applied (not)
    @param model is the name of the file containing the model
    @param cgtrname is the name of component to whose center this force will be pinned
*/
class UserResPot : public ExternalForce
{
private:
  
  string com_name;
  Component *c0;

  AxiSymModel  *halo_model;
  Perturbation *pert;
  ResPot *respot;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double ton, toff, delta, toffset, omega, Omega;
  
  int LMAX, NMAX, NUMR, L0, M0, L1, L2, tag;
  double rmin, rmax, scale, drfac;
  
  double MASS, LENGTH, AMP, COROT, A21, A32;

  int NUMX, NUME, RECS;
  int ITMAX;

  bool first;
  double domega, t0, phase, tlast, omlast;

  int *bcount;

  string ctr_name;
  string model_file;
  string filename;

  void userinfo();

public:

  //! Constructor
  UserResPot(string &line);

  //! Destructor
  ~UserResPot();

};


#endif
