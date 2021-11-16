#ifndef _UserResPoT_H
#define _UserResPoT_H

#include <massmodel.H>
#include <biorth.H>

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
    @param NUME is the grid in energy space for the force evaluation
    @param NUMK is the grid in kappa space for the force evaluation
    @param model is the name of the file containing the model
    @param ctrname is the name of component to whose center this force will be pinned
*/
class UserRPtest : public ExternalForce
{
private:
  
  string com_name;
  Component *c0;

  Eigen::VectorXcd bcoef;
  AxiSymModPtr halo_model;
  AxiSymBioPtr halo_ortho;
  std::shared_ptr<ResPot> respot;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  int LMAX, NMAX, NUMR, L0, M0, L1, L2;
  double rmin, rmax, scale, drfac;
  
  int NUMX, NUME, RECS;
  
  bool with_ps;
  unsigned npart;

  bool first;
  double tlast;

  string ctr_name;
  string model_file;
  string filename;

  void userinfo();

public:

  //! Constructor
  UserRPtest(const YAML::Node& conf);

  //! Destructor
  ~UserRPtest();

};


#endif
