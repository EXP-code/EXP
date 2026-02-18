
#ifndef _Shells_H
#define _Shells_H

#include <vector>
#include <string>
#include <set>

#include "expand.H"
#include "interp.H"

/**
   Computes the potential, acceleration and density using spherical,
   Lagrangian shells

   Settable parameters:
   @param nsample is the number of particles per sample for the shell 
   boundary computation.
   @param nselect is the number of particle intervals in computing the
   cumulative distribution

*/

/* provide an extended spherical model for point mass */
#include "AxisymmetricBasis.H"
#include "massmodel.H"

class Shells : public PotAccel
{
  //! Pull in base-class overrides
  using PotAccel::determine_coefficients;

private:

  bool self_consistent, firstime_accel;
  int nsample, nselect, used1;
  vector<int> usedT;

  typedef pair<double, double> Dpair;
  vector<Dpair> grid;

  std::vector< std::vector<double> >  rgridT, mgridT;
  std::vector< std::vector<int>    >  igridT, igrid;
  std::vector<int>                    igrid1;
  vector<double>                      rgrid1, mgrid1;
  vector<double>                      rgrid0, mgrid0, pgrid0;
  std::vector<std::map<int, double>>  rgrid, mgrid;
  std::vector<std::vector<int>>       update_fr, update_to, update_ii;

  vector<int> snumbr, rnumbr, sdispl, rdispl, sndcts, rcvcts;

  void initialize();

  void determine_coefficients(void);
  void determine_acceleration_and_potential(void);

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! The constructor
  //! \param c0 is the instantiating caller (a Component)
  //! \param conf passes in any explicit parameters
  Shells(Component* c0, const YAML::Node& conf);

  //! The constructor
  virtual ~Shells();

  //! The main force call
  void get_acceleration_and_potential(Component*);

  //! Execute to begin level shifts for particles
  void multistep_update_begin();

  //! Implementation of level shifts
  void multistep_update(int cur, int next, Component* c, int i, int id);

  //! Execute to finish level shifts for particles
  void multistep_update_finish();

};

#endif

