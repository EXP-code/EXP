
#ifndef _Shells_H
#define _Shells_H

#include "expand.h"
#include "interp.h"

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
#include <AxisymmetricBasis.H>
#include <massmodel.h>

class Shells : public PotAccel
{
private:

  bool self_consistent, firstime_accel;
  int nsample, nselect, used1;
  vector<int> usedT;

  typedef pair<double, double> Dpair;
  vector<Dpair> grid;

  vector<double>   *rgridT, *mgridT;
  vector<int>      *igridT, igrid1, *igrid;
  vector<double>   rgrid1, mgrid1;
  vector<double>   rgrid0, mgrid0, pgrid0;
  map<int, double> *rgrid, *mgrid;
  vector<int>      *update_fr, *update_to, *update_ii;

  vector<int> snumbr, rnumbr, sdispl, rdispl, sndcts, rcvcts;

  void initialize();

  void determine_coefficients(void);
  void determine_acceleration_and_potential(void);

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

public:

  //! The constructor
  /** \param line passes in any explicit parameters */
  Shells(string& line);

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

