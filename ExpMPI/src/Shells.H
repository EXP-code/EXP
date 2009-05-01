
#ifndef _Shells_H
#define _Shells_H

/**
   Computes the potential, acceleration and density using spherical,
   Lagrangian shells

   Settable parameters:
   @param nsample is the number of particles per sample for the shell 
   boundary computation.

*/

/* provide an extended spherical model for point mass */
#include <AxisymmetricBasis.H>
#include <massmodel.h>

class Shells : public PotAccel
{
private:

  int nsample;
  vector<double> *rgridT, *mgridT;
  vector<double> rgrid, mgrid, pgrid;
  vector<double> rgrid1, mgrid1;

  vector<int> pnumber;
  vector<double> model;
  vector< pair<double, double> > grid;

  void initialize();

  void create_grid();

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

};

#endif

