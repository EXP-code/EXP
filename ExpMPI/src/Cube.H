#ifndef _Cube_H
#define _Cube_H

#include <kevin_complex.h>

#include <PotAccel.H>

//! Periodic cube basis
class Cube : public PotAccel
{

private:

  Complex **expccof;
  double *expreal, *expimag;
  double *expreal1, *expimag1;

  int imx, imy, imz, jmax, use1, use0;
  double  dfac;
  Complex kfac;

protected:

  void initialize(void);
				// Parameters
  //! Update me
  int nmaxx, nmaxy, nmaxz;

  //! Update me
  int nminx, nminy, nminz;

public:

  //! Id string
  string id;

  //! Constructor
  Cube(string& line);

  //! Destructor
  ~Cube();

  void determine_coefficients(void);
  void get_acceleration_and_potential(Component*);

  // Threading

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

};


#endif
