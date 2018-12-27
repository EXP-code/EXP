#ifndef _Cube_H
#define _Cube_H

#include <kevin_complex.h>

#include <PotAccel.H>

//! Periodic cube basis
class Cube : public PotAccel
{

private:

  KComplex **expccof;
  double *expreal, *expimag;
  double *expreal1, *expimag1;

  int imx, imy, imz, jmax, use1, use0;
  double  dfac;
  KComplex kfac;

protected:

  void initialize(void);
				// Parameters

  //! Minimum wave numbers in each dimesion
  //@{
  int nmaxx, nmaxy, nmaxz;
  //@}

  //! Maximum wave numbers in each dimesion
  //@{
  int nminx, nminy, nminz;
  //@}

public:

  //! Id string
  string id;

  //! Constructor
  Cube(const YAML::Node& conf);

  //! Destructor
  virtual ~Cube();

  //! Compute the coefficients
  void determine_coefficients(void);

  //! Compute the force
  void get_acceleration_and_potential(Component*);

  // Threading

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

};


#endif
