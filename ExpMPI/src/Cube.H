#ifndef _Cube_H
#define _Cube_H

#include <kevin_complex.h>

#include <PotAccel.H>

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
  int nmaxx, nmaxy, nmaxz;
  int nminx, nminy, nminz;

public:

  //! Id string
  string id;

  Cube(string& line);
  ~Cube();

  void determine_coefficients(void);
  void get_acceleration_and_potential(vector<Particle>*);

  // Threading

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

};


#endif
