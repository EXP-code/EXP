#ifndef _Slab_H
#define _Slab_H

#include <kevin_complex.h>
#include <biorth1d.h>

#include <PotAccel.H>

/*! This routine computes the potential, acceleration and density 
  using periodic box expansion in X & Y and slab in Z */
class Slab : public PotAccel
{

struct SlabCoefHeader {
  double time;
  double zmax;
  int nmaxx, nmaxy, nmaxz;
  int jmax;
};

private:

  Complex **expccof;
  double *expreal, *expimag;
  double *expreal1, *expimag1;
  OneDTrig **trig;

  int nmaxx, nmaxy, nmaxz;
  int nminx, nminy;
  double zmax;

  int imx, imy, imz, jmax, nnmax;
  double  dfac;
  Complex kfac;

  Vector *zfrc, *zpot;

  SlabCoefHeader coefheader;

  void determine_coefficients(void);
  void get_acceleration_and_potential(vector<Particle>*);

  // Threading

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

protected:

  //! Parse parameters and initialize on first call
  void initialize(void);

public:

  //! Id string
  string id;

  //! Constructor
  Slab(string& line);

  //! Destructor
  virtual ~Slab() {};

  //! Print coefficients to output stream
  void dump_coefs(ostream& out);
};


#endif
