#ifndef _SlabSL_H
#define _SlabSL_H

#include <kevin_complex.h>
#include <SLGridMP2.h>
#include <biorth1d.h>

#include <PotAccel.H>

/*! This routine computes the potential, acceleration and density
  using expansion periodic in X & Y and outgoing vacuum boundary
  condtions in Z */
class SlabSL : public PotAccel
{

//! Header structure for Sturm-Liouville slab expansion
struct SlabSLCoefHeader {
  double time;
  double zmax;
  double h;
  int type;
  int nmaxx, nmaxy, nmaxz;
  int jmax;
};

private:

  SLGridSlab *grid;

  KComplex **expccof;
  double *expreal, *expimag;
  double *expreal1, *expimag1;
  OneDTrig **trig;

  int nminx, nminy;
  int nmaxx, nmaxy, nmaxz;
  double zmax, hslab;

  int imx, imy, imz, jmax, nnmax;
  double  dfac;
  KComplex kfac;

  Vector *zfrc, *zpot;

  SlabSLCoefHeader coefheader;

  int NGRID;

  // Usual evaluation interface

  void determine_coefficients(void);
  void get_acceleration_and_potential(Component*);

  // Threading

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

  // Biorth ID
  static const int ID=1;

protected:

  //! Parse parameters and initialize on first call
  void initialize(void);

public:

  //! Id string
  string id;

  //! Constructor
  SlabSL(string& line);

  //! Destructor
  virtual ~SlabSL();

  //! Print coefficients to output stream
  void dump_coefs(ostream& out);
};


#endif
