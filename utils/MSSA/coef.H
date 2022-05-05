#ifndef _coef_H
#define _coef_H 1

//! Header structure for cylindrical expansion
struct CylCoefHeader
{
  //! Time for this dump
  double time;
  //! Maximum polar harmonic order
  int mmax;		
  //! Radial order
  int nmax;
};

//! Header structure for spherical expansion
struct SphCoefHeader
{
  char   id[64];
  double tnow;
  double scale;
  int    nmax;
  int    Lmax;
};

#endif /* _coef_H */
