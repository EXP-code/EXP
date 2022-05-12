#ifndef _coef_H
#define _coef_H 1

//! Header structure for cylindrical expansion
struct CylCoefHeader
{
  //! Time for this dump
  double time;

  //! Number of azimuthal orders kept
  int mmax;

  //! Number of functions per azimuthal order
  int nmax;
};

//! Header structure for spherical expansion
struct SphCoefHeader
{
  char   id[64];		// Force id
  double tnow;			// Dump time
  double scale;			// Radial scaling
  int    nmax;			// Max radial order
  int    Lmax;			// Max harmonic order
};


#endif /* _coef_H */
