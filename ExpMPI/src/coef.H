#ifndef _coef_H
#define _coef_H 1

//! Header structure for cylindrical expansion
struct CoefHeader {
  //! Time for this dump
  double time;
  //! Maximum polar harmonic order
  int mmax;		
  //! Order for each eof
  int nord;
  //! Radial order (needs updating?)
  int nmax;
};


//! Header structure for cylindrical expansion
struct CoefHeader2 {
  //! Time for this dump
  double time;
  //! Number of azimuthal orders kept
  int mmax;
  //! Number of functions per azimuthal order
  int nmax;
};

#endif /* _coef_H */
