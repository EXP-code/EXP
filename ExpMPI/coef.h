
#ifndef _coef_h
#define _coef_h 1

struct CoefHeader {
  double time;
  int mmax;
  int nord;
  int nmax;
};

#ifdef IS_MAIN

CoefHeader coefheader;

#else

extern CoefHeader coefheader;

#endif /* IS_MAIN */

#endif /* _coef_h */
