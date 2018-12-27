#ifndef _CBrock_H
#define _CBrock_H

#include <Vector.h>
#include <SphericalBasis.H>

//! Biorthonormal pairs based on Clutton-Brock's ultraspherical series
class CBrock : public SphericalBasis
{

private:

  Vector uu, duu;

  void get_dpotl(int lmax, int nmax, double r, Matrix& p, Matrix& dp, int tid);
  void get_potl(int lmax, int nmax, double r, Matrix& p, int tid);
  void get_dens(int lmax, int nmax, double r, Matrix& p, int tid);
  void get_potl_dens(int lmax, int nmax, double r, Matrix& p, Matrix& d,
		     int tid);

  void initialize(void);

  double knl(int, int);
  double norm(int, int);

  double r_to_rh(double r);
  double rh_to_r(double rh);
  double d_r_to_rh(double);

public:

  //! Constructor
  /*! \param line is the parameter string */
  /*! \param m is the MixtureBasis for a multicenter expansion */
  CBrock(const YAML::Node& conf, MixtureBasis* m=0);

};

#endif
