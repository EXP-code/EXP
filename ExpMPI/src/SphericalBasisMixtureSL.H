#ifndef _SphericalBasisMixtureSL_H
#define _SphericalBasisMixtureSL_H

#include <string>
#include <SphericalBasis.H>

class SphereTwoCenter;

//! Used to indicate center type for expansion
enum CenterType {com, ej};

//! Spherical basis with force evaluation mixture function
class SphericalBasisMixtureSL : public SphericalBasis
{

private:

  SphereTwoCenter* A;
  CenterType ctr;

  void get_dpotl(int lmax, int nmax, double r, Matrix& p, Matrix& dp, int tid);

  void get_potl(int lmax, int nmax, double r, Matrix& p, int tid);

  void get_dens(int lmax, int nmax, double r, Matrix& p, int tid);

  void get_potl_dens(int lmax, int nmax, double r, Matrix& p, Matrix& d, int tid);

protected:

  /**
     Thread method for coefficient accumulation 
     (gets its center from SphereTwoCenter)
  */
  virtual void * determine_coefficients_thread(void * arg);

  /** 
      Thread method for accerlation compuation
      (get its center from SphereTwoCenter)
  */
  virtual void * determine_acceleration_and_potential_thread(void * arg);

public:

  /** Constructor 
      \param line contains the paramter list
      \param A is the caller needed to access the center
      \param ctr defines the center type (<code>COM</code> or otherwise)
  */
  SphericalBasisMixtureSL(string& line, SphereTwoCenter* A, CenterType ctr);

};

#endif
