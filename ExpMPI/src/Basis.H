#ifndef _Basis_H
#define _Basis_H

#include <PotAccel.H>
#include <Vector.h>

//! Defines a basis-based potential and acceleration class
class Basis : public PotAccel
{

protected:

public:

  //! Constructor
  Basis(string& line);

  //! All derived classes must provide a method to compute the coefficients
  virtual void determine_coefficients(void) = 0;

  //! All derived classes must provide a method to compute the force
  virtual void determine_acceleration_and_potential(void) = 0;
  
  

  //! Save coefficients to file (need type marker to id dump, component id?)
  virtual void dump_coefs(ostream& out) = 0;

  //! Utility
  // @{

  //! Compute Associated Legendre Polynomials, return Matrix type
  void legendre_R(int lmax, double x, Matrix& p);
  /** Compute Associated Legendre Polynomials and derivitives, 
      return Matrix type */
  void dlegendre_R(int lmax, double x, Matrix &p, Matrix &dp);
  //! Compute vectors of sines and cosines by recursion
  void sinecosine_R(int mmax, double phi, Vector& c, Vector& s);

  // @}

};


#endif // _Basis_H


