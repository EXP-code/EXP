#ifndef EllipForce_H
#define EllipForce_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <interp.h>
#include <gaussQ.h>

/** Defines monopole mass and potential for homogeneous ellipsoid
    
    @param a is the major axis
    @param b is the submajor (middle) axis
    @param c is the minor axis
    @param mass is the total mass for the ellipsoid
    @param num is the number of integration knots
    @param numr is the number of grid points in the radial table
*/    
class EllipForce {
private:

  double a, b, c, mass;
  int num, numr;
  LegeQuad *lq;

  vector<double> r, m, p;

public:

  //! Constructor
  EllipForce(double a, double b, double c, double mass, 
	     int num=100, int numr=100);
  //! Destructor
  ~EllipForce();

  //! Returns the mass inside of radius x
  double getMass(double x);
  //! Returns the monopole potential at radius x
  double getPot(double x);

  //! Prints the value table (for debugging)
  void PrintTable();
  void TestTable();

};

#endif // EllipForce
