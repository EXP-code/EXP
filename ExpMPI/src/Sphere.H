#ifndef _Sphere_H
#define _Sphere_H

#include <SLGridMP2.h>
#include <SphericalBasis.H>

/** Computes the potential, acceleration and density using the Sturm
    Liouville direct solution

    \param rmin is the minimum value in the table for the 
    radial basis functions (default is 0.001)

    \param rs is the halo scale length (default is 0.067*rmax)

    \param numr is the number radial coordinate knots in the table 
    (default is 1000)
*/
class Sphere : public SphericalBasis 
{

private:

  SLGridSph* ortho;

  void initialize(void);

  void get_dpotl(int lmax, int nmax, double r, Matrix& p, Matrix& dp, int tid);

  void get_potl(int lmax, int nmax, double r, Matrix& p, int tid);

  void get_dens(int lmax, int nmax, double r, Matrix& p, int tid);

  void get_potl_dens(int lmax, int nmax, double r, Matrix& p, Matrix& d, int tid);

				// Parameters
  double rsphSL;
  double rmin;
  double rs;
  int numr;

public:
				// Global parameters
  //! Constructor
  /*! <param> line passes any parameters (rsphSL, rmin and numr) */
  Sphere(string& line);

  //! Destructor
  ~Sphere();
};

#endif


