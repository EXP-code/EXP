#ifndef _Sphere_H
#define _Sphere_H

#include <SLGridMP2.h>
#include <SphericalBasis.H>

#if HAVE_LIBCUDA==1
#include <cudaUtil.cuH>
#endif

/** Computes the potential, acceleration and density using the Sturm
    Liouville direct solution

    \param rmin is the minimum value in the table for the 
    radial basis functions (default is 0.001)

    \param rs is the halo scale length (default is 0.067*rmax).  This
    is used to for the coordinate mapping  (see below).

    \param cmap should be 1 for algebraic coordinate scaling, 2 for log
    coordinate scaling, 0 for no scaling as follows: 
    \f[ 
    x = {r/r_s-1 \over r/r_s+1} 
    \f]
    when <code>cmap=1</code> and 
    \f[
    x =  \log(r)
    \f]
    when <code>cmap=2</code>.  Algebraic coordinate scaling is a good
    general choice.  For spherical models with a very large dynamic
    range (e.g. deep, cuspy halos), logarithmic scaling will be a
    better choice.  Logarithmic scaling will insist that the
    <code>rmin>0</code>.

    \param numr is the number radial coordinate knots in the table 
    (default is 1000)

    \param diverge sets the power law extrapolation of a cusp profile

    \param dfac is the inverse slope of the power law
*/
class Sphere : public SphericalBasis 
{

private:

  SLGridSph* ortho;

  void initialize(void);

  void get_dpotl(int lmax, int nmax, double r, Matrix& p, Matrix& dp, int tid);

  void get_potl(int lmax, int nmax, double r, Matrix& p, int tid);

#if HAVE_LIBCUDA==1
  virtual void initialize_cuda()
  {
    cudaDeviceReset();
    ortho->initialize_cuda(cuInterpArray, tex);
  }

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    return ortho->getCudaMappingConstants();
  }

#endif

  void get_dens(int lmax, int nmax, double r, Matrix& p, int tid);

  void get_potl_dens(int lmax, int nmax, double r, Matrix& p, Matrix& d, int tid);

				// Parameters
  double rsphSL;
  double rmin;
  double rs;
  int numr;
  int cmap;
  int diverge;
  double dfac;
  string model_file;
  string cache_file;

public:
				// Global parameters
  /** Constructor
      @param line passes any parameters to basis instance (rsphSL, rmin and numr
      @param m allows the spherical basis to be used for multiple center expansions 
  */
  Sphere(string& line, MixtureBasis* m=0);

  //! Destructor
  virtual ~Sphere();
};

#endif


