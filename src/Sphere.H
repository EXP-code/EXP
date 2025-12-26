#ifndef _Sphere_H
#define _Sphere_H

#include <memory>
#include <map>

#include "SLGridMP2.H"
#include "SphericalBasis.H"

typedef std::shared_ptr<SLGridSph> SLGridSphPtr;

#if HAVE_LIBCUDA==1
#include <cudaUtil.cuH>
#endif

/** @class Sphere
    @brief Computes the potential, acceleration and density using the Sturm
    Liouville direct solution

    @details **YAML configuration**

    @param rmapping is the halo scale length (default is 0.067*rmax).  This
    is used to for the coordinate mapping  (see below).

    @param cmap should be 1 for algebraic coordinate scaling, 2 for log
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

    @param numr is the number radial coordinate knots in the table 
    (default: 2000)

    @param nums is the number of bins in the empirical mass array
    (default: 2000)

    @param noff is the offset into the empirical mass array for
    computing upper and lower grid boundaries (default: 32)

    @param diverge set to true means assume a cuspy profile

    @param dfac is the power of the "divergent" cusp

    @param modelname is the file containing the input background model profile

    @param cachename is the name for the SL grid cache file

    @param dtime is the interval between basis recomputations (<=0 for never)

    @param diverge sets the power law extrapolation of a cusp profile

    @param dfac is the inverse slope of the power law

    @param logr set to true uses a logarithmic spacing for the
    empirical radial grid (default: false)

    @param plummer set to true uses a Plummer-model fit to the
    empirical mass array to make basis model (default: true)
*/
class Sphere : public SphericalBasis 
{

  //! Pull in base-class overrides
  using PotAccel::determine_coefficients;

private:

  SLGridSphPtr ortho;

  void initialize(void);

  void get_dpotl(int lmax, int nmax, double r, Eigen::MatrixXd& p, Eigen::MatrixXd& dp, int tid);

  void get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid);

  double mapIntrp(const std::map<double, double> &data, double x);
  double mapDeriv(const std::map<double, double> &data, double x);

  void make_model() {
    if (plummer) make_model_plummer();
    else         make_model_bin();
  }
  
  void make_model_bin();
  void make_model_plummer();

#if HAVE_LIBCUDA==1
  virtual void initialize_cuda()
  {
    ortho->initialize_cuda(cuInterpArray, tex);
  }

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    return ortho->getCudaMappingConstants();
  }

#endif

  void get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid);

  void get_potl_dens(int lmax, int nmax, double r,
		     Eigen::MatrixXd& p, Eigen::MatrixXd& d, int tid);

				// Parameters
  double rsphSL;
  double rmap;
  double tnext, dtime;
  int    numr;
  int    nums;
  int    noff;
  int    cmap;
  int    diverge;
  double dfac;
  string model_file;
  string cache_file;
  bool   recompute;
  bool   plummer;
  bool   logr;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

  //@{
  //! Coordinate mapping
  double r_to_xi(double r)  { return ortho->r_to_xi(r); }
  double xi_to_r(double x)  { return ortho->xi_to_r(x); }
  double d_r_to_xi(double r) { return 1.0/ortho->d_xi_to_r(r_to_xi(r)); }
  double d_xi_to_r(double x) { return ortho->d_xi_to_r(x); }
  //@}

public:
				// Global parameters
  /** Constructor
      @param c0 is the instantiating caller (Component)
      @param conf passes any parameters to basis instance (rsphSL, rmin and numr
      @param m allows the spherical basis to be used for multiple center expansions 
  */
  Sphere(Component* c0, const YAML::Node& conf, MixtureBasis* m=0);

  //! Destructor
  virtual ~Sphere();

  //! Override to redetermine basis, as desired
  virtual void determine_coefficients(void)
  {
    if (recompute and tnow >= tnext) make_model();
    SphericalBasis::determine_coefficients();
  }

};

#endif


