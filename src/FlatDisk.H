#ifndef _FlatDisk_H
#define _FlatDisk_H

#include <memory>
#include <map>

#include <BiorthCyl.H>
#include <PolarBasis.H>

typedef std::shared_ptr<BiorthCyl> CylPtr;

#if HAVE_LIBCUDA==1
#include <cudaUtil.cuH>
#endif

/** Computes the potential, acceleration and density using an empirical
    function analysis for a flat Bessel basis and a target density

    \param acyltbl is the scale length of target density (default is
    0.6)

    \param rcylmin is the minimum value in the table for the 
    radial basis functions (default is 0.0)

    \param rcylmax is the maximum value in the table for the 
    radial basis functions in acyltbl units (default is 10.0)

    \param scale is the expansion factor to get the physical scale from
    the internal scale (default is 0.01)

    \param numx is the number of grid points in the scaled radial direction

    \param numy is the number of grid points in the scaled vertical direction

    \param numr is the number radial coordinate knots in the table 
    (default is 1000)

    \param nfid is the number of basis elements used in the EOF computation

    \param knots is the number of quadrature points for EmpCyl2d

    \param logr scales the EmpCyl2d grid logarithmically

    \param model name for EmpCyl2d (default: expon)

    \param biorth set for EmpCyl2d (default: bess)
*/
class FlatDisk : public PolarBasis 
{

private:

  CylPtr ortho;

  void initialize(void);

  void get_dpotl(double r, double z,
		 Eigen::MatrixXd& p, Eigen::MatrixXd& dpr, Eigen::MatrixXd& dpz, int tid);

  void get_potl(double r, double z, Eigen::MatrixXd& p, int tid);

  virtual double getRtable() { return ortho->getRtable(); }

#if HAVE_LIBCUDA==1
  virtual void initialize_cuda()
  {
    sampT = floor(sqrt(component->CurTotal()));
    ortho->initialize_cuda(cuInterpArray, tex);
  }

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    return ortho->getCudaMappingConstants();
  }
#endif

  //! Background evaluation
  virtual std::tuple<double, double, double>
  get_pot_background(double r, double z)
  {
    auto [p, dr, dz] = ortho->background(r, z);
    return {p, -dr, -dz};
  }

  void get_dens(double r, double z, Eigen::MatrixXd& p, int tid);

  void get_potl_dens(double r, double z,
		     Eigen::MatrixXd& p, Eigen::MatrixXd& d, int tid);

				// Parameters
  double acyltbl, scale, rcylmin, rcylmax;
  int    numr;
  int    mmax;
  int    knots;
  string model;
  string biorth;
  bool   logr;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

protected:

  //! Use BiorthCyl to evaluate the basis at r, z
  void get_pot(Eigen::MatrixXd& Vc, Eigen::MatrixXd& Vs, double r, double z)
  { ortho->get_pot(Vc, Vs, r, z); }

public:
				// Global parameters
  /** Constructor
      @param c0 is the instantiating caller (Component)
      @param conf passes any parameters to basis instance (rsphSL, rmin and numr
      @param m allows the spherical basis to be used for multiple center expansions 
      
      Input line parameters include:
      @param rs is the radius for coordinate scaling
      @param numr is the number of radial grid points
      @param cmap set to true for scaling coordinates from the semi-infinite to finite segment
      @param diverge set to true means assume a cuspy profile
      @param dfac is the power of the "divergent" cusp
      @param modelname is the file containing the input background model profile
      @param dtime is the interval between basis recomputations (<=0 for never)
  */
  FlatDisk(Component* c0, const YAML::Node& conf, MixtureBasis* m=0);

  //! Destructor
  virtual ~FlatDisk();
};

#endif


