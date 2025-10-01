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

/** @class FlatDisk

    @brief Computes the potential, acceleration and density using an
    empirical function analysis for a flat Bessel basis and a target
    density

    @details **YAML configuration**

    @param rcylmin is the minimum value in the table for the 
    radial basis functions (default is 0.0)

    @param rcylmax is the maximum value in the table for the 
    radial basis functions in acyltbl units (default is 10.0)

    @param scale is the expansion factor to get the physical scale from
    the internal scale (default is 0.01)

    @param mmax is the maximum azimuthal order

    @param numx is the number of grid points in the scaled radial direction

    @param numy is the number of grid points in the scaled vertical direction

    @param knots is the number of quadrature points for EmpCyl2d

    @param logr scales the EmpCyl2d grid logarithmically

    @param model name for conditioning the new basis using EmpCyl2d
    (default: expon)

    @param biorth is the biorthogonal basis set used by EmpCyl2d for
    conditioning (default: bess)

    @param background sets a fixed monopole potential and force for
    the model.  Used for stability analysis of disks with active and
    inactive tapers, following Zang

    @param nmaxfid is the maximum radial order for the basis used by
    EmpCyl2d for conditioning

    @param numr is the number of grid points in the radial direction
    for the numerical basis in EmpCyl2d

    @param NQDHT is the number of grid points for the numerical Hankel transform

    @param diskconf is a YAML configuration for the target disk model
    in EmpCyl2d

    @param cachename is the name of the cache file for the newly
    constructed basis

    @param dumpbasis provides the user with a ascii table of the basis
    potential-density pairs
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
    ortho->initialize_cuda(disk, cuInterpArray, tex);
  }

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    return ortho->getCudaMappingConstants();
  }
#endif

  //! Background instance
  using Disk2d = std::shared_ptr<EmpCyl2d::ModelCyl>;
  Disk2d disk;

  //! Set background from YAML
  void setBackground();

  //! Background evaluation
  virtual std::tuple<double, double, double>
  background(double r, double z)
  { return {disk->pot(r), -disk->dpot(r), 0.0}; }
  
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

  bool   dump_basis;

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


