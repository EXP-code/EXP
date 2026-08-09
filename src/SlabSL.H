#ifndef _SlabSL_H
#define _SlabSL_H

#include <complex>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include "Coefficients.H"
#include "SLGridMP2.H"
#include "biorth1d.H"
#include "Basis.H"

#if HAVE_LIBCUDA==1
#include <thrust/complex.h>
#include <cudaParticle.cuH>
#include <cudaMappingConstants.cuH>
#endif

/** @class SlabSL
    @brief This routine computes the potential, acceleration and density
    using expansion periodic in X & Y and outgoing vacuum boundary
    conditions in Z

    @details **YAML configuration**

    @param nmaxx is the maximum order of the expansion in x (default 6)

    @param nmaxy is the maximum order of the expansion in y (default 6)

    @param nmaxz is the maximum order of the expansion in z (default 6)

    @param nminx is the minimum order of the expansion in x (default 0)

    @param nminy is the minimum order of the expansion in y (default 0)

    @param hslab is the scale height of the slab (default 0.2)

    @param zmax is the maximum z for the slab (default 10.0)

    @param ngrid is the number of grid points in z for the
    Sturm-Liouville solver (default 1000)

    @param cmap is the coordinate mapping for the SL solver (must be one of "tanh", "sech", or "linear", default "linear")

    @param type is the type of slab to solve for (default
    "isothermal", must be "isothermal", "parabolic", or "constant")

    @param no_even is a flag to exclude even modes in the expansion or
    K!=0 (default false)

    @param no_odd is a flag to exclude odd modes in the expansion for
    K!=0 (default false)

    @param self_consistent set to true allows the particles to evolve
    under the time-dependent basis expansion. For a basis fixed in time to
    the initial time: set to false.

    @param zero_coefs set to true will zero the coefficients at the
    start of the simulation for all k!=0 subspaces.  This is useful
    for testing the basis expansion with a fixed potential.

    @param nint integer number of steps between subsample computations
    (default: 0, no subsampling)

    @param samplesz integer number of subsample bins (default: 0,
    no subsampling). Must be greater than 0 when `nint>0`.

    @param orthochk boolean flag to compute orthogonality matrix for the
    basis (default false).  For debugging only; set to false for production.

    @param fixed0 boolean flag to set the k=0, n=0 to a background
    value with all other k=0 terms set to zero.

    @param value0 is the coefficient of the k=0, n=0 background for
    fixed0=true.

    @param cachename is the name of the basis cache file.  This is a
    required parameter.
*/
class SlabSL : public Basis
{

  //! Pull in base-class overrides
  using PotAccel::determine_coefficients;

  //! Header structure for Sturm-Liouville slab expansion
  //! Used for deprecated native coefficient files
  struct SlabSLCoefHeader {
    double time;
    double zmax;
    double h;
    int type;
    int nmaxx, nmaxy, nmaxz;
    int jmax;
  };
  
private:

  std::shared_ptr<SLGridSlab> grid;

  //! Coefficients are a 3-tensor
  using coefType = Eigen::Tensor<std::complex<double>, 3>;

  //! Current coefficient tensor
  std::vector<coefType> expccof, expccofP;

  //! Coefficient tensor for frozen potential (if self_consistent=false)
  coefType expcofF;

  int nminx, nminy, sampT=0, nint=0;
  int nmaxx, nmaxy, nmaxz;
  double zmax, hslab, value0=8.23766340;

  bool no_even = false, no_odd = false, zero_coefs = true, fixed0 = false;

  int imx, imy, imz, jmax, nnmax;
  double  dfac, last=-std::numeric_limits<double>::max();
  std::complex<double> kfac;
  std::string cachename;

  std::vector<Eigen::VectorXd> zfrc, zpot, zden;

  SlabSLCoefHeader coefheader;

#if HAVE_LIBCUDA==1
  virtual void determine_coefficients_cuda();
  virtual void determine_acceleration_cuda();
  virtual void multistep_update_cuda();

  thrust::host_vector<thrust::complex<cuFP_t>> host_coefs;
  thrust::device_vector<thrust::complex<cuFP_t>> dev_coefs;

  //! Move coefficients from host to device
  virtual void HtoD_coefs();

  //! Move coefficients from device to host
  virtual void DtoH_coefs(unsigned);

  //! Assign constants on the device
  virtual void initialize_constants();

  //! Deallocate storage
  virtual void destroy_cuda();

  //@{

  //! Working device vectors
  thrust::device_vector<cudaTextureObject_t> t_d;

  //! Host vectors
  std::vector<cudaArray_t> cuInterpArray;
  thrust::host_vector<cudaTextureObject_t> tex;

  //! Helper struct to hold device data
  struct cudaStorage
  {
    thrust::device_vector<thrust::complex<cuFP_t>> dN_coef;
    thrust::device_vector<thrust::complex<cuFP_t>> dc_coef;
    thrust::device_vector<thrust::complex<cuFP_t>> dw_coef;
    thrust::device_vector<thrust::complex<cuFP_t>> df_coef;

    thrust::device_vector<cuFP_t> u_d;
    
    std::vector<thrust::device_vector<thrust::complex<cuFP_t>>> T_samp;

    void resize_coefs(int N, int osize, int gridSize, int sampT);
  };

  //! A storage instance
  cudaStorage cuS;

  //! Only initialize once
  bool initialize_cuda_slab = true;

  //! Initialize cuda extra
  void cuda_initialize();

  //! Initialize the cuda streams
  void slab_initialize_cuda()
  {
    grid->initialize_cuda(cuInterpArray, tex);
  }

  //! Zero the coefficient output vectors
  void cuda_zero_coefs();
  //@}

#endif

  //! Flag self_consistency
  bool self_consistent = true;

  //! Flag whether coefficients have been initialized for the first time
  bool firstime_coef = true;

  //! Default number of grid points for SLGridSlab
  int ngrid = 1000;

  //! Default coordinate mapping for SLGridSlab (must be "tanh", "sech", or "linear")
  std::string cmap = "linear";

  //! Default slab type (must be "isothermal", "parabolic", or "constant")
  std::string type = "isothermal";

  //@{
  //! Covariance structures.  First index is T, second is the
  //! flattened 3-d k vector
  std::vector<Eigen::VectorXcd> meanV;
  std::vector<Eigen::MatrixXcd> covrV;

  //@{
  //! Per thread covariance structures
  std::vector<std::vector<Eigen::VectorXcd>> meanV1;
  std::vector<std::vector<Eigen::MatrixXcd>> covrV1;
  std::vector<Eigen::VectorXcd> workV1;
  std::vector<Eigen::VectorXi>  countV1;
  std::vector<Eigen::VectorXd>  massV1;
  //@}

  void init_covariance();
  void zero_covariance();
  //@}

  //@{
  //! Usual evaluation interface
  void determine_coefficients(void);
  void determine_acceleration_and_potential(void);
  void get_acceleration_and_potential(Component*);
  //@}

  //! Swap coefficients
  void swap_coefs(std::vector<coefType>& from, std::vector<coefType>& to)
  {
    if (from.size() != to.size()) {
      std::ostringstream sout;
      sout << "swap_coefs: size(from)=" << from.size() << " != "
	   << "size(to)=" << to.size();
      throw std::runtime_error(sout.str());
    }

    std::vector<coefType> tmp(from);
    from = to;
    to   = tmp;
  }

  //@{
  //! Threading
  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);
  //@}

  /** Extrapolate and sum coefficents per multistep level to get
      a complete set of coefficients for force evaluation at an
      intermediate time step
  */
  void compute_multistep_coefficients();

  //! For updating levels
  //@{
  std::vector< std::vector<coefType> > differ1;
  std::vector< std::complex<double> > pack, unpack;
  //@}

  //@{
  //! Interpolation arrays
  using coefTypePtr = std::shared_ptr<coefType>;
  std::vector<coefTypePtr> expccofN;
  std::vector<coefTypePtr> expccofL;
  //@}

  /** Update the multi time step coefficient table when moving particle 
      <code>i</code> from level <code>cur</code> to level 
      <code>next</code>
  */
  //@{
  virtual void multistep_update_begin();
  virtual void multistep_update(int cur, int next, Component* c, int i, int id);
  virtual void multistep_update_finish();
  //@}

  //! Coefficient container instance for writing HDF5
  CoefClasses::SlabCoefs slabCoefs;

  // Biorth ID
  static const int ID=1;

  //! Write basis-specific parameters to HDF5 covariance file
  void writeCovarH5Params(HighFive::File& file);

protected:

  //! Parse parameters and initialize on first call
  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

  //! Print orthogonality matrix for debugging
  bool orthochk = false;

public:

  //! Id string
  std::string id;

  //! Constructor
  SlabSL(Component* c0, const YAML::Node& conf);

  //! Destructor
  virtual ~SlabSL();

  //! Coefficient output
  void dump_coefs_h5(const std::string& file);

  /** Return a subsample of the basis coefficients and covariance
      elements for matrices for analysis. The default implementation
      returns empty vectors. This parallels the implementation in
      expui. */
  CovarData getSubsample();

  //! Indicate whether or not subsampling is supported
  virtual bool hasSubsample() {	return true; }

  //! Return density, potential, and potential gradient at a point in
  //! Cartesian coordinates
  void determine_fields_at_point(double x, double y, double z,
				 double *tdens0, double *tpotl0, 
				 double *tdens, double *tpotl, 
				 double *tpotx, double *tpoty, 
				 double *tpotz);

  //! Returns the potential, its derivatives & density at a point in
  //! spherical coordinates
  void determine_fields_at_point_sph
  (double r, double theta, double phi,
   double *tdens0, double *dpotl0, double *tdens, double *tpotl, 
   double *tpotr, double *tpott, double *tpotp) { 
    throw std::runtime_error("SlabSL::determine_fields_at_point_sph: Not implemented for this basis");
  }

  //! Returns the potential, its derivatives & density at a point in
  //! cylindrical coordinates
  void determine_fields_at_point_cyl
  (double r, double z, double phi,
   double *tdens0, double *dpotl0, double *tdens, double *tpotl,
   double *tpotr, double *tpotz, double *tpotp) {
    throw std::runtime_error("SlabSL::determine_fields_at_point_cyl: Not implemented for this basis");
  }
};


#endif
