#ifndef _SlabSL_H
#define _SlabSL_H

#include <complex>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include "Coefficients.H"
#include "SLGridMP2.H"
#include "biorth1d.H"
#include "PotAccel.H"

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

    @param type is the type of slab to solve for (default
    "isothermal", must be "isothermal", "parabolic", or "constant")

    @param self_consistent set to true allows the particles to evolve
    under the time-dependent basis expansion. For a basis fixed in time to
    the initial time: set to false.

*/
class SlabSL : public PotAccel
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

  int nminx, nminy;
  int nmaxx, nmaxy, nmaxz;
  double zmax, hslab;

  int imx, imy, imz, jmax, nnmax;
  double  dfac;
  std::complex<double> kfac;

  std::vector<Eigen::VectorXd> zfrc, zpot;

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

    void resize_coefs(int N, int osize, int gridSize, int stride);
  };

  //! A storage instance
  cudaStorage cuS;

  //! Only initialize once
  bool initialize_cuda_slab = true;

  //! Initialize cuda extra
  void cuda_initialize();

  //! Initialize the cuda streams
  void initialize_cuda()
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

  //! Default slab type (must be "isothermal", "parabolic", or "constant")
  std::string type = "isothermal";

  //@{
  //! Usual evaluation interface
  void determine_coefficients(void);
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

protected:

  //! Parse parameters and initialize on first call
  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Id string
  std::string id;

  //! Constructor
  SlabSL(Component* c0, const YAML::Node& conf);

  //! Destructor
  virtual ~SlabSL();

  //! Coefficient output
  void dump_coefs_h5(const std::string& file);
};


#endif
