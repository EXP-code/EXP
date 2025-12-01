#ifndef _Cube_H
#define _Cube_H

#include <complex>
#include <vector>
#include <string>
#include <set>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include "Coefficients.H"
#include "PotAccel.H"

#if HAVE_LIBCUDA==1
#include <thrust/complex.h>
#include <cudaParticle.cuH>
#include <cudaMappingConstants.cuH>
#endif

/** @class Cube

    @brief Periodic cube basis

    @details **YAML configuration**

    @param nminx is the minimum absoulte value of the wave number in
    the x direction (default: 0)

    @param nminy is the minimum absoulte value of the wave number in
    the y direction (default: 0)

    @param nminz is the minimum absoulte value of the wave number in
    the z direction (default: 0)

    @param nmaxx is the maximum absoulte value of the wave number in
    the x direction (default: 16)

    @param nmaxy is the maximum absoulte value of the wave number in
    the y direction (default: 16)

    @param nmaxz is the maximum absoulte value of the wave number in
    the z direction (default: 16)

    @param method string sets the evaluation method for the CUDA
    kernel.  Available methods are: all, full, 3d which compute all
    dimensions at once and planes, axes, and 1d which compute wave
    vectors by plane.  This may be extended to specific algorithms in
    the future.  This does nothing for the CPU-only computation.
*/
class Cube : public PotAccel
{

private:

  //! Coefficients are a 3-tensor
  using coefType = Eigen::Tensor<std::complex<double>, 3>;

  //! Current coefficient tensor
  std::vector<coefType> expcoef0;

  //! Coefficient tensor
  coefType expcoef, expcoefP;

  //@{
  //! Variables
  int imx, imy, imz, osize;
  int use1, use0, sampT=0, nint=0;
  std::complex<double> kfac;
  double dfac, last=-std::numeric_limits<double>::max();
  bool wrap = true;		// Enforce periodic wrapping
  //@}

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

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
  bool initialize_cuda_cube;

  //! Initialize the cuda streams
  void cuda_initialize();

  //! Zero the coefficient output vectors
  void cuda_zero_coefs();
  //@}

#endif

  //! Plane method (default: true)
  bool byPlanes;

  //! Cuda batch method (string, default: planes
  std::string cuMethod;

  //! Time routines
  class exeTimer
  {
  private:
    Cube * c;
    std::string routine;
    std::chrono::high_resolution_clock::time_point
    start0, start1, finish0, finish1;

  public:

    exeTimer(Cube* c, const std::string& name) : c(c), routine(name)
    {
      start0 = std::chrono::high_resolution_clock::now();
    }

    void Start1()
    {
      start1  = std::chrono::high_resolution_clock::now();
    }
      
    void Stop1()
    {
      finish1  = std::chrono::high_resolution_clock::now();
    }
      
    ~exeTimer()
    {
      if (c->component->timers) {
	finish0 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> duration0 = finish0 - start0;
	std::chrono::duration<double> duration1 = finish1 - start1;
      
	std::cout << std::string(60, '=') << std::endl;
	std::cout << "== " << routine << " evaluation [Cube::" << c->cC->name
		  << "] level=" << c->mlevel << std::endl;
	std::cout << std::string(60, '=') << std::endl;
	std::cout << "Time in CPU: "
		  << duration0.count()-duration1.count() << std::endl;
#if HAVE_LIBCUDA==1
	if (c->cC->cudaDevice>=0) {
	  std::cout << "Time in GPU: " << duration1.count() << std::endl;
	}
#endif
	std::cout << std::string(60, '=') << std::endl;
      }
    }
  };

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
  std::vector<coefTypePtr> expcoefN;
  std::vector<coefTypePtr> expcoefL;
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


  //! Parameter and storage initialization
  void initialize(void);
				// Parameters

  //! Minimum wave numbers in each dimesion
  //@{
  int nmaxx, nmaxy, nmaxz;
  //@}

  //! Maximum wave numbers in each dimesion
  //@{
  int nminx, nminy, nminz;
  //@}

  //@{
  //! Threading
  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);
  //@}

  //! Do the work
  void determine_acceleration_and_potential();

  //! Coefficient container instance for writing HDF5
  CoefClasses::CubeCoefs cubeCoefs;

  //! Swap coefficients
  void swap_coefs(coefType* from, coefType* to)
  {
    if (from->size() != to->size()) {
      std::ostringstream sout;
      sout << "swap_coefs: size(from)=" << from->size() << " != "
	   << "size(to)=" << to->size();
      throw std::runtime_error(sout.str());
    }

    coefType* tmp(from);
    from = to;
    to   = tmp;
  }


  //! Write basis-specific parameters to HDF5 covariance file
  void writeCovarH5Params(HighFive::File& file);

public:

  //! Id string
  std::string id;

  //! Constructor
  Cube(Component* c0, const YAML::Node& conf);

  //! Destructor
  virtual ~Cube();

  //! Compute the coefficients
  void determine_coefficients(void);

  //! Compute the force
  void get_acceleration_and_potential(Component*);

  //! Coefficient output
  void dump_coefs_h5(const std::string& file);

  /** Return a subsample of the basis coefficients and covariance
      elements for matrices for analysis. The default implementation
      returns empty vectors. This parallels the implementation in
      expui. */
  CovarElement getSubsample();

  //! Indicate whether or not subsampling is supported
  virtual bool hasSubsample() {	return true; }
};


#endif
