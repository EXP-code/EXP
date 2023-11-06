#ifndef _Cube_H
#define _Cube_H

#include <complex>
#include <vector>
#include <string>
#include <set>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include <Coefficients.H>
#include <PotAccel.H>

#if HAVE_LIBCUDA==1
#include <thrust/complex.h>
#include <cudaParticle.cuH>
#include <cudaMappingConstants.cuH>
#endif

//! Periodic cube basis
class Cube : public PotAccel
{

private:

  //! Coefficients are a 3-tensor
  using coefType = Eigen::Tensor<std::complex<double>, 3>;

  //! Current coefficient tensor
  std::vector<coefType> expcoef, expcoefP;

  //@{
  //! Variables
  int imx, imy, imz, osize;
  int use1, use0;
  std::complex<double> kfac;
  double  dfac;
  //@}

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

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

  //! Plane method (default: true)
  bool byPlanes;

#endif

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
};


#endif
