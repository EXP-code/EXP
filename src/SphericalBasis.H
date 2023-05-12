#ifndef _SphericalBasis_H
#define _SphericalBasis_H

#include <memory>
#include <random>
#include <vector>
#include <string>
#include <set>

#include <AxisymmetricBasis.H>
#include <Coefficients.H>

#include <config_exp.h>

#if HAVE_LIBCUDA==1
#include <cudaParticle.cuH>
#include <cudaMappingConstants.cuH>
#endif
#include <NVTX.H>

#include <CoefContainer.H>

// This prints out the level switch matrix.  For deep debugging only.
// #define SPH_UPDATE_TABLE

class MixtureBasis;

//! Defines a spherical basis with members specific to spherical symmetry
class SphericalBasis : public AxisymmetricBasis
{

protected:

  //! Will cache values of the kernel required for recursion relations
  Eigen::MatrixXd krnl;

  //! Work vectors to make recursion evaluations thread safe
  std::vector<Eigen::VectorXd> u;

  //! Work vectors to make recursion evaluations thread safe
  std::vector<Eigen::VectorXd> du;

  //! Matrices per thread for obtaining potential field
  std::vector<Eigen::MatrixXd> potd;

  //! Matrices per thread for obtaining derivative of potential field
  std::vector<Eigen::MatrixXd> dpot;

  //! Matrices per thread for obtaining legendre coefficients
  std::vector<Eigen::MatrixXd> legs;

  //! Matrices per thread for obtaining derivatives legendre coefficients
  std::vector<Eigen::MatrixXd> dlegs;


  //@{
  //! Matrices per thread for obtaining expansion coefficients
  typedef std::shared_ptr<Eigen::VectorXd> VectorXdP;
  std::vector<std::vector<VectorXdP>> expcoef0;
  //@}

  //@{
  //! Mass accumulation for PCA
  vector<double> muse1;
  double muse0;
  //@}

  //! Time at last multistep reset
  double resetT;

  //! Work vectors for cosines for all values <code>m</code>
  std::vector<Eigen::VectorXd> cosm;

  //! Work vectors for sines for all values <code>m</code>
  std::vector<Eigen::VectorXd> sinm;

  //! Work matrix
  Eigen::MatrixXd dend;

  //! Cache factorial values
  Eigen::MatrixXd factorial;

  //! Initialize method
  virtual void initialize();
  //! Provided by derived classes to check consistency
  virtual void check_range();

  //! Should be called by derived classes after initialization
  void setup();

  //! Get potential and accel field values from coefficients
  void get_pot_coefs(int l, const Eigen::VectorXd& coef, double& p, double& dp);

  /** Get potential and field values from coefficients 
      (provide workspace for thread safety) */
  void get_pot_coefs_safe(int l, const Eigen::VectorXd& coef,
			  double& p, double& dp,
			  Eigen::MatrixXd& potd1, Eigen::MatrixXd& dpot1);

  //! Get density from coefficient vector
  void get_dens_coefs(int l, Eigen::VectorXd& coef, double& p);

  //! Flag whether or not work space has been initialized for coefficients
  bool firstime_coef;

  //! Flag whether or not work space has been initialized for acceleration
  bool firstime_accel;

  //! Thread method for coefficient accumulation
  virtual void * determine_coefficients_thread(void * arg);

  //! Compute the coefficients from particles
  virtual void determine_coefficients_particles(void);

  //! Compute the coefficients from an table
  virtual void determine_coefficients_playback(void);

  //! CUDA method for coefficient accumulation
#if HAVE_LIBCUDA==1
  virtual void determine_coefficients_cuda(bool compute_pca);
  virtual void determine_acceleration_cuda();
  virtual void multistep_update_cuda();

  thrust::host_vector<cuFP_t> host_coefs, host_massT;
  thrust::device_vector<cuFP_t> dev_coefs;
  cuFP_t host_mass_tot;

  virtual void HtoD_coefs(const std::vector<VectorXdP>& coef);
  virtual void DtoH_coefs(std::vector<VectorXdP>& coef);

  //! For debugging only
  void coef_cuda_compare();

  std::vector<cudaArray_t> cuInterpArray;
  thrust::host_vector<cudaTextureObject_t> tex;

  virtual void initialize_cuda() {}

  virtual void destroy_cuda();

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    cudaMappingConstants ret;

    ret.rscale = 1.0;
    ret.hscale = 0.0;
    ret.xmin   = 0.0;
    ret.xmax   = 1.0;
    ret.ymin   = 0.0;
    ret.ymax   = 0.0;
    ret.numr   = 100;
    ret.numx   = 0;
    ret.numy   = 0;
    ret.dxi    = 1.0/(ret.numr-1);
    ret.dxi    = 0.0;
    ret.cmapR  = 1;

    return ret;
  }

  virtual void initialize_mapping_constants();

  //! Compare CPU and GPU force evaluation for debugging
  void host_dev_force_compare();

  //@{
  //! Working device vectors
  thrust::device_vector<cudaTextureObject_t> t_d;

  //! Helper struct to hold device data
  struct cudaStorage
  {
    thrust::device_vector<cuFP_t> dN_coef;
    thrust::device_vector<cuFP_t> dc_coef;
    thrust::device_vector<cuFP_t> dw_coef;
    thrust::device_vector<cuFP_t> df_coef;
    thrust::device_vector<cuFP_t> dN_tvar;
    thrust::device_vector<cuFP_t> dW_tvar;
    thrust::device_vector<cuFP_t> dc_tvar;
    thrust::device_vector<cuFP_t> dw_tvar;
    thrust::device_vector<cuFP_t> df_tvar;
    thrust::device_vector<cuFP_t> plm1_d, plm2_d, r_d, m_d, u_d, a_d, p_d;
    thrust::device_vector<int   > i_d;
    
    std::vector<thrust::device_vector<cuFP_t>> T_covr;

    void resize_coefs(int nmax, int Lmax, int N, int gridSize, int stride,
		      int sampT, bool pca, bool pcaeof, bool subsamp);
    void resize_acc  (int Lmax, int Nthread);
  };

  //! A storage instance
  cudaStorage cuS;

  //! Only initialize once
  bool initialize_cuda_sph;

  //! Initialize the container
  void cuda_initialize();

  //! Zero the coefficient output vectors
  void cuda_zero_coefs();
  //@}

#endif

  /** Test change level counts for deep debugging enabled by setting
      SPH_UPDATE_TABLE at the top of this header file
  */
  std::vector<std::vector<unsigned>> occt;

  //! Print deep debugging data
  void occt_output();

  //! Thread method for accerlation compuation
  virtual void * determine_acceleration_and_potential_thread(void * arg);

  //! Compute rms coefficients
  void compute_rms_coefs(void);

  //! Add noise
  void update_noise(void);
  bool setup_noise;
				// Parameters
  //! The minimum radius
  double rmin;

  //! The maximum radius
  double rmax;

  //! Flag self_consitency
  bool self_consistent;

  //! Flag fixed monopole
  bool NO_L0;

  //! Flag drop l=1 term
  bool NO_L1;

  //! Flag use of even l harmonics only
  bool EVEN_L;

  //! Flag use of even m harmonics only
  bool EVEN_M;

  //! Use axisymmetric term s
  bool M0_only;

  //! Flag adding noise to coefficients (for testing diffusion)
  bool NOISE;  

  //! Number of "effective" points for RMS noise computation
  double noiseN;

  //! Model for RMS noise computation
  string noise_model_file;

  //! Expected coefficient values for RMS noise computation
  Eigen::VectorXd meanC;

  //! Expected squared coefficient values for RMS noise computation
  Eigen::MatrixXd rmsC;

  //@{
  //! Pseudorandom generation
  unsigned int seedN;
  std::mt19937 rgen;
  std::normal_distribution<> nrand;
  //@}

  //! Subsample fraction of particles
  double ssfrac;

  //! Flag to use subset
  bool subset;

  /** Extrapolate and sum coefficents per multistep level to get
      a complete set of coefficients for force evaluation at an
      intermediate time step
  */
  void compute_multistep_coefficients();

  //! For updating levels
  //@{
  vector< vector<Eigen::MatrixXd> > differ1;
  vector< double > pack, unpack;
  //@}

  /** Dump current coefficients (all multistep levels)
      For debugging . . .
  */
  void dump_coefs_all(ostream &);

  //! Coefficient magic number
  const unsigned int cmagic = 0xc0a57a2;

  //! For massive satellite simulations
  MixtureBasis *mix;

  /** Coefficient playback instance.  Using shared_ptr to avoid
      constructing a simple pointer or null instance
  */
  std::shared_ptr<CoefClasses::SphCoefs> playback;

  /** Master node ships coefficients to hosts.  True (default) implies
      that only the master node caches the coefficients for playback
      to save core memory.  This is set in the config input using the
      'coefMaster: bool' parameter. Once I am sure that there are no
      algorithmic issues, I will remove this as an option.
  */
  bool coefMaster;

  //! Last playback coefficient evaluation time
  double lastPlayTime;

  //! Coefficient container instance for writing HDF5
  CoefClasses::SphCoefs sphCoefs;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Use YAML header in coefficient file
  static bool NewCoefs;

  //! Constructor (<code>line</code> contains the paramter list)
  SphericalBasis(Component* c0, const YAML::Node& conf, MixtureBasis *m=0);

  //! Copy constructor
  SphericalBasis(const SphericalBasis& p);

  //! Destructor
  virtual ~SphericalBasis();
  
  //! Required member to compute accelration and potential
  /** This implemenation should be fine for most cases */
  virtual void get_acceleration_and_potential(Component* cC);

  //! Required member to compute coeifficients with threading
  /** The thread member must be supplied by the derived class */
  virtual void determine_coefficients(void);

  //! Compute coefficients with threading using the provided component
  virtual void determine_coefficients(Component *c) 
  { cC = c; determine_coefficients(); }

  //! Required member to compute accleration and potential with threading
  /** The thread member must be supplied by the derived class */
  virtual void determine_acceleration_and_potential(void);

  /** Update the multi time step coefficient table when moving particle 
      <code>i</code> from level <code>cur</code> to level 
      <code>next</code>
  */
  //@{
  virtual void multistep_reset();
  virtual void multistep_update_begin();
  virtual void multistep_update(int cur, int next, Component* c, int i, int id);
  virtual void multistep_update_finish();
  virtual void multistep_add_debug
  (const std::vector<std::vector<std::pair<unsigned, unsigned>>>& data)
  {
    for (int i=0; i<occt.size(); i++) {
      for (int j=0; j<occt[i].size(); j++)
	occt[i][j] += data[i][j].second - data[i][j].first;
    }
  }
  //@}

  /** Get derivative of potential
    \param lmax is the maximum harmonic order
    \param nmax is the maximum radial order
    \param r is the evaluation radius
    \param p will be returned array in harmonics l and radial order n
    for the potential
    \param dp will be returned array in harmonics l and radial order n
    for the derivative of the potential
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_dpotl(int lmax, int nmax, double r, Eigen::MatrixXd& p, Eigen::MatrixXd& dp, 
		 int tid) = 0;

  /** Get derivative of potential
    \param lmax is the maximum harmonic order
    \param nmax is the maximum radial order
    \param r is the evaluation radius
    \param p will be returned array in harmonics l and radial order n
    for the potential
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p,
		int tid) = 0;

  /** Get derivative of potential
    \param lmax is the maximum harmonic order
    \param nmax is the maximum radial order
    \param r is the evaluation radius
    \param d will be returned array in harmonics l and radial order n
    for the density
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& d,
		int tid) = 0;

  /** Get derivative of potential
    \param lmax is the maximum harmonic order
    \param nmax is the maximum radial order
    \param r is the evaluation radius
    \param p will be returned array in harmonics l and radial order n
    for the potential
    \param d will be returned array in harmonics l and radial order n
    for the density
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_potl_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, Eigen::MatrixXd& d,
		     int tid) = 0;

  //! Return the value for the fields in Cartesian coordinates
  void  determine_fields_at_point(double x, double y, double z,
				  double *tdens0, double *tpotl0, 
				  double *tdens, double *tpotl, 
				  double *tpotx, double *tpoty, 
				  double *tpotz);

  //! Returns the potential, its derivatives & density in spherical coords
  void determine_fields_at_point_sph(double r, double theta, double phi,
				     double *tdens0, double *tpotl0, 
				     double *tdens, double *tpotl, 
				     double *tpotr, double *tpott, double *tpotp);


  //! Returns the potential, its derivatives & density in cylindrical coords
  void determine_fields_at_point_cyl(double r, double z, double phi,
				     double *tdens0, double *tpotl0, 
				     double *tdens, double *tpotl, 
				     double *tpotr, double *tpotz, double *tpotp);


  //! Normalization for recursion relation basis
  virtual double norm(int, int) { return 1.0; }

  //! Kernel value for recursion relation basis
  virtual double knl(int, int) { return 1.0; }

  //! Dump current coefficients
  void dump_coefs(ostream &);

  //! Dump current coefficients into named HDF5 file
  void dump_coefs_h5(const std::string& file);
};

#endif
