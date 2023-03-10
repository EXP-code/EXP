#ifndef _PolarBasis_H
#define _PolarBasis_H

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
#include <BiorthCyl.H>

class MixtureBasis;

/** Defines a cylindrical basis with members specific to polar symmetry

    Configurable parameters:

    @param scale is the characterisic scale of the basis
    @param rmin is the minimum radius of the basis
    @param rmax is the maximum radius of the basis
    @param self_consistent true for updating the the coefficients at time steps
    @param NO_M0 true to omit the m=0 harmonic (default: false)
    @param NO_M1 true to omint the m=1 harmonic (default: false)
    @param EVEN_M true to include only even m harmonics (default: false)
    @param M0_ONLY true to include only the m=0 harmonic (default: false)
    @param ssfrac set > 0.0 to compute a fraction of particles only
    @param playback true to replay from a coefficient file
    @param coefMaster true to have only the root node read and distribute the coefficients

    Other parameters may be defined and passed to any derived classes
    in addition to these.
*/
class PolarBasis : public AxisymmetricBasis
{

public:
  
  using VectorM      = std::vector<Eigen::VectorXd>;
  using MatrixM      = std::vector<Eigen::MatrixXd>;
  using ContribArray = std::vector<Eigen::VectorXd>;

protected:

  //! Work vectors to make recursion evaluations thread safe
  std::vector<Eigen::VectorXd> u;

  //! Work vectors to make recursion evaluations thread safe
  std::vector<Eigen::VectorXd> du;

  //! Matrices per thread for obtaining potential field
  std::vector<Eigen::MatrixXd> potd;

  //! Matrices per thread for obtaining derivative of potential field
  std::vector<Eigen::MatrixXd> dpotR, dpotZ;

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

  //! Initialize method
  virtual void initialize();
  //! Provided by derived classes to check consistency
  virtual void check_range();

  //! Should be called by derived classes after initialization
  void setup();

  //! Initialize for PCA variance computation
  void init_pca();

  //! Clear arrays for accumulation
  void setup_accumulation(int mlevel);

  //! Get potential and accel field values from coefficients
  void get_pot_coefs(int m, const Eigen::VectorXd& coef, double& p, double& dpr, double& dprz);

  /** Get potential and field values from coefficients 
      (provide workspace for thread safety) */
  void get_pot_coefs_safe(int m, const Eigen::VectorXd& coef,
			  double& p, double& dpr, double& dprz,
			  Eigen::MatrixXd& potd1, Eigen::MatrixXd& dpotr1, Eigen::MatrixXd& dpotz1);

  //! Get density from coefficient vector
  void get_dens_coefs(int m, Eigen::VectorXd& coef, double& p);

  //! Model is flat
  bool is_flat;

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

  std::vector<std::vector<unsigned>> howmany1;
  std::vector<unsigned> howmany;

  //! Add single particle to coefficients
  void accumulate(double r, double z, double phi, double mass,
		  unsigned long seq, int id, int mlev=0, bool compute=false);

  //! Accumulate coefficients from particle distribution
  void accumulate(std::vector<Particle>& part, int mlevel, bool verbose, bool compute);

  //! Accumulate coefficients from particle distribution by thread.
  //! Used by external applications.
  void accumulate_thread(std::vector<Particle>& part, int mlevel, bool verbose);

  //! Threaded call for the above call
  void accumulate_thread_call(int id, std::vector<Particle>* p, int mlevel, bool verbose);

  std::vector<bool> coefs_made;

  //@{
  /** More syntactic sugar for arrays of arrays: define an operator to
      the object of the main nested arrays.  That is, for:

      MatrixArray arr(N);
      for (auto & v : arr) v.resize(M);

      One can use:

      arr(n, m)[j][k] = ...

      instead of

      (*arr[n][m])[j][k] = ...

      to access elements.
  */
  struct CoefVector : public std::vector<std::vector<VectorM>>
  {
    Eigen::VectorXd & operator()(int i, unsigned j, int m) { return ((*this)[i][j][m]); }
  };

  struct CoefMatrix : public std::vector<std::vector<MatrixM>>
  {
    Eigen::MatrixXd & operator()(int i, unsigned j, int m) { return ((*this)[i][j][m]); }
  };
  //@}


  // Test for eof trim
  std::vector<std::vector<Eigen::MatrixXd>> tvar;

  // Error analysis
  CoefVector  covV;
  CoefMatrix  covM;

  std::vector< std::vector<unsigned>  > numbT1;
  std::vector< std::vector<double>  > massT1;
  std::vector<unsigned> numbT;
  std::vector<double> massT;
  unsigned sampT, defSampT;

  std::vector<double> cylmass1;
  bool cylmass_made;
  double cylmass;
  double minSNR, maxSNR;

  std::vector<Eigen::MatrixXd> vc, vs;

  virtual void get_pot(Eigen::MatrixXd& Vc, Eigen::MatrixXd& Vs, double r, double z) = 0;

  //! For updating levels
  //@{
  std::vector< std::vector<Eigen::MatrixXd> > differS1, differC1;
  std::vector<double> workC1, workC, workS1, workS;
  //@}


  //! Get table scale
  virtual double getRtable() = 0;

  //! CUDA method for coefficient accumulation
#if HAVE_LIBCUDA==1
  virtual void determine_coefficients_cuda(bool compute_pca);
  virtual void determine_acceleration_cuda();
  virtual void multistep_update_cuda();

  thrust::host_vector<cuFP_t> host_coefs, host_massT;
  thrust::device_vector<cuFP_t> dev_coefs;
  std::vector<thrust::host_vector<cuFP_t>> host_coefsT, host_covarT;
  cuFP_t host_mass_tot;

  virtual void HtoD_coefs();
  virtual void DtoH_coefs(unsigned);

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
    thrust::device_vector<cuFP_t> m_d, u_d, X_d, Y_d, a_d, p_d;
    thrust::device_vector<int   > iX_d, iY_d;
    
    std::vector<thrust::device_vector<cuFP_t>> T_covr;

    void resize_coefs(int nmax, int Lmax, int N, int gridSize, int stride,
		      int sampT, bool pca, bool pcaeof, bool subsamp);
    void resize_acc  (int Lmax, int Nthread);
  };

  double get_coef(int m, int n, char c)
  {
    if (m >  Mmax)
      throw std::runtime_error("m>mmax");
    
    if (n >= nmax)
      throw std::runtime_error("n>=nmax");
    
    if (c == 'c') {
      if (m==0) return (*expcoef[0])[n];
      else      return (*expcoef[2*m-1])[n];
    } else {
      return (*expcoef[2*m])[n];
    }
  }
  
  double& set_coef(int mlevel, int m, int n, char c)
  {
    if (m >  Mmax)
      throw std::runtime_error("m>mmax");
    
    if (n >= nmax)
      throw std::runtime_error("n>=nmax");
    
    if (c == 'c') {
      if (m==0) return (*expcoef0[0][0])[n];
      else      return (*expcoef0[0][2*m-1])[n];
    } else {
      return (*expcoef0[0][2*m])[n];
    }
  }
    
  double& set_coefT(int T, int m, int n)
  {
    if (m >  Mmax)
      throw std::runtime_error("m>mmax");
    
    if (n >= nmax)
      throw std::runtime_error("n>=nmax");
    
    if (T >= sampT) {
      std::ostringstream sout;
      sout << "(T=" << T << ")>=(sampT=" << sampT << ")";
      throw std::runtime_error(sout.str());
    }
    
    return (*expcoefT[T][m])[n];
  }
  
  double& set_covrT(int T, int m, int n, int o)
  {
    if (m >  Mmax)
      throw std::runtime_error("m>mmax");
    
    if (n >= nmax)
      throw std::runtime_error("n>=nmax");
    
    if (o >= nmax)
      throw std::runtime_error("o>=nmax");
    
    if (T >= sampT) {
      std::ostringstream sout;
      sout << "(T=" << T << ")>=(sampT=" << sampT << ")";
      throw std::runtime_error(sout.str());
    }
    
    return (*expcoefM[T][m])(n, o);
  }
    
  double& set_tvar(int m, int i, int j)
  {
    if (m >  Mmax)
      throw std::runtime_error("m>mmax");

    if (i >= nmax or j >= nmax)
      throw std::runtime_error("n>nmax");

    return tvar[0][m](i, j);
  }

  double& set_massT(int T)
  {
    if (T >= sampT) {
      std::ostringstream sout;
      sout << "(T=" << T << ")>=(sampT=" << sampT << ")";
      throw std::runtime_error(sout.str());
    }

    return massT1[0][T];
  }

  //! A storage instance
  cudaStorage cuS;

  //! Only initialize once
  bool initialize_cuda_plr;

  //! Initialize the container
  void cuda_initialize();

  //! Zero the coefficient output vectors
  void cuda_zero_coefs();
  //@}

#endif

  //! Thread method for accerlation compuation
  virtual void * determine_acceleration_and_potential_thread(void * arg);

  //! Maximum harmonic order restriction
  int mlim;

  //! Basis linear scale for mapping
  double scale;

  //! The minimum radius
  double rmin;

  //! The maximum radius
  double rmax;

  //! Flag self_consitency
  bool self_consistent;

  //! Flag fixed monopole
  bool NO_M0;

  //! Flag drop m=1 term
  bool NO_M1;

  //! Flag use of even m harmonics only
  bool EVEN_M;

  //! Use axisymmetric terms
  bool M0_only;

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
  std::shared_ptr<CoefClasses::CylCoefs> playback;

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
  CoefClasses::CylCoefs cylCoefs;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

  //! Mutexes for multithreading
  //@{
  static pthread_mutex_t used_lock, cos_coef_lock, sin_coef_lock;
  //@}

public:

  //! Use YAML header in coefficient file
  static bool NewCoefs;

  //! Constructor (<code>line</code> contains the paramter list)
  PolarBasis(Component* c0, const YAML::Node& conf, MixtureBasis *m=0);

  //! Copy constructor
  PolarBasis(const PolarBasis& p);

  //! Destructor
  virtual ~PolarBasis();
  
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
  //@}

  /** Get derivative of potential
    \param r is the evaluation radius
    \param z is the evaluation height
    \param p will be returned array in harmonics l and radial order n
    for the potential
    \param dpr will be returned array in harmonics l and radial order n
    for the derivative of the potential
    \param dpz will be returned array in harmonics l and radial order n
    for the derivative of the potential
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_dpotl(double r, double z, Eigen::MatrixXd& p,
		 Eigen::MatrixXd& dpr, Eigen::MatrixXd& dpz, int tid) = 0;

  /** Get derivative of potential
    \param r is the evaluation radius
    \param z is the evaluation height
    \param p will be returned array in harmonics l and radial order n
    for the potential
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_potl(double r, double z, Eigen::MatrixXd& p, int tid) = 0;

  /** Get derivative of potential
    \param r is the evaluation radius
    \param z is the evaluation height
    \param d will be returned array in harmonics l and radial order n
    for the density
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_dens(double r, double z, Eigen::MatrixXd& d, int tid) = 0;

  /** Get the potential and density
    \param r is the evaluation radius
    \param p will be returned array in harmonics l and radial order n
    for the potential
    \param d will be returned array in harmonics l and radial order n
    for the density
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_potl_dens(double r, double z,
		     Eigen::MatrixXd& p, Eigen::MatrixXd& d, int tid) = 0;

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
