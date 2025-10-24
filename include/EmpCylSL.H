#ifndef _EmpCylSL_H
#define _EmpCylSL_H

#include <functional>
#include <vector>
#include <memory>
#include <limits>
#include <cmath>

#include <Eigen/Eigen>

#include "gaussQ.H"
#include "interp.H"
#include "Timer.H"

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include "Particle.H"
#include "SLGridMP2.H"
#include "coef.H"

#if HAVE_LIBCUDA==1
#include <cudaParticle.cuH>
#include <cudaMappingConstants.cuH>
#endif

#include "libvars.H"
using namespace __EXP__;

/** Encapsulatates a SLGridSph (Sturm-Liouville basis) for use as
    force method
    
    This version breaks the spherical basis computation into even and
    odd subspaces if Nodd is selected.  This gives the user more
    control over the vertical parity, is numerically more accurate,
    and more efficient.
 */
class EmpCylSL
{
public:

  using SphModTblPtr = std::shared_ptr<SphericalModelTable>;
  using SLGridSphPtr = std::shared_ptr<SLGridSph>;
  using VectorM      = std::vector<Eigen::VectorXd>;
  using MatrixM      = std::vector<Eigen::MatrixXd>;
  using ContribArray = std::vector<Eigen::VectorXd>;

protected:

  struct CylCoefHeader coefheadercyl;

  int NMAX;
  int LMAX;
  int MMAX;
  int MMIN, MLIM;
  int NMIN, NLIM;
  int NORDER;
  int NKEEP;

  unsigned nbodstot;
  std::string hallfile;

  double YMIN, YMAX;
  double dX, dY;
  int M;
  std::vector<double> cylmass1;
  bool cylmass_made;
  double cylmass;
  double minSNR, maxSNR;

  std::vector<double> r, d, m, p;

  double ASCALE;
  double HSCALE;
  double pfac, dfac, ffac;

  std::vector<Eigen::MatrixXd> facC, facS;

  int rank2, rank3;

  //@{
  //! Storage buffers for MPI
  std::vector<double> MPIin, MPIout, MPIin2, MPIout2;
  std::vector<double> MPIin_eof, MPIout_eof;

  std::vector<double> mpi_double_buf2, mpi_double_buf3;
  int MPIbufsz, MPItable;
  MPI_Status status;
  //@}

  //@{
  //! EOF variance computation
  using VarMat = std::vector< std::vector< std::vector< std::vector<double> > > >;
  VarMat SC, SS, SCe, SCo, SSe, SSo;
  //@}

  std::vector<Eigen::MatrixXd> var, varE, varO;

  std::vector< std::vector<int> > lE, lO;
  Eigen::VectorXd ev, evE, evO;
  Eigen::MatrixXd ef, efE, efO;
  Eigen::MatrixXd potd, dpot, dend;
  std::vector<Eigen::VectorXd> cosm, sinm;
  std::vector<Eigen::MatrixXd> legs, dlegs;

  SLGridSphPtr ortho;

  double Rtable, XMIN, XMAX;

  std::vector< std::vector<Eigen::MatrixXd> > potC;
  std::vector< std::vector<Eigen::MatrixXd> > densC;
  std::vector< std::vector<Eigen::MatrixXd> > rforceC;
  std::vector< std::vector<Eigen::MatrixXd> > zforceC;

  std::vector< std::vector<Eigen::MatrixXd> > potS;
  std::vector< std::vector<Eigen::MatrixXd> > densS;
  std::vector< std::vector<Eigen::MatrixXd> > rforceS;
  std::vector< std::vector<Eigen::MatrixXd> > zforceS;

  std::vector<Eigen::MatrixXd> table;

  std::vector<Eigen::MatrixXd> tpot;
  std::vector<Eigen::MatrixXd> tdens;
  std::vector<Eigen::MatrixXd> trforce;
  std::vector<Eigen::MatrixXd> tzforce;
  
  typedef std::vector<std::vector<Eigen::VectorXd>> VectorD2;
  typedef std::shared_ptr<VectorD2> VectorD2ptr;

  /** Some syntactic sugar for array of shared pointers: define an
      operator to the object of the shared pointer.  That is, for:

      MstepArray step(N);

      One can use:

      step(n)[j][k] = ...

      instead of

      (*step[n])[j][k] = ...

      to access elements.
  */
  struct MstepArray : public std::vector<VectorD2ptr>
  {
    VectorD2 & operator()(int M) { return (*(*this)[M]); }
  };

  MstepArray cosL, cosN, sinL, sinN;

  std::vector<std::vector<unsigned>> howmany1;
  std::vector<unsigned> howmany;

  std::vector<Eigen::VectorXd> accum_cos;
  std::vector<Eigen::VectorXd> accum_sin;

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


  std::vector<Eigen::MatrixXd> vc, vs;

  Eigen::MatrixXd tabp, tabf, tabd;

  //! Basis magic number
  inline static const unsigned int hmagic = 0xc0a57a1;

  //! Coefficient magic number
  inline static const unsigned int cmagic = 0xc0a57a3;

  std::vector<bool> coefs_made;
  bool eof_made;

  SphModTblPtr make_sl();

  void make_grid();
  void send_eof_grid();
  void receive_eof     (int request_id, int m);
  void compute_eof_grid(int request_id, int m);
  void compute_even_odd(int request_id, int m);
  void eigen_problem   (int request_id, int M, Timer& timer);

  void setup_eof_grid(void);
  void parityCheck(const std::string& prefix);

  //! Write HDF5 cache
  void WriteH5Cache();

  //! Read HDF5 cache
  bool ReadH5Cache();

  //! Cache versioning
  inline static const std::string Version = "1.0";

  //! The cache file name
  std::string cachefile;
				// 1=write, 0=read
				// return: 0=failure
  int    cache_grid(int, std::string file);		
  double integral(int, int, int, int);
  void   get_pot(Eigen::MatrixXd&, Eigen::MatrixXd&, double, double);
  double massR(double R);
  double densR(double R);

  Linear1d densRg, massRg;

  //! Read the cache and return YAML db
  static YAML::Node getHeader(const std::string& cachefile);
  static YAML::Node getHeader_hdf5(const std::string& cachefile);
  static YAML::Node getHeader_binary(const std::string& cachefile);

  //! Data for each harmonic subspace
  class PCAelement
  {
  public:
    //@{
    //! All the public data
    Eigen::VectorXd evalJK;
    Eigen::VectorXd meanJK;
    Eigen::VectorXd ratio;
    Eigen::MatrixXd covrJK;
    Eigen::MatrixXd evecJK;
    //@}
    
    //! Constructor
    PCAelement(int n) {
      meanJK.resize(n);
      ratio .resize(n);
      covrJK.resize(n, n);
      evecJK.resize(n, n);
    }
      
    //! Zero all data
    void reset() {
      meanJK.setZero();
      ratio .setZero();
      covrJK.setZero();
      evecJK.setZero();
    }
    
  };

  typedef std::shared_ptr<PCAelement> PCAelemPtr;

  //! PCA basis structure for caching and diagnostics
  class PCAbasis : public std::map<int, PCAelemPtr>
  {
  public:

    //! Total number in accumulation
    unsigned Tnumb;

    //! Mass in the accumulation
    double Tmass;

    //! Constructor
    PCAbasis(int M, int n)
    {
      for (int m=0; m<=M; m++) {
	(*this)[m] = PCAelemPtr(new PCAelement(n));
      }
      reset();
    }

    //! Reset all variables to zero for accumulation
    void reset()
    {
      for (auto v : *this) v.second->reset();
      Tnumb = 0;
      Tmass = 0.0;
    }

  };

  typedef std::shared_ptr<PCAbasis> PCAbasisPtr;

  //! Cache PCA information between calls
  PCAbasisPtr pb;

  pthread_mutex_t used_lock;

  //! Thread body for coef accumulation
  void accumulate_thread_call(int id, std::vector<Particle>* p, int mlevel, bool verbose);

  //! Thread body for eof accumulation
  void accumulate_eof_thread_call(int id, std::vector<Particle>* p, bool verbose);


  //! Suppress odd modes
  bool EVEN_M;

  //! Use EvenOdd partition
  bool EvenOdd;

  //! Number of even and odd terms per subspace
  int Neven, Nodd;

  //! MPI is active
  bool use_mpi;

  //! Allow older cache files with the Eigen3 row-major bug to be read
  //! with the new Eigen3 API
  bool allow_old_cache = false;

public:

  /*! Enum listing the possible selection algorithms for coefficient
    selection */
  enum TKType {
    Hall,             /*!< Tapered signal-to-noise as defined by Hall      */
    Truncate,         /*!< Truncated signal-to-noise                       */
    None              /*!< Compute the S/N but do not modify coefficients  */
  };

  //! Type of density model to use
  enum EmpModel {
    Exponential,
    Gaussian, 
    Plummer,
    Power,
    Deproject,
  };

  //! Axisymmetric disk density function for deprojection
  class AxiDisk
  {
  protected:
    std::string ID;
    std::vector<double> params;
    double M;

  public:

    //! Constructor
    AxiDisk(double M=1, std::string id="AxiDisk") : M(M), ID(id) {}

    //! Density function
    virtual double operator()(double R, double z, double phi=0.) = 0;

    //! Get ID
    std::string getID() { return ID; }

    //! Get Mass
    int getMass() { return M; }

    //! Get parameters
    std::vector<double> getParams() { return params; }
  };
  
  typedef std::shared_ptr<AxiDisk> AxiDiskPtr;

  //! TRUE if signal-to-noise methods are on (default: false)
  static bool PCAVAR;

  //! TRUE if VTK diagnostics are on (default: false)
  static bool PCAVTK;

  //! TRUE if EOF diagnostics are on (default: false)
  static bool PCAEOF;

  //! Compute taper but to not apply to coefficients (default: true)
  static bool PCADRY;

  //! VTK diagnostic frequency (default: false)
  static unsigned VTKFRQ;

  //! TRUE if we are using coordinate mapping (0=no, 1
  //! rational-function mapping (default)
  static int CMAPR;

  //! TRUE if we are using coordinate mapping (0=no, 1 for
  //! vertical hyberbolic mapping (default), 2 for power mampping)
  static int CMAPZ;

  //! TRUE if mapping is logarithmic (default: false)
  static bool logarithmic;

  //! No extrapolating beyond grid (default: false)
  static bool enforce_limits;

  //! Density model type
  static EmpModel mtype;
  
  //! Radial basis grid in radial direction
  static int NUMX;
 
  //! Radial basis grid in vertical direction
  static int NUMY;

  //! Number of bases to print (for debug)
  static int NOUT;

  //! Number of entries in radial basis table
  static int NUMR;

  //! Selector output freq (this only affects diagnostic output).
  //! Current default is to perform Hall on every step when selected
  static int HALLFREQ;

  //! Hall smoothing exponent (default: 1.0)
  static double HEXP;

  //! Minimum radial value for basis
  static double RMIN;

  //! Maximum radial value for basis
  static double RMAX;

  //! Power exponent for EOF conditioning density function
  static double PPOW;

  //! Use YAML header in coefficient file
  static bool NewCoefs;

  //! Convert EmpModel to ascii label
  static std::map<EmpModel, std::string> EmpModelLabs;

  //! Fraction of table range for basis images (for debug)
  static double HFAC;

  /** Verbose level flags
      bit   Action
      ---   ------
      0   = quiet
      1   = SL model output
      2   = EOF diagnostics
      4   = Accumulation, PCA Hall, and all other diagnostics
      8   = Debug-level communication details
      16  = Timing info
  */
  static unsigned VFLAG;		// Default=0


  //! Constructor (reset must called later)
  EmpCylSL();

  /** Constructor with parameters

      \par Parameters:

      @param multistep is the number of coefficient levels for n-body

      @param numr is the spherical radial order of the input basis

      @param lmax is the spherical angular order of the input basis

      @param mmax is the output aximuthal order for the EOF basis

      @param nord is the output radial order for the EOF basis

      @param ascale is the target disk scale LENGTH

      @param hscale is the target disk scale HEIGHT

      @param nodd is the number of vertically odd parity basis
      functions.  If unspecified, you get eigenvalue order.
      
   */
  EmpCylSL(int numr, int lmax, int mmax, int nord,
	   double ascale, double hscale, int Nodd,
	   std::string cachename);

  //! Construct from cache file
  EmpCylSL(int mlim, const std::string cache);

  //! Destructor
  ~EmpCylSL(void);

  //! Reconstruct basis with new parameters
  void reset(int numr, int lmax, int mmax, int nord,
	     double ascale, double hscale, int Nodd,
	     std::string cachename);

  //! Read EOF basis header from saved file
  int read_eof_header(const std::string& eof_file);

  //! Read EOF basis from saved file
  int read_eof_file(const std::string& eof_file);

  //! Dump the EOF basis in ascii format
  void dump_eof_file(const std::string& eof_file, const std::string& dump_file);

  //! Read basis from cache file
  int read_cache(void);

  //! Parameter access: get mmax
  int get_mmax(void) {return MMAX;}

  //! Set radial order range for evaluation
  void set_nrange(int nmin, int nmax) {NMIN=nmin; NLIM = nmax;}

  //! Set azimuthal order range for evaluation
  void set_mrange(int mmin, int mmax) {MMIN=mmin; MLIM = mmax;}

  //! Set limit to mmax for testing
  void set_mmin(int mmin) {MMIN = mmin;}

  //! Set limit to mmax for testing
  void set_mlim(int mlim) {MLIM = mlim;}

  //! Get limit to mmax for testing
  int get_mmin() { return MMIN; }

  //! Get limit to mmax for testing
  int get_mlim() { return MLIM; }

  //! Parameter access: get norder
  int get_order(void) {return NORDER;}

  //! Compute non-dimensional vertical coordinate from Z
  double z_to_y(double z);

  //! Compute Z from non-dbimensional vertical coordinate
  double y_to_z(double y);

  //! For measure transformation
  double d_y_to_z(double y);

  /** Compute deprojection of axisymmetric disk for and use this
      generate the EOF spherical basis.  The scale length must be O(1)
      with scale height H in scale length units.
   */
  void create_deprojection(double H, double Rfactor,
			   int numR, int numI, AxiDiskPtr func);

  /** Generate EOF by direct integration conditioned on a user
      supplied function

      @param numr is the number of radial knots
      @param nump is the number of azimuthal knots
      @param numt is the number of inclination knots
      @param func is the user-supplied density target function

      If func is axisymmetric, you may use nump=1 to save computation time.
  */
  void generate_eof(int numr, int nump, int numt,
		    std::function<double(double, double, double, int)> f);
  
  //! Get a single basis function values for a phase-space point
  void get_all(int m, int n, double r, double z, double phi,
	       double& p, double& d, double& fr, double& fz, double& fp);

  //! Setup for accumulated coefficients
  void setup_accumulation(int toplev=0);

  //! For PCAVAR: set subsample size
  void setSampT(int N) { defSampT = N; }

  //! For EOF
  void setup_table(void);

  //! For EOF computation
  void setup_eof(void);

  //! Clear mass counter
  void reset_mass(void);
  //@}

  //! Make coefficients from accumulated data
  //@{
  //! All levels
  void make_coefficients(bool compute=false);

  //! Single level
  void make_coefficients(unsigned mlevel, bool compute=false);

  //! Make empirical orthgonal functions
  void make_eof(void);

  //! Compute PCA
  void pca_hall(bool compute, bool subsamp);

  //! Minimum SNR coefficient value
  double getMinSNR(void) { return minSNR; }

  //! Maximum SNR coefficient value
  double getMaxSNR(void) { return maxSNR; }

  //! True if coefficients are made at all levels
  bool coefs_made_all() 
  {
    for (unsigned M=0; M<=multistep; M++) 
      if (!coefs_made[M]) return false;
    return true;
  }
  //@}

  //! Initialize PCA work space
  void init_pca();

  //! Necessary member function currently unused (change design?)
  void determine_coefficients() {};

  //! Necessary member function currently unused (change design?)
  void determine_acceleration_and_potential() {};

  //! Accumulate coefficients from particle distribution
  void accumulate(vector<Particle>& p, int mlev=0,
		  bool verbose=false, bool compute=false);

  //! Accumulate coefficients from particle distribution by thread.
  //! Used by external applications.
  void accumulate_thread(vector<Particle>& p, int mlev=0, bool verbose=false);

  //! Make EOF from particle distribution
  void accumulate_eof(vector<Particle>& p, bool verbose=false);

  //! Make EOF from particle distribution by thread.  Used by external
  //! applications.
  void accumulate_eof_thread(vector<Particle>& p, bool verbose=false);

  //! Add single particle to coefficients
  void accumulate(double r, double z, double phi, double mass,
		  unsigned long seq, int id, int mlev=0, bool compute=false);

  //! Add single particle to EOF coefficients
  void accumulate_eof(double r, double z, double phi, double mass, int id, int mlev=0);


  //! Evaluate potential and force field 
  void accumulated_eval(double r, double z, double phi, double& p0,
			double& p, double& fr, double& fz, double& fp);

  //! Evaluate density field
  double accumulated_dens_eval(double r, double z, double phi, double& d0);

  //! Evaluate peak density height within two limits
  double accumulated_midplane_eval(double r, double zmin, double zmax, double phi, int num=40);

  /** Extrapolate and sum coefficents per multistep level to get
      a complete set of coefficients for force evaluation at an
      intermediate time step
  */

  //! Set coefficients from Vectors
  void set_coefs(int mm, const Eigen::VectorXd& cos1, const Eigen::VectorXd& sin1,
		 bool zero=true);

  //! Set coefficients from std::vectors
  void set_coefs(int mm,
		 const std::vector<double>& cos1,
		 const std::vector<double>& sin1, bool zero=true);

  //! Get coefficients from Vectors
  void get_coefs(int mm, Eigen::VectorXd& cos1, Eigen::VectorXd& sin1);

  //! Get coefficients from std::vectors
  void get_coefs(int mm,
		 std::vector<double>& cos1,
		 std::vector<double>& sin1);

  //! Set cylmass manually
  void set_mass(double mass) {
    cylmass_made = true;
    cylmass = mass;
  }

  //! Dump out coefficients to output stream
  void dump_coefs(std::ostream& out);

  //! Dump out coefficients to stream in bianry format
  void dump_coefs_binary(std::ostream& out, double time);

  //! Plot basis
  void dump_basis(const std::string& name, int step, double Rmax=-1.0);

  //! Plot full fields for debugging
  void dump_images(const std::string& OUTFILE,
		   double XYOUT, double ZOUT, int OUTR, int OUTZ,
		   bool logscale);

  //! Plot basis images for debugging
  void dump_images_basis(const std::string& OUTFILE,
			 double XYOUT, double ZOUT, 
			 int OUTR, int OUTZ, bool logscale,
			 int M1, int M2, int N1, int N2);

  //! Plot PCA basis images for debugging
  void dump_images_basis_pca(const std::string& runtag,
			     double XYOUT, double ZOUT, 
			     int OUTR, int OUTZ, int M, int N, int cnt);

  //! Plot EOF basis images for debugging
  void dump_images_basis_eof(const std::string& runtag,
			     double XYOUT, double ZOUT, 
			     int OUTR, int OUTZ, int M, int N, int cnt,
			     const Eigen::VectorXd& tp);

  //! Compare current basis with another for debugging
  void compare_basis(const EmpCylSL *p);

  //! Restrict order
  void restrict_order(int n);

  /** @name Utility functions */
  // @{

  //! Compute Associated Legendre Polynomials, return Matrix type
  void legendre_R(int lmax, double x, Eigen::MatrixXd& p);
  /** Compute Associated Legendre Polynomials and derivitives, 
      return Matrix type */
  void dlegendre_R(int lmax, double x, Eigen::MatrixXd &p, Eigen::MatrixXd &dp);
  //! Compute vectors of sines and cosines by recursion
  void sinecosine_R(int mmax, double phi, Eigen::VectorXd& c, Eigen::VectorXd& s);

  // @}

  //! Convert from non-dimensional to dimensional radial coordinate
  double xi_to_r(double);

  //! Convert from dimension to non-dimension radial coordinate
  double r_to_xi(double);

  //! Jacobian
  double d_xi_to_r(double);

  //! Return current value of disk scale length
  double get_ascale(void) { return ASCALE; }

  //! Return current value of disk scale height
  double get_hscale(void) { return HSCALE; }

  //! Return current table radius
  double get_rtable(void) { return Rtable; }

  //! Set even modes only
  void setEven(bool even=true) { EVEN_M = even; }

  //! Set file name for EOF analysis and sample size for subsample
  //! computation
  inline void setHall(std::string file, unsigned tot)
  {
    hallfile = file;
    nbodstot = tot;
    
    init_pca();

    if (myid==0) {

      if (PCAVAR) {

	const std::string types[] =
	  {
	   "Hall",
	   "Truncate",
	   "None"
	  };
	
	const std::string desc[] =
	  {
	   "Compute the S/N but do not modify coefficients",
	   "Tapered signal-to-noise power defined by Hall",
	   "Compute the S/N but do not modify coefficients"
	  };
      
	std::cout << "EmpCylSL: using PCA type: " << types[tk_type]
		  << "====>" << desc[tk_type] << std::endl;
      }
      if (PCAEOF) {
	std::cout << "EmpCylSL: using PCA EOF" << std::endl;
      }
    }
  }

  //! Read and print the cache and return the header parameters as a
  //! map/dictionary
  static std::map<std::string, std::string>
  cacheInfo(const std::string& cachefile, bool verbose=true);

  /** Cache versioning override.  Allow read of old cache files
      without triggering a cache recomputation. */
  void AllowOldCache(bool allow=true) { allow_old_cache = allow; }


  //@{
  //! Get potential grid interpolated entries
  void getPotSC(int m, int j, double R, double z, double& pC, double& pS);

  //! Get density grid interpolated entries
  void getDensSC(int m, int j, double R, double z, double& dC, double& dS);

  using TableArray = std::vector<Eigen::MatrixXd>;

  //! Return density and potential matrices
  //! Indexed (M, n) pair as id = M*NORDER + n
  TableArray getDensC()
  {
    TableArray ret((MMAX+1)*NORDER);
    for (int M=0; M<=MMAX; M++) {
      for (int n=0; n<NORDER; n++) {
	int id = M*NORDER+n;
	ret[id].resize(NUMX+1, NUMY+1);
	for (int i=0; i<=NUMX; i++)
	  for (int j=0; j<=NUMY; j++)
	    ret[id](i, j) = densC[M][n](i, j);
      }
    }
    return ret;
  }

  TableArray getDensS()
  {
    TableArray ret((MMAX+1)*NORDER);
    for (int M=1; M<=MMAX; M++) {
      for (int n=0; n<NORDER; n++) {
	int id = M*NORDER+n;
	ret[id].resize(NUMX+1, NUMY+1);
	for (int i=0; i<=NUMX; i++)
	  for (int j=0; j<=NUMY; j++)
	    ret[id](i, j) = densS[M][n](i, j);
      }
    }
    return ret;
  }

  //! Return density and potential matrices
  TableArray getPotlC()
  {
    TableArray ret((MMAX+1)*NORDER);
    for (int M=0; M<=MMAX; M++) {
      for (int n=0; n<NORDER; n++) {
	int id = M*NORDER + n;
	ret[id].resize(NUMX+1, NUMY+1);
	for (int i=0; i<=NUMX; i++)
	  for (int j=0; j<=NUMY; j++)
	    ret[id](i, j) = potC[M][n](i, j);
      }
    }
    return ret;
  }

  TableArray getPotlS()
  {
    TableArray ret((MMAX+1)*NORDER);
    for (int M=1; M<=MMAX; M++) {
      for (int n=0; n<NORDER; n++) {
	int id = M*NORDER + n;
	ret[id].resize(NUMX+1, NUMY+1);
	for (int i=0; i<=NUMX; i++)
	  for (int j=0; j<=NUMY; j++)
	    ret[id](i, j) = potS[M][n](i, j);
      }
    }
    return ret;
  }


  void getPotParticle(double x, double y, double z,
		      ContribArray& vc, ContribArray& vs);

  //! Get the coefficients trimmed by a SNR value using the defined algorithm
  void get_trimmed
  (double snr,
   std::vector<Eigen::VectorXd>& ac_cos,   std::vector<Eigen::VectorXd>& ac_sin,
   std::vector<Eigen::VectorXd>* rt_cos=0, std::vector<Eigen::VectorXd>* rt_sin=0,
   std::vector<Eigen::VectorXd>* sn_rat=0);
  
  //! Set the coefficients trimmed by a SNR value using the defined algorithm
  void set_trimmed(double snr, double rem);

  //! Set number of bodies for subsample computation
  inline void setTotal(unsigned tot) {
    nbodstot = tot;
  }

  void setTK(const std::string& tk)
  {
    if      (tk == "Hall")     tk_type = Hall;
    else if (tk == "Truncate") tk_type = Truncate;
    else if (tk == "None")     tk_type = None;
    else {
      if (myid==0) {
	cout << "EmpCylSL: no such TK type <" << tk << ">"
	     << " using None type\n";
      }
    }
  }

  vector<double> sanity() { 
    vector<double> ret;
    for (int m=0; m<=MMAX; m++) ret.push_back(accum_cos[0][m]);
    return ret;
  }

  //! Check orthogonality for basis (debugging)
  void ortho_check(std::ostream& out);

  //! Check orthogonality for basis (pyEXP style)
  std::vector<Eigen::MatrixXd> orthoCheck();

  //! Return a vector of tuples of basis coefficients and the
  //! covariance matrix for subsamples of particles
  using CoefCovarType = std::tuple<Eigen::VectorXd, Eigen::MatrixXd>;
  std::vector<std::vector<CoefCovarType>> getCoefCovariance();

  //! Covariance statistics
  std::tuple<Eigen::VectorXi, Eigen::VectorXd> getCovarSamples();
  

#if HAVE_LIBCUDA==1
  cudaMappingConstants getCudaMappingConstants();

  void initialize_cuda(std::vector<cudaArray_t>& cuArray,
		       thrust::host_vector<cudaTextureObject_t>& tex);

  double get_coef(int m, int n, char c)
  {
    if (m >  MMAX)
      throw std::runtime_error("m>mmax");

    if (n >= rank3)
      throw std::runtime_error("n>=norder");

    if (c == 'c')
      return accum_cos[m][n];
    else
      return accum_sin[m][n];
  }

  double& set_coef(int mlevel, int m, int n, char c)
  {
    if (m >  MMAX)
      throw std::runtime_error("m>mmax");

    if (n >= rank3)
      throw std::runtime_error("n>=norder");

    if (c == 'c')
      return cosN(mlevel)[0][m][n];
    else
      return sinN(mlevel)[0][m][n];
  }

  double& set_coefT(int T, int m, int n)
  {
    if (m >  MMAX)
      throw std::runtime_error("m>mmax");

    if (n >= rank3)
      throw std::runtime_error("n>=norder");

    if (T >= sampT) {
      std::ostringstream sout;
      sout << "(T=" << T << ")>=(sampT=" << sampT << ")";
      throw std::runtime_error(sout.str());
    }
    
    return covV(0, T, m)[n];
  }

  double& set_covrT(int T, int m, int n, int o)
  {
    if (m >  MMAX)
      throw std::runtime_error("m>mmax");

    if (n >= rank3)
      throw std::runtime_error("n>=norder");

    if (o >= rank3)
      throw std::runtime_error("o>=norder");

    if (T >= sampT) {
      std::ostringstream sout;
      sout << "(T=" << T << ")>=(sampT=" << sampT << ")";
      throw std::runtime_error(sout.str());
    }

    return covM(0, T, m)(n, o);
  }

  double& set_tvar(int m, int i, int j)
  {
    if (m >  MMAX)
      throw std::runtime_error("m>mmax");

    if (i >= rank3 or j >= rank3)
      throw std::runtime_error("n>norder");

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

#endif

private:
  TKType tk_type;

};

extern void legendre_R(int lmax, double x, Eigen::MatrixXd& p);

extern void dlegendre_R(int lmax, double x,
			Eigen::MatrixXd &p, Eigen::MatrixXd& dp);

extern void sinecosine_R(int mmax, double phi,
			 Eigen::VectorXd& c, Eigen::VectorXd& s);

typedef std::shared_ptr<EmpCylSL> EmpCylSLptr;

#endif

// -*- C++ -*-

