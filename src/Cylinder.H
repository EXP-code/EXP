#ifndef _Cylinder_H
#define _Cylinder_H

#include <memory>

#include <Orient.H>
#include <Basis.H>
#include <EmpCylSL.H>
#include <CylindricalCoefs.H>

#include "config.h"

#if HAVE_LIBCUDA==1
#include <cudaParticle.cuH>
#include <cudaMappingConstants.cuH>
#endif
#include <NVTX.H>

class MixtureBasis;

/** This routine computes the potential, acceleration and density 
    using the Cylindrical biorthogonal expansion 

    Configuration parameters:

    @param rcylmin is the minimum cylindrical grid radius in scale-length units

    @param rcylmax is the maximum cylindrical grid radius in scale-length units

    @param acyl is the scale length

    @param hcyl is the scale height

    @param lmax is the maximum spherical harmonic index for EOF construction

    @param mmax is the maximum azimuthal order for the resulting basis

    @param mlim applies a limit <= mmax without basis reconstruction

    @param ncylnx is the grid size for the cylindrical radius axis

    @param ncylny is the grid size for the vertical axis

    @param ncylr radial grid size for basis computation (NUMR in EmpCylSL)

    @param ncylorder is the maximum order in for each harmonic m for the new basis

    @param ncylodd is the number of vertically antisymmetric basis functions out of the total ncylorder

    @param ncylrecomp is the frequency of basis recompution during a running simulation

    @param npca is the number steps between PCA variance/error analyses

    @param npca0 is the number steps to skip before the first PCA variance/error analysis

    @param nvtk is the number steps between vtk visualization output for the PCA variance/error analyses

    @param eof_file is the name of the basis cache file

    @param override boolean ignores the specificed cache file name and performes a recomputation

    @param samplesz is the default particle number in PCA subsampling partitions (default is 1).  The value 0 sets the sample size to sqrt(N).

    @param rnum if the size of the radial quadrature grid for analytic EOF computation

    @param pnum if the size of the azimuthal quadrature grid for analytic EOF computation (set to 1 if model is axisymmetric)

    @param tnum is the size of the equitorial quadrature grid for analytic EOF computation

    @param ashift offsets the center of the analytic model for basis function computation

    @param expcond boolean true turns on basis function reconditioning

    @param logr boolean turns on logarithmic radial basis gridding in EmpCylSL

    @param pcavar turns on variance analysis

    @param pcaeof turns on basis conditioning based on variance analysis

    @param pcadiag logs diagnostic output for variance analysis

    @param pcavtk turns on VTK visualization diagnostic output for variance analysis

    @param hexp is the Hall sharpness parameter

    @param snr is the signal-to-noise prefactor for sharpness computation

    @param evcut is the cumulative trim remainder on the covariance eigenvalues (default: -1, which implies no trimming).  Set to small positive value for trimming.
    
    @param try_cache false suppresses cache reading on restart (default: true)

    @param density builds grid with density fields (e.g. for post processing)

    @param EVEN_M true uses even harmonic orders only

    @param cmap is the coordinate mapping type (deprecated but kept for backward consistency)

    @param cmapr is the radial coordinate mapping type

    @param cmapz is the vertical coordinate mapping type

    @param self_consistent set to false turns off potential expansion; only performed the first time

    @param playback file reads a coefficient file and uses it to compute the basis function output for resimiulation

*/
class Cylinder : public Basis
{
private:

  bool expcond, EVEN_M, subsamp;
  int rnum, pnum, tnum;
  double ashift;
  unsigned int vflag;

  std::shared_ptr<EmpCylSL> ortho;
  int eof, ncompcyl;
  double cylmass, resetT;

  std::vector<Eigen::Vector3d> pos, frc;

  std::vector<double> cylmass0;
  std::vector<int> offgrid;

  void initialize(void);

				// Parameters
  double rcylmin, rcylmax, zmax, acyl;
  int nmax, lmax, mmax, mlim;
  int ncylnx, ncylny, ncylr;
  double hcyl, hexp, snr, rem;
  int ncylorder, ncylodd, ncylrecomp, npca, npca0, nvtk, cmapR, cmapZ;
  string eof_file;
  bool self_consistent, logarithmic, density, pcavar, pcainit, pcavtk, pcadiag, pcaeof, eof_over;
  bool try_cache, firstime, dump_basis, compute, firstime_coef;

  // These should be ok for all derived classes, hence declared private

  void determine_coefficients();

  void determine_coefficients_playback();

  void determine_coefficients_particles();

  void determine_coefficients_eof();

  void determine_acceleration_and_potential();

  void * determine_coefficients_thread(void * arg);

  void * determine_acceleration_and_potential_thread(void * arg);

  /** Extrapolate and sum coefficents per multistep level to get
      a complete set of coefficients for force evaluation at an
      intermediate time step
  */
  void compute_multistep_coefficients() 
  {
    if (play_back and not play_cnew) return;
    
    nvTracerPtr tPtr;
    if (cuda_prof)
      tPtr = nvTracerPtr(new nvTracer("Cylinder::compute_multistep_coefficients"));
    ortho->compute_multistep_coefficients(mfirst[mstep]);
  }

  //! Reset used particle counter and mass registers
  virtual void multistep_reset();

  //! For massive satellite simulations
  MixtureBasis *mix;

  //! Manual computation of on-grid mass
  void compute_grid_mass();

protected:

  //! For variance subsampling
  int sampT, defSampT;

  //! CUDA method for coefficient accumulation
#if HAVE_LIBCUDA==1
  virtual void determine_coefficients_cuda(bool compute_pca);
  virtual void determine_acceleration_cuda();
  virtual void multistep_update_cuda();

  thrust::host_vector  <cuFP_t> host_coefs, host_covar, host_massT;
  thrust::device_vector<cuFP_t> dev_coefs;
  std::vector<thrust::host_vector<cuFP_t>> host_coefsT, host_covarT;

  virtual void HtoD_coefs();
  virtual void DtoH_coefs(int mlevel);

  std::vector<cudaArray_t> cuInterpArray;
  thrust::host_vector<cudaTextureObject_t> tex;

  bool initialize_cuda_cyl;

  virtual void initialize_cuda()
  {
    sampT = floor(sqrt(component->CurTotal()));
    ortho->initialize_cuda(cuInterpArray, tex);
  }

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    return ortho->getCudaMappingConstants();
  }

  virtual void destroy_cuda();

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
    thrust::device_vector<cuFP_t> X_d, Y_d, m_d, u_d, p_d;
    thrust::device_vector<int   > iX_d, iY_d;

    std::vector<thrust::device_vector<cuFP_t>> T_covr;

    thrust::device_vector<cuFP_t> dN_tvar;
    thrust::device_vector<cuFP_t> dW_tvar;
    thrust::device_vector<cuFP_t> dc_tvar;
    thrust::device_vector<cuFP_t> dw_tvar;
    thrust::device_vector<cuFP_t> df_tvar;

    void resize_coefs(int ncylorder, int mmax, int N, int gridSize, int stride,
		      int sampT, bool pcavar, bool pcaeof, bool subsamp);
  };

  //! Cyclical ring interator for storage
  cudaStorage cuS;

  //! Initialize the container
  void cuda_initialize();

  //! Zero the coefficient gpu vectors
  void cuda_zero_coefs();

  //@}


#endif

  /** Coefficient playback instance.  Using shared_ptr to avoid
      constructing a simple pointer or null instance
  */
  std::shared_ptr<CylindricalCoefs> playback;

  /** Master node ships coefficients to hosts.  True (default) implies
      that only the master node caches the coefficients for playback
      to save core memory.  This is set in the config input using the
      'coefMaster: bool' parameter. Once I am sure that there are no
      algorithmic issues, I will remove this as an option.
  */
  bool coefMaster;

  //! Last playback coefficient evaluation time
  double lastPlayTime;

  //! Save data pointers
  CylindricalCoefs::DataPtr P, P1;

  //! Assign storage for coefficients
  CylindricalCoefs::DataPtr getDataPtr()
  {
    CylindricalCoefs::DataPtr ret = std::make_shared<CylindricalCoefs::Dpair>();
    std::get<0>(*ret).resize(mmax+1);
    std::get<1>(*ret).resize(mmax+1);

    for (auto & v : std::get<0>(*ret)) {
      v.resize(nmax);
      std::fill(v.begin(), v.end(), 0.0);
    }

    for (auto & v : std::get<1>(*ret)) {
      v.resize(nmax);
      std::fill(v.begin(), v.end(), 0.0);
    }
    
    return ret;
  }

  //! Store coefficients from EmpCylSL
  void getCoefs(CylindricalCoefs::DataPtr C)
  {
    for (int m=0; m<=mmax; m++) {
      ortho->get_coefs(m, std::get<0>(*C)[m], std::get<1>(*C)[m]);
    }
  }

  //! Installed stored coefficients in EmpCylSL
  void setCoefs(CylindricalCoefs::DataPtr C)
  {
    for (int m=0; m<=mmax; m++) {
      bool zero = false;
      if (m==0) zero = true;
      ortho->set_coefs(m, std::get<0>(*C)[m], std::get<1>(*C)[m], zero);
    }
  }


public:

  //! Mutexes for multithreading
  //@{
  static pthread_mutex_t used_lock, cos_coef_lock, sin_coef_lock;
  //@}

  //! Constructor
  //! \param c0 is the instantiating caller (a Component)
  //! \param conf contains the parameters passed to the method
  //! \param m is the MixtureBasis for a multicenter expansion
  //! 
  //! Parameters settable in the configuration line include:
  //! \param acyl is the disk scale length for the basis
  //! \param hcyl is the disk scale height for the basis
  //! \param nmax is the radial order for the underlying spherical basis
  //! \param lmax is the harmonic order for the underlying spherical basis
  //! \param mmax is the azimuthal harmonic order for final cylindrical expansion
  //! \param mlim restricts the azimuthal harmonic order for testing
  //! \param ncylnx is the radial grid in the basis function table
  //! \param ncylny is the vertical grid in the basis function table
  //! \param ncylr is the size of the spherical radial grid used to construct the cylindrical basis
  //! \param ncylorder is the order of the radial expansion
  //! \param ncylodd is the number of terms with vertically antisymmetric parity out of ncylorder.  If unspecified, you will get the original variance order.
  //! \param ncylrecomp is the number of steps between basis recomputation (default: -1 which means NEVER)
  //! \param npca is the number of steps between Hall coefficient recomputaton 
  //! \param npca0 is the first step for Hall coefficient computaton 
  //! \param nvtk is the number of step VTK output 
  //! \param pcadiag set to true enables PCA output diagnostics (default: false)
  //! \param eof_file is an override for the default EOF cache file name
  //! \param vflag sets verbosity (see EmpCylSL.cc)
  //! \param rnum is the number of Legendre radial knots for numerical basis computation
  //! \param pnum is the number of azimuthal knots for numerical basis computation
  //! \param tnum is the number of Legendre polar knots for numerical basis computation
  //! \param ashift is the shift applied in the x-axis relative to the original for basis conditioning
  //! \param self_consistent set to false for fixed potential
  //! \param expcond set to true for analytic basis function conditioning (default: true)
  //! \param logr tabulate basis in logarithmic coordinates (default: false)
  //! \param pcavar set to true for real-time Hall analysis (default: false)
  //! \param samplesz is the default particle number in PCA subsampling partitions (default is 1).  The value 0 sets the sample size to sqrt(N).
  //! \param pcavtk set to true dumps PCA functions in VTK format for diagnostics (default: false)
  //! \param try_cache set to true means try to read basis from cache (default: true)
  //! \param dump_basis set to true outputs basis into file
  //! \param density set to true means compute density basis in addition to potential (default: false)
  //! \param cmapR selects the radial coordinate mapping (default: 2, power mapping)
  //! \param cmapZ selects the vertical coordinate mapping (default: 1, sinh mapping)
  //!
  Cylinder(Component* c0, const YAML::Node& conf, MixtureBasis *m=0);

  //! Destructor
  virtual ~Cylinder();
  
  //! The main force call
  void get_acceleration_and_potential(Component*);

  //! Return the value for the fields in spherical polar coordinates
  void 
  determine_fields_at_point_sph(double r, double theta, double phi,
				double *tdens0, double *tpotl0, 
				double *tdens, double *tpotl, 
				double *tpotr, double *tpott, 
				double *tpotp);

  //! Return the value for the fields in cylindrical coordinates
  void 
  determine_fields_at_point_cyl(double r, double z, double phi,
				double *tdens0, double *tpotl0, 
				double *tdens, double *tpotl, 
				double *tpotr, double *tpotz, 
				double *tpotp);

  //! Return the value for the fields in cylindrical coordinates with limits
  void 
  determine_fields_at_point_cyl(double r, double z, double phi,
				double *tdens0, double *tpotl0, 
				double *tdens, double *tpotl, 
				double *tpotr, double *tpotz, 
				double *tpotp, int mmin, int mlim)
  {
    int MMIN=-1, MLIM=-1;
    // Set new values and cache current values
    if (mmin>=0) {
      MMIN = ortho->get_mmin();
      ortho->set_mmin(mmin);
    }
    if (mlim>=0) {
      MLIM = ortho->get_mlim();
      ortho->set_mlim(mlim);
    }

    determine_fields_at_point_cyl(r, z, phi,
				  tdens0, tpotl0, tdens, tpotl, 
				  tpotr, tpotz, tpotp);

    // Reset limits to original values
    if (MMIN>=0) ortho->set_mmin(MMIN);
    if (MLIM>=0) ortho->set_mlim(MLIM);
  }

  //! Return the value for the fields in Cartesian coordinates
  void 
  determine_fields_at_point(double x, double y, double z,
			    double *tdens0, double *tpotl0, 
			    double *tdens, double *tpotl, 
			    double *tpotx, double *tpoty, 
			    double *tpotz);

  /** Update the multi time step coefficient table when moving particle 
      <code>i</code> from level <code>cur</code> to level 
      <code>next</code>.
  */
  //@{
  virtual void multistep_update_begin() { 
    if (play_back and not play_cnew) return;
    if (self_consistent) ortho->multistep_update_begin();
  }
  virtual void multistep_update(int cur, int next, Component* c, int i, int id);
  virtual void multistep_update_finish() { 
    if (play_back and not play_cnew) return;
    if (self_consistent) ortho->multistep_update_finish();
  }
  //@}

  //! Print debug info
  virtual void multistep_debug();

  //! Save coefficients to file (need type marker to id dump, component id?)
  void dump_coefs(ostream& out);

  //! Sanity check on grid: dumps SM-style images of initial field
  void dump_mzero(const string& name, int step);

};


#endif
