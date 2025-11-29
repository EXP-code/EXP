#ifndef _PotAccel_H
#define _PotAccel_H

#include <cstdlib>
#include <string>
#include <chrono>
#include <list>
#include <map>

#include <memory>
#include <yaml-cpp/yaml.h>
#include <Eigen/Eigen>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>
#include <highfive/eigen.hpp>

#include "Particle.H"
#include "StringTok.H"
#include "YamlCheck.H"
#include "localmpi.H"

#include "config_exp.h"

class Component;
struct thrd_pass_PotAccel;

using namespace std;

//! Defines the basic components of the Force class
/*!
  This class provides the threading members.  Derived classes
  must provide the per thread members for evaluating coeffients
  (if any) and acceleration. 
*/
class PotAccel
{

public:
  //! For timing
  typedef std::vector<std::time_t> TList;

private:

  // Threading stuff
  thrd_pass_PotAccel *td;
  pthread_t *t;

protected:

  //! Contains parameter database
  const YAML::Node& conf;

  //@{
  /** Call to return <code>value</code> for requested parameter
      <code>name</code> */

  //! Number of bodies
  int nbodies;
  
  //! "Back" pointer to component that instantiated this force
  Component *component;

  //! Current component pointer
  Component *cC;

  //! Used by derived class to initialize any storage and parameters
  virtual void initialize(void) = 0;

  //! Number of particles used on last coefficient evaluation
  int used;

  //! Set flag to add to external potential field
  int use_external;

  //! Coefficient dump flag (true if a useful dump_coefs member is defined)
  bool coef_dump;

  //! Coefficient play back flag (true cached coefficients in use)
  bool play_back;

  //! Compute new coefficients during play back (false by default)
  bool play_cnew;

  //! Current level for multistepping
  unsigned mlevel;

  //! For thread timing
  TList timer_list;

  //! Length scaling (not used for all methods)
  double scale;

  //! Cuda implementation is available for this source method.
  //! Otherwise, particle structures will be copied between the GPU
  //! and CPU at significant expanse
  bool cuda_aware;

  //! Current YAML keys to check configuration
  std::set<std::string> current_keys;

  //! @name Subsampling
  //@{
  /** Define the type returned by getSubsample

      The outermost vectors contain subsamples of mass per sample,
      counts per sample, and arrays of coefficients and covariance
      elements per sample indexed by harmonic
  */
  using CovarElement =
    std::tuple< std::vector<double>, std::vector<unsigned>,
    std::vector<std::vector<std::tuple<Eigen::VectorXcd, Eigen::MatrixXcd>>> >;

  /** Return a subsample of the basis coefficients and covariance
      elements for matrices for PCA analysis. The default
      implementation returns empty vectors. */
  virtual CovarElement getSubsample()
  {
    std::vector<double>   mass;
    std::vector<unsigned> count;
    std::vector<std::vector<std::tuple<Eigen::VectorXcd,
				       Eigen::MatrixXcd>>> covar;
    
    return std::make_tuple(mass, count, covar);
  }

  //! Reset subsample data
  virtual void subsampleReset() { subsampleComputed = false; }

  //! Extend the subsample data
  virtual void extendSubsample(const std::string& fname, double time);

  //! Write basis-specific parameters to HDF5 covariance file
  virtual void writeCovarH5Params(HighFive::File& file)
  {
    if (myid==0)
      std::cout << "Basis::writeCovarH5Params: "
		<< "not implemented for basis type <" << id << ">"
		<< std::endl;
  }

  //! Write basis-specific parameters to HDF5 covariance file
  virtual unsigned writeCovarH5(HighFive::Group& group,
				unsigned count, double time);

  //! Variance file versioning (Version 1.0 has more granularity and
  //! poor compression)
  inline static const std::string CovarianceFileVersion = "1.1";

  //! Use 32-bit floating point type in covariance file
  bool floatType = false;

  //! Save covariance subsamples to file
  bool fullCovar = false;

  //! Is subsample computed and ready to write?
  bool subsampleComputed = false;

  //@}

  //@{
  //! HDF5 compression settings
  unsigned H5compress  = 5;
  unsigned H5chunk     = 1024*1024; // (1 MB)
  bool     H5shuffle   = true;
  bool     H5szip      = false;
  //@}


public:

  //! For timing data
  struct NData {
    int node;
    int tid;
    string name;
  };

  //! Id string
  string id;

  //! Possible geometries
  enum Geometry {sphere, cylinder, cube, slab, table, other};

  //! Map from enum to string
  static std::map<Geometry, std::string> geoname;

  //! sphere, cylinder, cube, slab, other
  Geometry geometry;

  //! Dimension of the force vector (degrees of freedom)
  int dof;

  //! Constructor
  //! \param c0 is the instantiating caller (a Component)
  //! \param line is parsed and info is made available to derived
  //! classes on request
  PotAccel(Component* c0, const YAML::Node& conf);

  //! Destructor
  virtual ~PotAccel();

  //! The main force evaluation routine to be called by clients
  virtual void get_acceleration_and_potential(Component *c) = 0;

  /** Used by get_acceleration_and_potential to request that the
      expansion be performed. */
  //@{
  virtual void determine_coefficients(void) = 0;
  virtual void determine_coefficients(Component *c) 
  { cC = c; determine_coefficients(); }
  //@}

  //! Multithreading implementation of the expansion computation
  virtual void * determine_coefficients_thread(void * arg) = 0;

  //! Multithreading implementation of the force computation
  virtual void * determine_acceleration_and_potential_thread(void * arg) = 0;

  /*! Mutex used to be used in the threading implementations to lock
    PCA variance matrix for updating */
  static pthread_mutex_t cc_lock;

  //! Flag to indicate that PCA variance matrix should be recomputed
  bool compute;

  //! For access to parent component
  virtual void RegisterComponent(Component *c) { component = c; }

  //! Thread particle counter
  std::vector<int> use;

  //! Make some threads
  void exp_thread_fork(bool coef);

  //! Make a mutex
  void make_mutex(pthread_mutex_t *m, const char *caller, const char *name);

  //! Kill a mutex
  void kill_mutex(pthread_mutex_t *m, const char *caller, const char *name);

  //! Return umber of particles used on last coefficient evaluation
  int Used() { return used; }

  //! Set external flag for current force call 
  void SetExternal() { use_external = true; }

  //! Clear external flag
  void ClearExternal() { use_external = false; }

  //! Dump coefficients for this force?
  bool HaveCoefDump() { return coef_dump; }

  //! Dump current coefficients
  virtual void dump_coefs(ostream &out) {};
  virtual void dump_coefs_h5(const std::string &file) {};

  //@name Subsampling methods
  //@{
  //! Write the subsample to file
  virtual void writeSubsample();

  //! Indicate whether or not subsamples are available
  virtual bool subsampleReady() { return subsampleComputed; }

  //! Indicate whether or not subsampling is supported
  virtual bool hasSubsample() {	return false; }
  //@}

  //! HDF5 compression settings
  void setCovarH5Compress(unsigned level, unsigned chunksize, bool shuffle, bool szip=false)
  {
    H5compress = level;
    H5chunk    = chunksize;
    H5shuffle  = shuffle;
    H5szip     = szip;
  }

  /** Update the multi time step force algorithm when moving particle 
      <code>i</code> from level <code>cur</code> to level 
      <code>next</code>
  */
  
  //! Execute to begin level shifts for particles
  virtual void multistep_update_begin() {}

  //! Implementation of level shifts
  virtual void multistep_update(int cur, int next, Component* c, int i, int id) {}

  //! Execute to finish level shifts for particles
  virtual void multistep_update_finish() {}

#if HAVE_LIBCUDA==1
  //! Cuda implementation of level shifts
  virtual void multistep_update_cuda() {}
#endif

  //! Set the level for evaluating acceleration at level n
  virtual void set_multistep_level(unsigned n) { mlevel=n; }

  //! Reset data for multistep
  virtual void multistep_reset() {}

  //! Print debug info
  virtual void multistep_debug() {}

  //! Set new length scale (default: 1.0)
  virtual void setScale(double s) { scale = s; }

  //! Get current length scale
  virtual double getScale() { return scale; }

  //! Check whether coeffcients are not available
  bool NoCoefs() { return play_back and !play_cnew; }

  //! For timing the fraction spent in threaded force methods
  //@{
  void print_timings(const string&, TList& tlist);
  void print_timings(const string& label);
  void thread_timing_beg(int id);
  void thread_timing_end(int id);
  //@}

  //! For parameter parsing
  bool atol(string& val) { return atoi(val.c_str()) ? true : false; }

  //! Cuda aware
  bool cudaAware() { return cuda_aware; }

  //! Get unaccounted keys
  std::set<std::string> unmatched() { return current_keys; }

  //! Round time key to emulated fixed-point arithmetic
  double roundTime(double time);

};

//! Helper class to pass info to threaded member
struct thrd_pass_PotAccel 
{
  //! Pointer to force
  PotAccel *t;
  //! Call coefficient thread member (true) or acceleration/potential (false)
  bool coef;
  //! Thread counter id
  int id;
};


//! Used for thread timing only (VERBOSE>5)
bool ltdub(const pair<double, PotAccel::NData>& A, 
	   const pair<double, PotAccel::NData>& B);

#endif // _PotAccel_H
