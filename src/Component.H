#ifndef _Component_H
#define _Component_H

#include <mpi.h>
#include <utility>

#include <yaml-cpp/yaml.h>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

#include <EXPException.H>
#include <global.H>
#include <Orient.H>
#include <header.H>
#include <localmpi.H>
#include <ParticleFerry.H>
#include <CenterFile.H>
#include <PotAccel.H>
#include <Circular.H>
#include <Timer.H>

#include <config_exp.h>

#if HAVE_LIBCUDA==1
#include <cudaParticle.cuH>  
#endif


//! Structure used to sort old & new indices for load balancing
struct loadb_datum
{
  int top;		///< Top particle index for this process
  int indx;		///< Process ID
  unsigned short s;	///< 0-->Old partiion, 1-->New partiion
};

//! Class associating the phase space with a biorthogonal basis
/*!
  The Component provides "interactors" for specific particle fields
  <br>
  Each field in the Component config file format contains four fields,
  each delimited by a colon (":"), and new lines are not significant:
  <br>
  <ol>
  <li> name of the potential solver (force methods) for this phase space
  <li> component-specific parameters [see Parameter section below]
  <li> file name containing phase space to be read
  <li> comma-delimited parameter list passed to the potential solver
  </ol>
  <br>
  Currently available force methods:
  \par bessel 
  %Bessel function basis, see Bessel 

  \par sphereSL 
  Weinberg's adaptive spherical basis derived by numerical solution to the Sturm-Liouville equation, see Sphere

  \par cube 
  Trigonmetric basis with periodic boundaries, see Cube

  \par slabSL 
  Trigonmetric basis with periodic boundaries the horizontal
  directions and vacuum boundary conditions in the vertical dirction,
  see Slab

  \par slab 
  Trigonmetric basis with periodic boundaries the horizontal
  directions and vacuum boundary conditions in the vertical dirction
  based on a direct solution of the Sturm-Liouville equation
  perpendicular to the slab, see SlabSL

  \par cylinder 
  Weinberg's adaptive cylindrical basis derived by numerical solution
  to the Sturm-Liouville equation, see Cylinder

  \par direct
  The usual n<sup>2</sup> force calculation, see Direct 

  \par noforce
  No force (presumably, you are supplying your own force)

  @param EJ enables energy weighed center-of-mass.  The paramter value
  is a two-bit flag: AXIS=1 and CENTER=2.
    
  @param nEJkeep is the number of time steps in running box car average for
  EJ-determined parameters

  @param nEJwant is the target number of particles to keep in the
  ensemble used to determine the paramters.  The energy cutoff is
  dynamically adjusted to achieve this target.

  @param nEJwant is the size of the past states used to estimate the
  the acceleration of the expansion frame

  @param EJx0 is the initial EJ center x-coordinate (default: 0)

  @param EJy0 is the initial EJ center y-coordinate (default: 0)

  @param EJz0 is the initial EJ center z-coordinate (default: 0)

  @param EJu0 is the initial EJ center u-velocity   (default: 0)

  @param EJv0 is the initial EJ center v-velocity   (default: 0)

  @param EJw0 is the initial EJ center w-velocity   (default: 0)

  @param EJdT is the spacing in time between samples   (default: 0)

  @param EJkinE set true includes kinetic energy in the energy
  determination, otherwise potential energy alone is used

  @param EJext set true includes the external potential energy in the
  energy computation

  @param EJdiag turns on verbose EJ reporting

  @param EJdryrun set to true does not apply EJ centering or rotation but
  computes the transformation

  @param EJlinear set true moves center with constant velocity.  For
  testing only!! (default: false)

  @param EJdamp multiplies the least squares slope to prevent instability
  (default: 1.0, no damping)

  @param rmax is the maximum radius particle.  This is determined
  dynamically from the initial conditions but can be overridden as a
  parameter.

  @param rtrunc is the radius from the <code>center</code> beyond
  which it does not contribute to the force expansion (but does
  continues to be affected by it).  The default value is
  <code>1.0e20</code>.

  @param rcom is the radius beyond which a particle is not used
  to compute the center of mass coordinate system

  @param ton is the midpoint of the adiabatic turn on for this
  component

  @param twid is the width for the adiabatic turn on.  Unless
  <code>ton</code> or <code>toff</code> is specified, the component
  begins and remains fully on.

  <br>
  Reference frames:
  <ol>

  <li> <em>Inertial</em> is the rest reference frame whose coordinate
  system is centered at (0, 0, 0)

  <li> <em>Local</em> is the coordinate system of specific to the
  particles of an individual Component instance.  If com_system is
  defined, this is defined to have its origin at the center of mass
  and the center of velocity of the component's partcles at all times.
  If <code>com_system</code> is off, this is the same as
  <em>Inertial</em>

  <li> <em>Centered</em> offsets the particles coordinates using the
  current centering scheme for the expansion.  At this point, this is
  the EJ method selectable by the user.  The SphereTwoCenter method
  indirectly manipulates this to achieve two separate expansion
  centers.

  <li> <em>keyPos</em> is the value of the species key in the particle
  integer attribute vector

  <li> <em>pbufsiz</em> is the number of particles in each particle
  buffer sent to disk when using MPI-IO with OutPSP

  <li> <em>blocking</em> set to true chooses blocking writes when
  using MPI-IO with OutPSP

  <li> <em>buffered</em> set to true chooses buffered writes when
  using POSIX binary

  <li> <em>ignore</em> the PSP info stanza on restart (i.e. for
  specifying alternative parameters, using an old-style PSP file, or
  starting a new job with a previous output from another simulation)
  </ol>
*/
class Component
{

  friend class ScatterMFP;
  friend class OrbTrace;
  friend class ComponentContainer;

private:

  //! Parallel distribute and particle io
  void load_balance(void);
  void update_indices(void);
  void read_bodies_and_distribute_ascii(void);
  void read_bodies_and_distribute_binary_out(istream *);
  void read_bodies_and_distribute_binary_spl(istream *);

  //! Helper for split-PSP-file particle consumer
  void openNextBlob(std::ifstream& in,
		    std::list<std::string>::iterator& fit, int& N);


  //! For magic number checking
  const static unsigned long magic = 0xadbfabc0;
  const static unsigned long mmask = 0xf;
  const static unsigned long nmask = ~mmask;
  unsigned rsize;


  //@{
  //! Internal data for sequence checking
  bool indexing, aindex, umagic;
  int seq_beg, seq_end, seq_cur;
  //@}

  //! Particle buffer count for MPI-IO writing
  int pBufSiz;

  //! Using blocking writes for MPI-IO
  bool blocking;

				// Used in computing target rates for
				// load balancing
  int round_up(double dnumb);

  // Particle buffer
  std::vector<PartPtr> pbuf;

  // Adiabatic mass
  bool adiabatic;
  double ton, toff, twid;

  // For load balancing
  vector <loadb_datum> loadb;
  void add_particles(int from, int to, std::vector<PartPtr>& plist);

  // Compute initial com position and velocity from phase space
  void initialize_com_system();
  vector<double> com_lev, cov_lev, coa_lev, com_mas, angmom_lev;
  vector<double> comE_lev, covE_lev, comE_mas;

  // Momentum tracking
  bool consp;
  int tidal;

  // For exchanging particles
  ParticleFerryPtr pf;

protected:

  //! Set configuration and force
  void configure(void);

  /** Initialize parameters (called by ComponentContainer after component
      instantiation */
  void initialize(void);

  //! Enumerate and initialize cuda devices
  void initialize_cuda(void);

  //! Slab parameters
  double azcm_slab;

  //! Bodies on this node
  unsigned int nbodies;

  //! Bodies on all nodes
  unsigned int nbodies_tot;

  /** Used by gather and distribution routines to define a particle
      structure to MPI */
  MPI_Datatype Particletype;

  //! Setup distribution vectors
  void setup_distribution(void);

  //! Vector holding current particle partition
  //@{
  vector<unsigned int> nbodies_table, nbodies_index;
  //@}

  //! Vectors for holding the rates per node for load balancing
  //@{
  vector<double> orates, trates;
  //@}

  //! Binary read flag
  bool binary;

  //! COM log file
  string comfile;

  //! COM flag
  int com_restart;

  //! COM restart function
  void restart_com_system(void);

  //! Points to associated particles using sequence number
  PartMap particles;

  //! Level occupation output
  int nlevel;

  //! Use center from this component
  Component *c0;

  //! Locagte centering component
  void find_ctr_component();

  //! Name of the centering component
  std::string ctr_name;

  //! Center file playback
  std::shared_ptr<CenterFile> CF;

  //! Write binary component phase-space structure using MPI (blocking)
  void write_binary_mpi_b(MPI_File& out, MPI_Offset& offset, bool real4 = false);

  //! Write binary component phase-space structure using MPI (non-blocking)
  void write_binary_mpi_i(MPI_File& out, MPI_Offset& offset, bool real4 = false);

  //! Set default, unset parameters in config
  void set_default_values();
  
  //! Particles waiting for a new index
  std::list<PartPtr> new_particles;

  //! Maximum sequence value so far
  unsigned long top_seq;

  //! Update sequence number of added particles
  void seq_new_particles();

  //! Particle number update is needed
  unsigned modified=0;

  //! No multistep switching
  bool noswitch;

  //! Compute new levels on the first step only
  bool freezeLev;

  //! Select level from criteria over last time step
  bool dtreset;

  //@{
  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys_top;
  static const std::set<std::string> valid_keys_parm;
  static const std::set<std::string> valid_keys_force;
  //@}

  
  //@{
  //! HDF5 compression parameters
  int H5compress = 0;
  int H5chunk = 4096;
  //@}

public:

  //! Describe the phase space coordinates
  enum PSTYPE {Inertial=1, Local=2, Centered=4};

  //! Configuration for this component
  YAML::Node conf;

  //! Identify this component
  std::string name;

  //! Force method
  std::string id;

  //! Parameters for component
  YAML::Node cconf;

  //! Parameters for force
  YAML::Node fconf;

  //! Particle file
  std::string pfile;

  //! Slab geometry
  bool slab;

  //! Periodic cube
  bool cube;

  //! Point masses
  bool point_mass;

  //! Cylindrical geometry
  bool disk;
  
  //! Spherical geometry
  bool sphere;
  
  //! Maximum initial particle radius
  double rmax;

  //! Species key position in particle integer attribute vector
  int keyPos;

  /**
     These are pointers to the center of mass, center of velocity,
     angular momentum and the phase-space 6 vector.
  */
  //@{
  //! Center of mass
  double *com;
  //! Center of velocity
  double *cov;
  //! Center of acceleration
  double *coa;
  //! Center (e.g. set by Orient)
  double *center;
  //! Angular momentum vector
  double *angmom;
  //! Phase space vector
  double *ps;
  //@}

  //@{
  //! Pseudo-accleration (set by Orient)
  Eigen::Vector3d accel, omega, domdt, pseudo;
  Eigen::Vector3d& getPseudoAccel(double* pos, double* vel);
  //@}

  /**
     @name Center of mass system coordinates
  */
  //@{
  //! System mass, center of mass, center of velocity and acceleration
  double mtot0, *com0, *cov0, *acc0;
  //! Maximum radius for contribution to center-of-mass computation
  double rcom;

  //@}

  //! Total mass
  double mtot;

  //! Dimension of the phase space
  int dim;

  /** Particle list per level
      This will be a sorted list for each level.  We use a vector rather
      than (e.g.) a set for the random insertion
  */
  std::vector< vector<int> > levlist;

  //! Multstep dt type counter
  std::vector< vector<unsigned> > mdt_ctr;

  //! Points to associated force instance
  PotAccel *force;

  //! Attribute array dimensions for integer parameters
  int niattrib;

  //! Attribute array dimensions for float-valued parameters
  int ndattrib;

  //! Binding energy orientation flag
  int EJ;

  //! Number to keep in orientation stack
  int nEJkeep;

  //! Target number of particles for average
  int nEJwant;

  //! Target number of states for pseudo-acceleration estimation
  int nEJaccel;

  //! Initial EJ center
  //@{
  //! x-coord
  double EJx0;
  //! y-coord
  double EJy0;
  //! z-coord
  double EJz0;
  //@}

  //! Initial EJ velocity
  //@{
  //! x-coord
  double EJu0;
  //! y-coord
  double EJv0;
  //! z-coord
  double EJw0;
  //@}

  //! Spacing between EJ samples
  double EJdT;

  //! Print verbose diagnostics
  bool EJdiag;

  //! Do not apply transformation
  bool EJdryrun;

  //! Use particle KE for energy computation
  bool EJkinE;

  //! Use external potential for energy computation
  bool EJext;

  //! Move center at constant velocity specified by (EJu0, EJv0, EJw0)
  bool EJlinear;

  //! Damping factor for center determination (< 1.0 for damping)
  double EJdamp;

  //! Use center of mass coordinate system
  bool com_system;

  //! Log the center of mass coordinate system at each step
  bool com_log;

  //! Maximum radius from center for force evaluation
  double rtrunc;

  //! CUDA timing info
  bool timers;

  //! Use buffered binary writes
  bool buffered;

  //! Orientation cache
  Orient *orient;

  //! The constructor
  /*! 
    Takes the id string id and phase space file name as arguments
  */
  Component(YAML::Node& conf);

  //! Initialize from file
  Component(YAML::Node&, istream *in, bool SPL);

  //! Destructor
  ~Component();

  //! Parallel gather
  PartPtr * get_particles(int* number);
  
  //! Retrieve particle count on a particular node (one node at a time)
  int particle_count(int node) { 
    if (node<0 || node>=numprocs) return 0;
    return nbodies_table[node]; 
  }

  //! Retrieve particle count on a particular node (all nodes at once)
  vector<unsigned int> particle_count() { return nbodies_table; }

  //! Write binary component phase-space structure
  void write_binary(ostream *out, bool real4 = false);
  
  //! Write HDF5 phase-space structure for real or float
  template <typename T>
  void write_HDF5(HighFive::Group& group, bool masses, bool IDs);

  //! Write header for per-node writes
  void write_binary_header(ostream* out, bool real4, std::string prefix, int nth=1);

  //! Write particles for per-node writes
  void write_binary_particles(std::ostream* out, bool real4);

  //! Write particles for per-node writes (multithreaded version)
  void write_binary_particles(std::vector<std::shared_ptr<std::ofstream>>& out, bool real4);

  //! Write particles for per-node writes with multithreading
  void write_binary_particles(std::ostream* out, int threads, bool real4);

  //! Write binary component phase-space structure using MPI
  void write_binary_mpi(MPI_File& out, MPI_Offset& offset, bool real4 = false)
  {
    if (blocking)
      write_binary_mpi_b(out, offset, real4);
    else
      write_binary_mpi_i(out, offset, real4);
  }
  
  //! Write ascii component phase-space structure
  void write_ascii(ostream *out, bool accel = false);
  
  //! Redestribute this component
  void redistribute_particles(void);

  //! Compute center of mass and center of velocity (CPU version)
  void fix_positions_cpu(unsigned mlevel=0);

#if HAVE_LIBCUDA==1
  //! Compute center of mass and center of velocity (GPU version)
  void fix_positions_cuda(unsigned mlevel=0);
#endif

  //! Compute center of mass and center of velocity
  void fix_positions(unsigned mlevel=0)
  {
#if HAVE_LIBCUDA==1
    if (use_cuda) fix_positions_cuda(mlevel);
    else
#endif
      fix_positions_cpu(mlevel);
  }

  /** Redistribute particle by list (deprecated?)

      Serialize redistribution list into an integer array
      of numprocs stanzas, each with the format
      n -- current node
      M -- number to redistribute (may be zero)
      index_1
      tonode_1
      index_2
      tonode_2
      \.
      \.
      index_M
      tonode_M

      so each stanza has 2(M+1) integers
  */
  void redistributeByList(vector<int>& redist);

  //! Subtract mean acceleration from each particle accel
  void update_accel(void);

  //! Update angular momentum values
  void get_angmom(unsigned mlevel=0);

  //! Adiabatic turn on factor, range in [0, 1]
  double Adiabatic(void);

  //! Is particle out of bounds?
  bool freeze(unsigned indx);

  //! Is particle beyond COM system
  bool escape_com(const Particle&);

  //! Access to particle vector
  PartMap& Particles() {
    return particles;
  }

  //! Access to particle as a pointer
  Particle *Part(unsigned long i) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    return tp->second.get();
  }

  //! Access to particle via the shared pointer
  PartPtr partPtr(unsigned long i) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    return tp->second;
  }

  //! Check if particle index is in the component
  bool partExists(unsigned long i) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      return false;
    }
    return true;
  }

  //! Generate new particle: this adds a new particle and returns a
  //! pointer to the same particle.  It is the caller's responsbility
  //! to populate the phase space and particle attributes with sane
  //! values.
  Particle *GetNewPart();

  //! Add a particle to the component
  void AddPart(PartPtr p);

  //! Remove a particle from the component
  void DestroyPart(PartPtr p);

  //! Erase a particle with no level list check
  void ErasePart(PartPtr p)
  {
    particles.erase(p->indx);
    nbodies = particles.size();
  }
  
  //! Particle vector size
  unsigned Number() {
    return particles.size();
  }

  //! Return the last particle total (no MPI)
  unsigned CurTotal()
  {
    return nbodies_tot;
  }

  //! Update total number of particles in component.  Should be called
  //! after adding or killing particles.  Uses MPI.
  unsigned NewTotal() {
    MPI_Allreduce(MPI_IN_PLACE, &modified, 1, MPI_UNSIGNED, MPI_SUM,
		  MPI_COMM_WORLD);
    if (modified) seq_new_particles();
    return nbodies_tot;
  }

  //! Access to mass
  inline double Mass(int i) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    return tp->second->mass;
  }

  //! Access to positions
  inline double Pos(int i, int j, unsigned flags=Inertial)
  {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    double val = tp->second->pos[j];
    if (com_system and flags & Local) val -= com0[j];
    if (flags & Centered) val -= center[j];
    return val;
  }

  //! Get rectified component center
  std::vector<double> getCenter(unsigned flags=Inertial)
  {
    std::vector<double> ret(3, 0.0);
    for (int j=0; j<3; j++) {
      if (com_system and flags & Local) ret[j] += com0[j];
      if (flags & Centered) ret[j] += center[j];
    }

    return ret;
  }

  //! Access to velocities
  inline double Vel(int i, int j, unsigned flags=Inertial) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    double val = tp->second->vel[j];
    if (com_system and flags & Local) val -= cov0[j];
    return val;
  }
  
  //! Get positions
  inline void Pos(double *pos, int i, unsigned flags=Inertial) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    for (int k=0; k<3; k++) {
      pos[k] = tp->second->pos[k];
      if (com_system and flags & Local) pos[k] -= com0[k];
      if (flags & Centered) pos[k] -= center[k];
    }
  }

  //! Get velocities
  inline void Vel(double *vel, int i, unsigned flags=Inertial) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    for (int k=0; k<3; k++) {
      vel[k] = tp->second->vel[k];
      if (com_system and flags & Local) vel[k] -= cov0[k];
    }
  }
  
  //! Convert coordinates from inertial to local
  inline void ConvertPos(double* pos, unsigned flags=Local | Centered) {
    if (flags & Inertial) return;
    for (int k=0; k<3; k++) {
      if (com_system and flags & Local) pos[k] -= com0[k]; 
      if (flags & Centered) pos[k] -= center[k];
    }
  }
  
  //! Convert velocities from inertial to local
  inline void ConvertVel(double* vel, unsigned flags=Local | Centered) {
    if (flags & Inertial) return;
    for (int k=0; k<3; k++) {
      if (com_system and flags & Local) vel[k] -= cov0[k]; 
    }
  }
  //! Convert coordinates from local to inertial
  inline void PosConvert(double* pos, unsigned flags=Local | Centered) {
    if (flags & Inertial) return;
    for (int k=0; k<3; k++) {
      if (com_system and flags & Local) pos[k] += com0[k]; 
      if (flags & Centered) pos[k] += center[k];
    }
  }
  
  //! Convert velocities from local inertial
  inline void VelConvert(double* vel, unsigned flags=Local | Centered) {
    if (flags & Inertial) return;
    for (int k=0; k<3; k++) {
      if (com_system and flags & Local) vel[k] += cov0[k]; 
    }
  }

  //! Access to acceleration
  inline double Acc(int i, int j, unsigned flags=Inertial) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    double val = tp->second->acc[j];
    if (com_system and flags & Inertial) val += acc0[j];
    return val;
  }
  
  //! Add to position (by component)
  inline void AddPos(int i, int j, double val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    tp->second->pos[j] += val;
  }
  
  //! Add to position (by array)
  inline void AddPos(int i, double* val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    for (int k=0; k<3; k++) tp->second->pos[k] += val[k];
  }
  
  //! Add to position (by vector)
  inline void AddPos(int i, vector<double>& val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    for (int k=0; k<3; k++) tp->second->pos[k] += val[k];
  }
  
  //! Add to velocity (by component)
  inline void AddVel(int i, int j, double val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    tp->second->vel[j] += val;
  }
  
  //! Add to velocity (by array)
  inline void AddVel(int i, double* val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    for (int k=0; k<3; k++) tp->second->vel[k] += val[k];
  }
  
  //! Add to velocity (by vector)
  inline void AddVel(int i, vector<double>& val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    for (int k=0; k<3; k++) tp->second->vel[k] += val[k];
  }
  
  //! Add to accerlation (by component)
  inline void AddAcc(int i, int j, double val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    auto acc = getPseudoAccel(tp->second->pos, tp->second->vel);
    tp->second->acc[j] += val - acc[j];
  }
  
  //! Add to accerlation (by component)
  inline void AddAccExt(int i, int j, double val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    tp->second->acc[j] += val;
  }
  
  //! Add to acceleration (by array)
  inline void AddAcc(int i, double *val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    auto acc = getPseudoAccel(tp->second->pos, tp->second->vel);
    for (int k=0; k<3; k++)
      tp->second->acc[k] += val[k] - acc[k];
  }
  
  //! Add to acceleration (by array)
  inline void AddAccExt(int i, double *val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    for (int k=0; k<3; k++) tp->second->acc[k] += val[k];
  }
  
  //! Add to accerlation (by vector)
  inline void AddAcc(int i, vector<double>& val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    auto acc = getPseudoAccel(tp->second->pos, tp->second->vel);
    for (int k=0; k<3; k++)
      tp->second->acc[k] += val[k] - acc[k];
  }
  
  //! Add to accerlation (by vector)
  inline void AddAccExt(int i, vector<double>& val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    for (int k=0; k<3; k++) tp->second->acc[k] += val[k];
  }
  
  //! Add to potential
  inline void AddPot(int i, double val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    tp->second->pot += val;
  }
  
  //! Add to external potential
  inline void AddPotExt(int i, double val) {
    PartMap::iterator tp = particles.find(i);
    if (tp == particles.end()) {
      throw BadIndexException(i, particles.size(), __FILE__, __LINE__);
    }
    tp->second->potext += val;
  }
  
  //! Reset the level lists
  void reset_level_lists();

  //! Print out the level lists to stdout for diagnostic purposes
  void print_level_lists(double T);

  //! Running clock on the current potential/force evaluation
  Timer time_so_far;

  //! Time in potential/force computation so far
  double get_time_sofar() { return time_so_far.getTime(); }

  //! Check for indexing
  inline bool Indexing() { return indexing; }

  //! Check for multistep restriction (this is still a test option)
  inline bool NoSwitch() { return noswitch; }

  //! Only compute new levels on first step
  inline bool FreezeLev() { return freezeLev; }

  //! Compute level from minimum requested time step from last master step
  inline bool DTreset()  { return dtreset;  }

  //! Set H5 compression parameters; in principle, compression can be
  //! set per component, although it is set to be the same for all
  //! components currently.
  void setH5(int compress, int chunk)
  {
    H5compress = compress;
    H5chunk    = chunk;
  }

#if HAVE_LIBCUDA==1
  //@{
  //! CUDA utilities for handling host <===> device exchange

  //! Index of the device to use (-1 means ignore)
  int cudaDevice;

  typedef thrust::host_vector<cudaParticle> hostPartType;
  typedef thrust::host_vector<cudaParticle>::iterator hostPartItr;

  //! Cuda particles in host memory
  hostPartType host_particles;

  //! The cuda particles in CUDA device memory
  //@{

  //! Helper class to hold data and create/destroy the stream
  class cudaStreamData
  {
  private:
    static unsigned totalInstances;

  public:
    //! Instance id for debugging
    unsigned instance, id;

    //! Host particle beginning and ending iterator
    hostPartItr first, last;

    //! Device particles
    thrust::device_vector<cudaParticle> cuda_particles;

    //! Level list
    thrust::device_vector<int> levList;

    //! Level list pairs
    thrust::device_vector<thrust::pair<int, int>> levPair;

    //! Sorted indices
    thrust::device_vector<int> indx1, indx2;

    //! My stream
    cudaStream_t stream;

    //! Constructor
    cudaStreamData();

    //! Destructor
    ~cudaStreamData();
  };

  //! Initialize the stream(s)
  void cuda_initialize();

  typedef std::shared_ptr<cudaStreamData> cuSharedStream;
  cuSharedStream cuStream;
  //@}

  //! Copy particles host to device for stream entry
  void HostToDev(cuSharedStream cr);
  
  //! Copy particles device to host for stream entry
  void DevToHost(cuSharedStream cr);

  //! Copy particles from Component C to the cuda device structure
  void ParticlesToCuda(PartMap::iterator beg, PartMap::iterator fin);

  //! Convert all to device particles and copy to device
  void ParticlesToCuda()
  {
    // Copy from Particle to cudaParticle
    //
    ParticlesToCuda(particles.begin(), particles.end());

    // Push to device
    //
    HostToDev(cuStream);
  }

  //! Copy cuda device particles for this component back Component C
  void CudaToParticles(hostPartItr beg, hostPartItr end);

  //! Copy from device and convert to host particles
  void CudaToParticles()
  {
    DevToHost(cuStream);	// Pull from device
    CudaToParticles(cuStream->first, cuStream->last);
  }

  //! Make a new level list from copied particles
  void MakeLevlist();

  //! Sort the device vector by level and return level boundaries
  std::pair<unsigned, unsigned> CudaGetLevelRange(int lo, int hi);

  using PII   = std::pair<unsigned int, unsigned int>;
  using I1vec = std::vector<PII>;
  using I2vec = std::vector<I1vec>;

  //! Sort the device vector by level
  void CudaSortByLevel();

  //! Sort the device vector by level changes
  I2vec CudaSortLevelChanges();

  //! Zero the acceleration
  void ZeroPotAccel(int mlevel);

  //! Default bunch size for CUDA
  size_t bunchSize;

  //! Print out the level lists to stdout for diagnostic purposes
  void print_level_lists_cuda(double T);

  //! Get the level lists from cuda particles for diagnostic purposes
  std::vector<int> get_level_lists_cuda();

  //@{
  //! Time step vectors
  thrust::device_vector<int>   loLev, hiLev;
  thrust::device_vector<float> minDT, maxDT;
  //@}

  //@}
#endif

};

#endif
