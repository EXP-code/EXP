#ifndef _Component_H
#define _Component_H

#include <mpi.h>
#include <PotAccel.H>
#include <Orient.H>
#include <header.H>
#include <localmpi.h>
#include <utility>


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

  \par c_brock 
  Clutton-Brock's spherical basis, see CBrock

  \par c_brock_disk 
  Clutton-Brock's two-dimensional disk basis, see CBrockDisk

  \par hernq 
  %Hernquist's spherical basis, see Hernquist

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

  @param eEJ0 is the initial guess for the largest energy particle to
  keep in the ensemble

  @param EJx0 is the initial EJ center x-coordinate (default: 0)

  @param EJy0 is the initial EJ center y-coordinate (default: 0)

  @param EJz0 is the initial EJ center z-coordinate (default: 0)

  @param EJu0 is the initial EJ center u-velocity   (default: 0)

  @param EJv0 is the initial EJ center v-velocity   (default: 0)

  @param EJw0 is the initial EJ center w-velocity   (default: 0)

  @param EJkinE set true includes kinetic energy in the energy
  determination, otherwise potential energy alone is used

  @param EJext set true includes the external potential energy in the
  energy computation

  @param EJdiag turns on verbose EJ reporting

  @param EJdryrun set to true does not apply EJ centering or rotation but
  computes the transformation

  @param EJlinear set true moves center with constant velocity.  For
  testing only!! (default: false)

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

  </ol>
  

*/
class Component
{

  friend class ScatterMFP;
  friend class OrbTrace;
  friend class ComponentContainer;

public:

  //! Maximum number of *extra* integer components
  const static int nimax = 5;
  //! Maximum number of *extra* double components
  const static int ndmax = 12;

  //! For interprocess communication, need static structure for MPI mapping
  struct Partstruct
  {
				///< Type:    Offsets:
    double mass;		///< 0        0
    double pos[3];		///< 1        1,2,3
    double vel[3];		///< 2        4,5,6
    double acc[3];		///< 3        7,8.9
    double pot;			///< 4        10
    double potext;		///< 5        11
    unsigned level;		///< 6        12
    unsigned long indx;		///< 7        13
    int    iatr[nimax];		///< 8        14,15,16,16,18 (e.g.)
    double datr[ndmax];		///< 9        19,20,21,22,23 (e.g.)
  };

private:

  //! Parallel distribute and particle io
  void read_bodies_and_distribute_ascii(void);
  void read_bodies_and_distribute_binary(istream *);
  void get_next_particle_from_file(Partstruct *onepart, istream *in);
  void part_to_Particle(Partstruct&, Particle&);
  void Particle_to_part(Partstruct&, Particle&);
  void load_balance(void);

  Partstruct *buf;
  static const int nbuf = 2000;

  bool seq_check, ordered;
  int seq_beg, seq_end, seq_cur;

				// Used in computing target rates for
				// load balancing
  int round_up(double dnumb);

  //! Error
  void bomb(const string& msg);

  // Particle handler variables
  bool npart;
  Partstruct *npart_p;
  int npart_tot, npart_cur;

  // Adiabatic mass
  bool adiabatic;
  double ton, toff, twid;

  // For particle retrival in the proper order
  Partstruct * get_particles_ordered(int* number);
  Partstruct * get_particles_unordered(int* number);

  // For load balancing
  vector <loadb_datum> loadb;
  void append_particles(int from, int to, int begin, int number);
  void insert_particles(int from, int to, int begin, int number, int loc);

  // Compute initial com position and velocity from phase space
  void initialize_com_system();

  // Momentum tracking
  bool consp;
  int tidal;

protected:

  /** Initialize parameters (called by ComponentContainer after component
      instantiation */
  void initialize(void);

  //! Slab parameters
  double azcm_slab;

  //! Bodies on this node
  int nbodies;

  /** Used by gather and distribution routines to define a particle
      structure to MPI */
  MPI_Datatype Particletype;

  //! Setup distribution vectors
  void setup_distribution(void);

  //! Vectors for holding current particle partition
  //@{
  vector<int> plist, ncount, nbodies_index, nbodies_table;
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

  //! Points to associated particles
  vector<Particle> particles;

public:

  //! Describe the phase space coordinates
  enum PSTYPE {Inertial=1, Local=2, Centered=4};

  //! Particle file
  string pfile;

  //! Identify this component
  string name;

  //! Force method
  string id;

  //! Parameter string for component
  string cparam;

  //! Parameter string for force
  string fparam;

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
  

  //! Maximum initial particle radious
  double rmax;

  /**
     These are pointers to the center of mass, center of velocity,
     angular momentum and the phase-space 6 vector.
  */
  //@{
  //! Center of mass
  double *com;
  //! Center of velocity
  double *cov;
  //! Center (e.g. set by Orient)
  double *center;
  //! Angular momentum vector
  double *angmom;
  //! Phase space vector
  double *ps;
  //@}

  /**
     @name Center of mass system coordinates
  */
  //@{
  //! System mass, center of mass, center of velocity and acceleration
  double mtot0, *com0, *cov0, *acc0;
  //! Initial center of mass and velocity
  double *comI, *covI;
  //! Maximum radius for contribution to center-of-mass computation
  double rcom;

  //@}

  //! Total mass
  double mtot;

  //! Dimension of the phase space
  int dim;

  //! Bodies on all nodes
  int nbodies_tot;

  /** Particle list per level
      This will be a sorted list for each level.  We use a vector rather
      than (e.g.) a set for the random insertion
  */
  vector< vector<int> > levlist;

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

  //! Initial guess for energy cutoff
  double eEJ0;

  //! Initial EJ center
  //@{
  //! x-coord
  double EJx0;
  //! y-coord
  double EJy0;
  //! z-coord
  double EJz0;
  //@}

  //! EJ center
  double *EJcen;

  //! Initial EJ velocity
  //@{
  //! x-coord
  double EJu0;
  //! y-coord
  double EJv0;
  //! z-coord
  double EJw0;
  //@}

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

  //! Use center of mass coordinate system
  bool com_system;

  //! Log the center of mass coordinate system at each step
  bool com_log;

  //! Maximum radius from center for force evaluation
  double rtrunc;

  //! Orientation cache
  Orient *orient;

  //! The constructor
  /*! 
    Takes the id string id and phase space file name as arguments
  */
  Component(string NAME, string ID, string CPARAM, string PS, string FPARAM);

  //! Initialize from file
  Component(istream *in);

  //! Destructor
  ~Component();

  //! Parallel gather
  Partstruct* get_particles(int* number);
  
  //! Particle handler
  int get_next_particle(Partstruct *onepart);
  
  //! Retrieve particle count on a particular node (one node at a time)
  int particle_count(int node) { 
    if (node<0 || node>=numprocs) return 0;
    return ncount[node]; 
  }

  //! Retrieve particle count on a particular node (all nodes at once)
  vector<int> particle_count() { return ncount; }

  //! Write binary component phase-space structure
  void write_binary(ostream *out);
  
  //! Write ascii component phase-space structure
  void write_ascii(ostream *out, bool accel = false);
  
  //! Redestribute this component
  void redistribute_particles(void);

  //! Compute center of mass and center of velocity
  void fix_positions(void);

  //! Control: assume particle on nodes remain unchanged
  void setOrdered() { ordered = true; }

  //! Control: assume particle on nodes may change
  void setUnordered() { ordered = false; }

  /** Redistribute particle by list

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
  void get_angmom(void);

  //! Adiabatic turn on factor, range in [0, 1]
  double Adiabatic(void);

  //! Is particle out of bounds?
  bool freeze(const Particle&);

  //! Is particle beyond COM system
  bool escape_com(const Particle&);

  //! Access to particle
  Particle *Part(int i) {
    return &particles[i];
  }

  //! Particle vector size
  unsigned Number() {
    return particles.size();
  }

  //! Access to mass
  inline double Mass(int i) {
    return particles[i].mass;
  }

  //! Access to positions
  inline double Pos(int i, int j, unsigned flags=Inertial) {
    double val = particles[i].pos[j];
    if (com_system) {
      if (flags & Local)
	val -= comI[j];
      else if (flags & Inertial) 
	val += com0[j] - comI[j];
    }
    if (flags & Centered) val -= center[j];
    return val;
  }

  //! Access to velocities
  inline double Vel(int i, int j, unsigned flags=Inertial) {
    double val = particles[i].vel[j];
    if (com_system) {
      if (flags & Local)
	val -= covI[j];
      else if (flags & Inertial) 
	val += cov0[j] - covI[j];
    }
    return val;
  }
  
  //! Get positions
  inline void Pos(double *pos, int i, unsigned flags=Inertial) {
    for (int k=0; k<3; k++) {
      pos[k] = particles[i].pos[k];
      if (com_system) {
	if (flags & Local) 
	  pos[k] -= comI[k];
	else if (flags & Inertial) 
	  pos[k] += com0[k] - comI[k];
      }
      if (flags & Centered) pos[k] -= center[k];
    }
  }

  //! Get velocities
  inline void Vel(double *vel, int i, unsigned flags=Inertial) {
    for (int k=0; k<3; k++) {
      vel[k] = particles[i].vel[k];
      if (com_system) {
	if (flags & Local)
	  vel[k] -= covI[k];
	else if (flags & Inertial) 
	  vel[k] += cov0[k] - covI[k];
      }
    }
  }
  
  //! Convert coordinates from inertial to local
  inline void ConvertPos(double* pos, unsigned flags=Local | Centered) {
    if (flags & Inertial) return;
    for (int k=0; k<3; k++) {
      if (com_system) {
	if (flags & Local)
	  pos[k] -= com0[k]; 
      }
      if (flags & Centered) pos[k] -= center[k];
    }
  }

  //! Convert velocities from inertial to local
  inline void ConvertVel(double* vel, unsigned flags=Local | Centered) {
    if (flags & Inertial) return;
    for (int k=0; k<3; k++) {
      if (com_system) {	
	if (flags & Local) 
	  vel[k] -= cov0[k]; 

      }
    }
  }

  //! Access to acceleration
  inline double Acc(int i, int j, unsigned flags=Inertial) {
    double val = particles[i].acc[j];
    if (com_system && (flags & Inertial)) val += acc0[j];
    return val;
  }
  
  //! Add to position (by component)
  inline void AddPos(int i, int j, double val) {
    particles[i].pos[j] += val;
  }
  
  //! Add to position (by array)
  inline void AddPos(int i, double* val) {
    for (int k=0; k<3; k++) particles[i].pos[k] += val[k];
  }
  
  //! Add to position (by vector)
  inline void AddPos(int i, vector<double>& val) {
    for (int k=0; k<3; k++) particles[i].pos[k] += val[k];
  }
  
  //! Add to velocity (by component)
  inline void AddVel(int i, int j, double val) {
    particles[i].vel[j] += val;
  }
  
  //! Add to velocity (by array)
  inline void AddVel(int i, double* val) {
    for (int k=0; k<3; k++) particles[i].vel[k] += val[k];
  }
  
  //! Add to velocity (by vector)
  inline void AddVel(int i, vector<double>& val) {
    for (int k=0; k<3; k++) particles[i].vel[k] += val[k];
  }
  
  //! Add to accerlation (by component)
  inline void AddAcc(int i, int j, double val) {
    particles[i].acc[j] += val;
  }
  
  //! Add to acceleration (by array)
  inline void AddAcc(int i, double *val) {
    for (int k=0; k<3; k++) particles[i].acc[k] += val[k];
  }
  
  //! Add to accerlation (by vector)
  inline void AddAcc(int i, vector<double>& val) {
    for (int k=0; k<3; k++) particles[i].acc[k] += val[k];
  }
  
  //! Add to potential
  inline void AddPot(int i, double val) {
    particles[i].pot += val;
  }
  
  //! Add to external potential
  inline void AddPotExt(int i, double val) {
    particles[i].potext += val;
  }
  
  //! Insert particle n into the level list at lev
  inline void LevInsert(int n, int lev) {
    vector<int>::iterator p = 
      lower_bound(levlist[lev].begin(), levlist[lev].end(), n);
    levlist[lev].insert(p, n);
  }

  //! Remove particle n from the level list at lev
  inline void LevRemove(int n, int lev) {
    vector<int>::iterator p = 
      find(levlist[lev].begin(), levlist[lev].end(), n);
    if (p!=levlist[lev].end()) levlist[lev].erase(p);
    else cerr << "Process " << myid << ": particle not found, n=" << n
	      << " for level=" << lev << endl;
  }

};

#endif
