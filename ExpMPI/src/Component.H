#ifndef _Component_H
#define _Component_H

#include <mpi.h>
#include <PotAccel.H>
#include <Orient.H>
#include <header.H>
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

  @param ton is the midpoint of the adiabatic turn on for this
  component

  @param twid is the width for the adiabatic turn on.  Unless
  <code>ton</code> or <code>toff</code> is specified, the component
  begins and remains fully on.

*/
class Component
{

  friend class ScatterMFP;
  friend class OrbTrace;
  friend class ComponentContainer;

  static const int ndmax = 5;
  static const int nimax = 5;

public:

  //! For interprocess communication
  struct Partstruct
  {				///< Type:    Offsets:
    double mass;		///< 0        0
    double pos[3];		///< 1        1,2,3
    double vel[3];		///< 2        4,5,6
    double acc[3];		///< 3        7,8.9
    double pot;			///< 4        10
    double potext;		///< 5        11
    int    iatr[nimax];		///< 6        12,13,14,15,16
    double datr[ndmax];		///< 7        17,18,19,20,21
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

  bool seq_ok;
  int seq_beg, seq_end;

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

  // For load balancing
  vector <loadb_datum> loadb;
  void append_particles(int from, int to, int begin, int number);
  void insert_particles(int from, int to, int begin, int number, int loc);

  // Compute initial com position and velocity from phase space
  void initialize_com_system();


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
  vector<int> plist, ncount, nbodies_index, nbodies_table;
  //! Vectors for holding the rates per node for load balancing
  vector<double> orates, trates;

  //! Binary read flag
  bool binary;

  //! COM log file
  string comfile;

public:

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
  double *com, *cov, *center, *angmom, *ps;

  /**
     Center of mass system coordinates
  */
  double mtot0;
  double *com0, *cov0, *acc0;
  double *comI, *covI;

  //! Total mass
  double mtot;

  //! Dimension of the phase space
  int dim;

  //! Bodies on all nodes
  int nbodies_tot;

  //! Points to associated force instance
  PotAccel *force;

  //! Points to associated particles
  vector<Particle> particles;

  //! Attribute array dimensions
  int niattrib, ndattrib;

  //! Binding energy orientation flag
  int EJ;

  //! Number to keep in orientation stack
  int nEJkeep;

  //! Target number of particles for average
  int nEJwant;

  //! Initial guess for energy cutoff
  double eEJ0;

  //! Initial EJ center
  double EJx0, EJy0, EJz0;

  //! EJ center
  double *EJcen;

  //! Initial EJ velocity
  double EJu0, EJv0, EJw0;

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
  
  //! Write binary component phase-space structure
  void write_binary(ostream *out);
  
  //! Write ascii component phase-space structure
  void write_ascii(ostream *out, bool accel = false);
  
  //! Redestribute this component
  void redistribute_particles(void);

  //! Compute center of mass and center of velocity
  void fix_positions(void);

  //! Subtract mean acceleration from each particle accel
  void update_accel(void);

  //! Update angular momentum values
  void get_angmom(void);

  //! Adiabatic turn on factor, range in [0, 1]
  double Component::Adiabatic(void);

  //! Is particle out of bounds?
  bool freeze(const Particle&);

};

#endif
