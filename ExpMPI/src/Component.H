#ifndef _Component_H
#define _Component_H

#include <mpi.h>
#include <PotAccel.H>
#include <Orient.H>
#include <header.H>


//! Class associating the phase space with a biorthogonal basis
/*!
  The Component provides "interactors" for specific particle fields
  <br>
  Component file format:
  <br>
  <ol>
  <li> name [string in characters (whole line, spaces ok)]
  <li> force id [string]
  <li> file name [string] containing phase space to be read
  <li> biorthogonal [enum] followed by parameters necessary for contruction
  <li> [blank line]
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

*/
class Component
{

  friend class ScatterMFP;
  friend class OrbTrace;
  friend class ComponentContainer;

  static const int ndmax = 5;
  static const int nimax = 5;

  //! For interprocess communication
  struct Partstruct
  {				// Type:    Offsets:
    double mass;		// 0        0
    double pos[3];		// 1        1,2,3
    double vel[3];		// 2        4,5,6
    double acc[3];		// 3        7,8.9
    double pot;			// 4        10
    double potext;		// 5        11
    int    iatr[nimax];		// 6        12,13,14,15,16
    double datr[ndmax];		// 7        17,18,19,20,21
  };

private:

  //! Parallel distribute and particle io
  void read_bodies_and_distribute_ascii(void);
  void read_bodies_and_distribute_binary(istream *);
  void get_next_particle_from_file(Partstruct *onepart, istream *in);
  void part_to_Particle(Partstruct&, Particle&);
  void Particle_to_part(Partstruct&, Particle&);

  Partstruct *buf;
  static const int nbuf = 2000;
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

  // Initialize com tie list
  void get_com_component(string name);

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
  vector<double> rates, orates, trates;

  //! Binary read flag
  bool binary;

  //! Current phase-space time
  double time;
  
  //! Component list for tying center of mass
  list<Component*> com_tie;


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

  //! Use com or cov to center grid
  bool use_com, use_cov;
  /**
     These are pointers to the center of mass, center of velocity,
     angular momentum and the phase-space 6 vector.
  */
  double *com, *cov, *center, *angmom, *ps;

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
  //! Print verbose diagnostics
  bool EJdiag;
  //! Use particle KE for energy computation
  bool EJkinE;
  //! Use external potential for energy computation
  bool EJext;

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

  //! Update angular momentum values
  void get_angmom(void);

};

#endif
