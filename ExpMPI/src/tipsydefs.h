/* tipsydefs.h */

#define MAXDIM 3
#define forever for(;;)

typedef float Real;

//! Tipsy particle structure for a gas particle
struct gas_particle {
  //! Particle mass
    Real mass;
  //! Particle position vector
    Real pos[MAXDIM];
  //! Particle velociy vector
    Real vel[MAXDIM];
  //! Density
    Real rho;
  //! Temperature
    Real temp;
  //! Smoothing scale
    Real hsmooth;
  //! Metal value
    Real metals ;
  //! Gravitational potential
    Real phi ;
} ;

struct gas_particle *gas_particles;

//! Tipsy particle structure for a dark particle
struct dark_particle {
  //! Mass
    Real mass;
  //! Position vector
    Real pos[MAXDIM];
  //! Velocity vector
    Real vel[MAXDIM];
  //! Smoothing
    Real eps;
  //! Gravitational potential
    Real phi ;
} ;

struct dark_particle *dark_particles;

//! Tipsy particle structure for a star particle
struct star_particle {
  //! Particle mass
    Real mass;
  //! Particle position vector
    Real pos[MAXDIM];
  //! Particle velocity vector
    Real vel[MAXDIM];
  //! Metallicty value
    Real metals ;
  //! Formation time
    Real tform ;
  //! Smoothing
    Real eps;
  //! Gravitational potential
    Real phi ;
} ;


struct star_particle *star_particles;

//! Tipsy file header
struct dump {
  //! Dump time
    double time ;

  //! Total number of bodies
    int nbodies ;

  //! Dimensions
    int ndim ;

  //! Number of Gas particles
    int nsph ;

  //! Number of Dark particles
    int ndark ;

  //! Number of Star particles
    int nstar ;
} ;

struct dump header ;
