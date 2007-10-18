/* tipsydefs.h */

#define MAXDIM 3

typedef float Real;

struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals ;
    Real phi ;
} ;

extern struct gas_particle *gas_particles;

struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
} ;

extern struct dark_particle *dark_particles;

struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
} ;

extern struct star_particle *star_particles;

struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;

extern struct dump header ;

void xdr_gas();
void xdr_dark();
void xdr_star();
int  xdr_init();
void xdr_quit();

