/*
  Defines for EXPANSION code

  $Id$
*/

#define BSD_SOURCE

#undef MPE_PROFILE
/* #define MPE_PROFILE */

#include "mpi.h"
#ifdef MPE_PROFILE
#include "mpe.h"
#endif

#if !defined(__cplusplus)
/* #define _BSD_SOURCE */
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <values.h>
/* #include <malloc.h> */
#include <stdio.h>
#include "cutil.h"
#endif

#if defined(__cplusplus)
extern "C" {
#endif

#ifdef MAIN

				/* Numerical parameters (nearly unchanged) */
unsigned int seed = 3;
int nbodmax = 20000;
int lmax = 4;
int nmax = 10;
int mmax2 = 4;
int nmax2 = 10;
int nmaxx = 10;
int nmaxy = 10;
int nmaxz = 10;
int nminx = 0;
int nminy = 0;
int nminz = 0;
int nfft = 6;			/* DFT power of two for cylindrical basis */
int ncylordz = 6;
int ncylnzof = -1;
int ncylnx = 32;		/* Radial grid for table for PCA bases */
int ncylny = 32;		/* Vertical grid for table for PCA bases */
int ncylrecomp = 50;		/* Number of steps for PCA basis recomp */
int ncylkeep = 100;		/* Number in average for center and axis */
int nsteps = 500;
int nscale = 20;
int nthrds = 2;			/* Number of POSIX threads */
double dtime = 0.1;
int nlog = 10;
int ncoef = 1;
int ncoefcyl = 1;
int nrelx = 1;
int nlist = 50;
int nscat = 20;
int ntipsy = 50;
int ndiag = 50;
int nchkpt = 100;		/* Number of steps bewteen checkpoints */
int ncylamom = 500;
int EJcyl = 0;			/* Use orientation calc for axis and com */
double rmax = -1.0;
double zmax = 20.0;		/* Vertical cutoff for cylindrical basis */
double ecyl0 = -0.5;		/* Initial energy cutoff for L-plane */
double rdiag = 10.0;
double scale = 1.0;
double rcylEM = 20.0;		/* Radial cutoff for empirical orthgonal */
				/* cylindrical basis */
double rcylSL = 100.0;		/* Radial cutoff for Sturm-Lioville grid */
				/* cylindrical basis */
double rsphSL = 100.0;		/* Radial cutoff for Sturm-Lioville grid */
				/* spherical basis */
double acyl = 1.0;		/* Radial scaling for cylindrical basis */
double hcyl = 1.0;		/* Vertical scaling for cylindrical basis */
double hslab = 0.1;		/* Slab scale height */
double rmax_tidal=1.e10;	/* radius for freezing particles. KL 5/27/92*/
double hills_omega = 0.0;	/* frequency of tidal field. KL 5/27/92 */
double hills_p = 1.0;		/* log derivative of tidal field. default is 
				constant rotation curve. KL 5/27/92 */

double tksmooth = 1.0;		/* Choose 1.0 for original Hall selector */
double tkcum = 0.95;		/* Cutoff for cumulative var */

double tauscat = 1.0;		/* rho*v*L for scattering */
			
				/* variables for mass loss */
double *initial_mass; 
double frac_mloss=0.0; 
double t_mloss=1.0;

				/* Load balancing tolerance */
double tbalance=0.1;

				/* Flags */
int bessel_sph = 1;
int c_brock = 0;
int c_brock_disk = 0;
int hernq = 0;
int sphereSL = 0;
int cube = 0;
int slab = 0;
int slabSL = 1;
int cylinder = 0;
int nulltest = 0;
int fixpos = 1;			/* 1=set origin for all, 2=for each comp'nt */
int fixvel = 1;
int fixacc = 1;
int selfgrav = 1;
int NO_L1 = 0;			/* 1 to suppress L=1 terms in force (only!) */
int adjust = 0;
int outbods = 1;
int outengy = 0;
int tipsy = 1;
int olist = 0;
int diag = 0;
int coef = 0;
int coefcyl = 0;
int relx = 0;
int scatter = 0;
int finish = 0;
int tides = 0;			/* is there a tidal force? KL 5/27/92 */
int self_consistent = 1;	/* do we update pot at each step? KL 3/20/92 */
int inertia = 0;		/* print inertia tensor?  KL 6/23/92 */
int shock = 0;			/* external shocking routine */
int halo = 0;			/* add spherical halo routine */
int disk_on_halo = 1;		/* effect of component 2 on component 1 */
int selector = 0;		/* use MISE to select coefficients */
				/* non-zero for Tarter & Kronmel */
int tk_type = 0;		/* 0=Hall, 1=TK, 2=cumulative variance */
int eigen_tk = 1;		/* Use maximum variance decomposition */
int npca = 10;			/* Number of steps per weight recomputation */
int pcaout = 0;			/* print diagnostics for pca decomposition */
int npcaout = 10;		/* Number of steps per diagnostic output */
int L_pca = 0;			/* L value for pca diag output */
int M_pca = 0;			/* M value for pca diag output */
int balance = 0;		/* Load balancing? */
int nbalance = 20;		/* Number of steps to check balance */
int NICE=0;			/* Nice value on remote machines */
int zerocom = 1;		/* Set com and com velocity for each comp */
int zerovel = 1;		/* initially */

int user = 0;			/* external perturbation to be added */

double pmsoft = 0.1;		/* Softening for point mass perturbers */

				/* User definable parameters */
double U1 = 0.0;
double U2 = 0.0;
double U3 = 0.0;
double U4 = 0.0;
double U5 = 0.0;
double U6 = 0.0;

				/* Files */
char homedir[100] = "./";
char logfile[100] = "LOG.FILE";
char infile[100] = "IN.FILE";
char outname[100] = "OUT";
char parmfile[100] = "PARMS.FILE";
char coeffile[100] = "COEF.FILE";
char coeffilecyl[100] = "COEFCYL.FILE";


				/* Global variables */
double *mass,*x,*y,*z,*rr,*vx,*vy,*vz,*ax,*ay,*az,*pot,*potext,*size;
double *esave, *mfp;
int *component;
double tpos,tvel,tnow;
double **expcoef, **expcoef1, **bb, ***cc, ***cc1, **ss, **normM;
double *com1, *com2;
int pmnum0, pmnum, *pmlist;
int this_step,nbodies,used,ninteract=0,nmerge=0,dof=3;
int restart;
int ncompcyl;
double cylmass,cyltime;
double zcm_slab=0.0, vzcm_slab=0.0, azcm_slab=0.0;

				/* MPI variables */
int is_init=1;
int numprocs, slaves, myid, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

int *nbodies_index;		/* final particle indices for each bunch */
int *nbodies_table;
double *rates, *trates, *orates;
int *rate_index;
MPI_Comm MPI_COMM_SLAVE;

#else

extern unsigned int seed;
extern int nbodmax;
extern int lmax;
extern int nmax;
extern int mmax2;
extern int nmax2;
extern int nmaxx;
extern int nmaxy;
extern int nmaxz;
extern int nminx;
extern int nminy;
extern int nminz;
extern int nfft;
extern int ncylnzof;
extern int ncylordz;
extern int ncylrecomp;
extern int ncylkeep;
extern int ncylnx;
extern int ncylny;
extern int nsteps;
extern int nscale;
extern int nthrds;
extern double dtime;
extern int nlog;
extern int ncoef;
extern int ncoefcyl;
extern int nrelx;
extern int nlist;
extern int nscat;
extern int ntipsy;
extern int ndiag;
extern int nchkpt;
extern int ncylamom;
extern int EJcyl;
extern double rmax;
extern double zmax;
extern double ecyl0;
extern double rdiag;
extern double scale;
extern double rcylEM;
extern double rcylSL;
extern double rsphSL;
extern double acyl;
extern double hcyl;
extern double hslab;
extern double rmax_tidal;	/* radius for freezing particles. KL 5/27/92*/
extern double hills_omega;	/* frequency of tidal field. KL 5/27/92 */
extern double hills_p;		/* log derivative of tidal field. KL 5/27/92 */
				/* default is constant rotation curve */

extern double tksmooth;		/* Choose 1.0 for original Hall selector */
extern double tkcum;

extern double tauscat;		/* rho*v*L for scattering */


				/* variables for mass loss */
extern double *initial_mass; 
extern double t_mloss; 
extern double frac_mloss;
extern double tbalance;

extern int bessel_sph;
extern int c_brock;
extern int c_brock_disk;
extern int hernq;
extern int sphereSL;
extern int cube;
extern int slab;
extern int slabSL;
extern int cylinder;
extern int nulltest;
extern int fixpos;
extern int fixvel;
extern int fixacc;
extern int adjust;
extern int selfgrav;
extern int NO_L1;
extern int outbods;
extern int outengy;
extern int tipsy;
extern int olist;
extern int diag;
extern int coef;
extern int coefcyl;
extern int relx;
extern int scatter;
extern int finish;
extern int self_consistent;		/* KL 3/20/92 */
extern int tides;			/* tidal field? KL 5/27/92 */
extern int inertia;			/* print inertia tensor? KL 6/23/92 */
extern int shock;
extern int halo;
extern int disk_on_halo;
extern int selector;
extern int tk_type;
extern int eigen_tk;
extern int npca;
extern int pcaout;
extern int npcaout;
extern int L_pca;
extern int M_pca;
extern int balance;
extern int nbalance;
extern int NICE;
extern int zerocom;
extern int zerovel;

extern int user;
extern double pmsoft;
extern double U1, U2, U3, U4, U5, U6;

extern char homedir[];
extern char logfile[];
extern char infile[];
extern char outname[];
extern char parmfile[];
extern char coeffile[];
extern char coeffilecyl[];

extern double *mass,*x,*y,*z,*rr,*vx,*vy,*vz,*ax,*ay,*az,*pot,*potext,*size;
extern double *esave, *mfp;
extern int *component;
extern double tpos,tvel,tnow;
extern double **expcoef, **expcoef1, **bb, ***cc, ***cc1, **ss, **normM;
extern double *com1, *com2;
extern int pmnum0, pmnum, *pmlist;

extern int this_step,nbodies,used,ninteract,nmerge,dof;
extern int restart;
extern int ncompcyl;
extern double cylmass,cyltime;
extern double zcm_slab, vzcm_slab, azcm_slab;

				/* MPI variables */
extern int is_init;
extern int numprocs, slaves, myid, proc_namelen;
extern char* processor_name;
extern int *nbodies_index;
extern int *nbodies_table;
extern double *rates, *trates, *orates;
extern int *rate_index;
extern MPI_Comm MPI_COMM_SLAVE;


#endif /* MAIN */

				/* Function declarations */

void read_bodies_and_init(void);
void compute_potential(void);
void fix_positions(void);
void fix_positions_by_component(void);
void fix_acceleration(void);
void init_velocity(void);
void begin_run(void);
void out_list(int);
void out_tipsy(int);
void out_chkpt(void);
void out_log(int);
void out_coef(int);
void out_coef_cyl(int);
void incr_position(void);
void incr_velocity(void);
void get_acceleration_and_potential(void);
void bessacp(void);
void set_radius(double);
void tidal_field(void);		/* tidal acceleration. KL 5/27/92 */
int freeze_particle(int);	/* freeze particle? KL 5/27/92 */
void mass_loss(void);		/* reduce masses. KL 7/3/92 */
void out_pca(int);
void out_relx(int);
void distribute_relx(void);
void distribute_mfp(void);
void get_acceleration_and_potential_pointmass();
void make_pointmass_list();


#ifndef MAX
#define MAX(A,B) (A>B ? A : B)
#endif
#ifndef MIN
#define MIN(A,B) (A<B ? A : B)
#endif

/*****Stuff for expand.c*****/

#define RNUM 1000

struct RGRID {
  double **rw;
  double **rw2;
  int nmax;
};

#ifdef MAIN
struct RGRID *dens_grid, *potl_grid;
double *r_grid, r_grid_del;
#else
extern struct RGRID *dens_grid, *potl_grid;
extern double *r_grid, r_grid_del;
#endif

double get_dens(double r, int l, double *coef);
double densi(double rr, int l, int n),potli(double rr, int l, int n);
void make_grid(double rmin, double rmax, int lmax, int nmax);
void get_potacc(double r, int l, double *coef, double *p, double *dp);



/******Mathematical utilities******/
void splint();
void spline();
double zbrent();
double dgammln();
double *sbessjz();
double sbessj();
void locate();
double plgndr(int, int, double);
double dplgndr(int l, int m, double x);
double factrl(int n);

/* Constants */

#define DSMALL  1.0e-8
#define DSMALL2 1.0e-8

/* parallel prototypes */

void setup_distribution(void);
void distribute_particles(void);
void gather_particles(void);
void test_mpi(void);
void parallel_gather_coefficients(void);
void parallel_distribute_coefficients(void);
void recompute_processor_rates(int *);
void MPL_reset_timer(void);
void MPL_start_timer(void);
void MPL_stop_timer(void);
double MPL_read_timer(int reset);

#if defined(__cplusplus)
}
#endif

