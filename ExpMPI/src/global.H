#ifndef _global_H
#define _global_H

#include <pthread.h>
#include <mpi.h>
#include <string>

#include <coef.H>

//! Maximum number of bodies per process
extern int nbodmax;
//! Maximum number of steps to execute
extern int nsteps;
//! Number of threads per process (e.g. one per processor)
extern int nthrds;
//! Time step
extern double dtime;
//! Ignore particles beyond this radius
extern double rmax_tidal;

//! Home directory for configuration files, etc.
extern string homedir;

//! Directory to look for loadable modules
extern string ldlibdir;


//! Input file (for restart)
extern string infile;

//! Parameter dump file
extern string parmfile;

extern bool use_cwd;

extern bool restart;

extern int NICE;

extern double tpos,tvel,tnow;
extern int this_step;

extern double mtot;
extern double *gcom, *gcov; 
extern bool global_cov;
extern bool fixacc;

				// MPI variables
extern int is_init;
extern int numprocs, slaves, myid, proc_namelen;
extern char processor_name[MPI_MAX_PROCESSOR_NAME];
extern MPI_Comm MPI_COMM_SLAVE;

extern CoefHeader coefheader;
extern struct CoefHeader2 coefheader2;

extern char threading_on;
extern pthread_mutex_t mem_lock;

// Particle components

#include <list>
#include <ComponentContainer.H>

extern ComponentContainer comp;

// External forces

#include <ExternalCollection.H>

extern ExternalCollection external;

// Output list

#include <OutputContainer.H>
extern OutputContainer output;

// Parameter database

#include <ParamParseMPI.H>

extern ParamParseMPI *parse;

#endif
