#ifndef _global_H
#define _global_H

#include <pthread.h>
#include <mpi.h>
#include <string>

#include <coef.H>

using namespace std;

/*! \file global.H
    \brief Global variables
    
    These variables can be changed/specified in the [global] stanza of
    the inputfile
*/

//! Maximum number of bodies per process
extern int nbodmax;

//! Maximum number of steps to execute
extern int nsteps;

//! Number of threads per process (e.g. one per processor)
extern int nthrds;

//! Number of steps between load balancing (use 0 for none)
extern int nbalance;

//! Load balancing threshold (larger difference initiates balancing)
extern double dbthresh;

//! Time step
extern double dtime;

//! Home directory for configuration files, etc.
extern string homedir;

//! Directory to look for loadable modules
extern string ldlibdir;

//! Input file (for restart)
extern string infile;

//! Parameter dump file
extern string parmfile;

//! Initial processor rate file
extern string ratefile;

//! Used for labeling report files
extern string runtag;

//! Toggle phase space advance (e.g. for use with externally supplied mapping).  On by default.
extern bool eqmotion;

				// Internal variables
//! Deprecated
extern bool use_cwd;

//! Global set on restart (to used by initializers and user modules)
extern bool restart;

//! Process priority set to NICE value
extern int NICE;

//! Times for each phase space slice (for Leap-Frog)
extern double tpos,tvel,tnow;

//! Current step number
extern int this_step;

//! Total mass in simulation (used by Orient and log diagnostics)
extern double mtot;

//! Global center of mass and center of velocity
extern double *gcom, *gcov; 

//! Resets total center of velocity to zero if true
extern bool global_cov;

//! Resets total accelerationto zero if true
extern bool fixacc;

				// MPI variables

//! Used to flag separate behavior on inintialization (unused)
extern int is_init;

//! Maintains number of processes, slaves, MPI id, etc.
extern int numprocs, slaves, myid, proc_namelen;

//! MPI returned processor name
extern char processor_name[MPI_MAX_PROCESSOR_NAME];

//! Internally defined slave communicator
extern MPI_Comm MPI_COMM_SLAVE;

//! Coefficient header on each output dump
// @{
extern CoefHeader coefheader;
extern struct CoefHeader2 coefheader2;
// @}

//! Theading variables
// @{
extern char threading_on;
extern pthread_mutex_t mem_lock;
// @}

// Particle components

#include <list>
#include <ComponentContainer.H>

//! One component container per simulation
extern ComponentContainer comp;

// External forces

#include <ExternalCollection.H>

//! One external force container per simulation
extern ExternalCollection external;

// Output list

#include <OutputContainer.H>
//! One output routine container per simulation
extern OutputContainer output;

// Parameter database

#include <ParamParseMPI.H>

//! Parameter database
extern ParamParseMPI *parse;

#endif
