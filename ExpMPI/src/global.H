#ifndef _global_H
#define _global_H

#include <pthread.h>
#include <mpi.h>

#include <vector>
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

//! Directory for output
extern string outdir;

//! Used for labeling report files
extern string runtag;

//! Toggle phase space advance (e.g. for use with externally supplied mapping).  On by default.
extern bool eqmotion;

//! Multistep levels (default: 0 means to multistepping)
extern unsigned multistep;

//! Synchronize multistep positions at each substep (default: true)
extern int posnsync;

//! Fraction of dynamical time for determining multistep level
//@{
//! Velocity coefficient (default: 0.01)
extern double dynfracV;
//! Acceleration coefficient (default: 0.03)
extern double dynfracA;
//@}


				// Internal variables

//! Global set on restart (to used by initializers and user modules)
extern bool restart;

//! Use Node 0's home dir for the working dir on all nodes
extern bool use_cwd;

//! Process priority set to NICE value
extern int NICE;

//! Output logging level
extern int VERBOSE;

//! Times for each phase space slice (for Leap-Frog)
extern double tnow;

//! Current step number
extern int this_step;

//! Total mass in simulation (used by Orient and log diagnostics)
extern double mtot;

//! Global center of mass
extern double *gcom;

//! Global center of velocity
extern double *gcov; 

//! Resets total center of velocity to zero if true
extern bool global_cov;

//! Maximum number of multistep steps
extern int Mstep;

//! Current point in multistepping
extern int mstep;

//! Multistep time step spacing
extern vector<int> mintvl;

//! Lowest synchronized level at a step
extern vector<int> mfirst;

//! Particles at each level
extern vector<int> levpop;

//! Sub steps for last coefficient evaluation
extern vector<int> stepL;

//! Sub steps for current coefficient evaluation
extern vector<int> stepN;

//! Multistep level flag: levels currently synchronized
extern vector< vector<bool> > mactive;

/// Helper class to pass info for incr_postion and incr_velocity
struct thrd_pass_posvel 
{
  //! Time step
  double dt;

  //! Levels flag
  int mlevel;

  //! Thread counter id
  int id;
};


//! Multithreding data structures for incr_position and incr_velocity
extern vector<thrd_pass_posvel> posvel_data;

//! Threads for incr_position and incr_velocity
extern vector<pthread_t> posvel_thrd;


/**
   Phase space dump on next step

   This flag is set by signalling the root process with SIGHUP
*/
extern unsigned char dump_signal;

/**
   Stop on next step

   This flag is set by signalling the root process with SIGTERM
*/

extern unsigned char stop_signal;

				// MPI variables

//! Used to flag separate behavior on inintialization (unused)
extern int is_init;

#ifndef _localmpi_h
//! Maintains number of processes, slaves, MPI id, etc.
//@{
extern int numprocs, slaves, myid, proc_namelen;
//@}

//! MPI returned processor name
extern char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif

//! Internally defined slave communicator
extern MPI_Comm MPI_COMM_SLAVE;

/** @name Coefficient header on each output dump */
//@{
//! Header type 1
extern CoefHeader coefheader;
//! Header type 2
extern struct CoefHeader2 coefheader2;
//@}

/** @name Theading variables */
//@{
//! signal that threading is in use
extern char threading_on;
//! mutex for memory
extern pthread_mutex_t mem_lock;
//@}

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
