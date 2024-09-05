#ifndef _global_H
#define _global_H

#include <pthread.h>
#include <mpi.h>

#include <vector>
#include <string>
#include <memory>
#include <random>

#include <list>
#include <coef.H>

#include <chkTimer.H>
#include <BarrierWrapper.H>

#include <yaml-cpp/yaml.h>

#include <localmpi.H>

using namespace std;

/*! \file global.H
    \brief Global variables
    
    These variables can be changed/specified in the [global] stanza of
    the inputfile
*/

//! Install error handler for MPI debugging
extern bool mpi_wait;

//! Install error handler for fpe debugging
extern bool fpe_trap;

//! Install error handler for fpe backtrace
extern bool fpe_trace;

//! Use gdb or native trace handling
extern bool gdb_trace;

//! Install error handler for fpe wait
extern bool fpe_wait;

//! Enable traceback on exceptions
extern bool traceback;

//! Maximum number of steps to execute
extern int nsteps;

//! Number of gpu devices per node
extern int ngpus;

//! Number of steps between particle number reports (use 0 for none)
extern int nreport;

//! Number of steps between load balancing (use 0 for none)
extern int nbalance;

//! Load balancing threshold (larger difference initiates balancing)
extern double dbthresh;

//! Particle ferry buffer size
extern unsigned PFbufsz;

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

//! Toggle phase space advance (e.g. for use with externally supplied mapping).  On by default.
extern bool eqmotion;

//! Constrain level changes per step (default: 0 means no constraint)
extern unsigned shiftlevl;

//! Lowest level for centering recomputation (default: multistep/2)
extern int centerlevl;

//! Fraction of dynamical time for determining multistep level
//@{

//! Scale coefficient (default: 1.00)
extern double dynfracS;

//! Drift coefficient (default: 1000.0)
extern double dynfracD;

//! Velocity coefficient (default: 0.01)
extern double dynfracV;

//! Acceleration coefficient (default: 0.03)
extern double dynfracA;

//! Escape time scale coefficient (default: 0.05)
extern double dynfracP;

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

//! Step timing
extern bool step_timing;

//! Times for each phase-space time slice (for Leap-Frog)
extern double tnow;

//! Last phase space dump step
extern int psdump;

//! Maximum fraction of particles that request a smaller than minimum timestep
extern double max_mindt;

//! Total mass in simulation (used by Orient and log diagnostics)
extern double mtot;

//! Global center of mass
extern double *gcom;

//! Global center of velocity
extern double *gcov; 

//! Resets total center of velocity to zero if true
extern bool global_cov;

//! True when initializing before first time step
extern bool initializing;

//! Current Nbody step id
extern int this_step;

//! Current point in multistepping for drift step
extern int mdrft;

//! Multistep sublevel info used for nbody
extern unsigned maxlev;

//! Current point in multistepping for velocity step
extern int mstep;

//! Maximum number of multistep steps
extern int Mstep;

//! Multistep time step spacing
extern vector<int> mintvl;

//! Lowest synchronized level at a step
extern vector<int> mfirst;

//! Multistep interpolation array
extern vector< vector<int> > dstepL, dstepN;

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

//! Suppress parsing of info fields on restart; use config specified
//! parameters instead
extern bool ignore_info;

//! Toggle interactions "on" or "off" by default.  If interactions are
//! "on" (the default), interactions listed in the 'Interaction' list
//! be turned "off".  Alternatively, if interactions are "off",
//! interactions listed in the 'Interaction' list be turned "on".
extern bool all_couples;

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

/**
   Quit now

   This flag is set by checkpoint itmer
*/
extern unsigned char quit_signal;

/** Implementation variables for handlers */
extern int stop_signal0;
extern int dump_signal0;

				// MPI variables

//! Used to flag separate behavior on inintialization (unused)
extern int is_init;

//! Internally defined worker communicator
extern MPI_Comm MPI_COMM_WORKER;

//! List of host names and ranks
extern std::map<std::string, std::vector<int> > nameMap;

//! List of sibling ranks
extern std::vector<int> siblingList;

/** @name Coefficient header on each output dump */
//@{
//! Cylindrical header type
extern struct CylCoefHeader  coefheadercyl;
//! Spherical header
extern struct SphCoefHeader  coefheadersph;
//@}


// Particle components

//! One component container per simulation
class ComponentContainer;
extern ComponentContainer *comp;

// External forces

//! One external force container per simulation
class ExternalCollection;
extern ExternalCollection *external;

// Output list

class OutputContainer;
//! One output routine container per simulation
extern OutputContainer *output;

//! Parameter database
extern YAML::Node parse;

//! Total alloted runtime (runtime < 0 turns off timer)
extern double runtime;

//! Last PS file name
extern string lastPS, lastPSQ, lastPSR;

//! Checkpoint timer
extern CheckpointTimer chktimer;

/** Restart command on successful completion but early termination
    according to the checkpoint timer */
extern string restart_cmd;

//! MPI datatype for pHOT key
extern MPI_Datatype MPI_EXP_KEYTYPE;

//! Barrier instance for debugging
extern BarrierWrapper *barrier;

//! Turn on Barrier check debugging
extern bool barrier_check;

//! Turn on Barrier check debugging
extern bool barrier_debug;

//! Turn on Barrier extra debugging
extern bool barrier_extra;

//! Barrier label checking
extern bool barrier_label;

//! Barrier label light weight method
extern bool barrier_light;

//! Do not report good barrier condtion if true
extern bool barrier_quiet;

//! Default memory limit value in GB.  If value is set to zero, use
//! the system default.  If value is less than zero, attempt to set
//! memory to unlimited.
extern int  rlimit_val;

//! Endless loop for setting gdb breakpoints
extern bool debug_wait;
extern bool main_wait;

//! Set to false to suppress cuda computation
extern bool use_cuda;

//! Can we leave phase space on GPUs or do we need to copy back to host?
extern bool leapfrog_cuda;

//! Mersenne random number generator provided by Boost.Random.  This
//! could be changed to another if needed.
extern unsigned int   random_seed;

//! exputils globals
#include <libvars.H>
using namespace __EXP__;

//! For toggling CUDA profiling
extern bool cuda_prof;

#include "Species.H"

#if HAVE_LIBCUDA==1
#include <cudaUtil.cuH>
extern thrust::device_vector<int> cuDstepL, cuDstepN;
extern int cudaGlobalDevice;
#endif

#endif
