#include <mpi.h>
#include <string>

#include <coef.H>
#include <ComponentContainer.H>
#include <ExternalForce.H>
#include <ExternalCollection.H>
#include <OutputContainer.H>
#include <BarrierWrapper.H>
#include <chkTimer.H>

				// Numerical parameters

int nsteps = 500;		// Number of steps to execute
int nscale = 20;		// Number of steps between rescaling
int ngpus  = 0;	                // Number of GPUs per node (0 means use all available)
int nbalance = 0;		// Steps between load balancing
int nreport = 0;		// Steps between particle reporting
double dbthresh = 0.05;		// Load balancing threshold (5% by default)
double dtime = 0.1;		// Default time step size
unsigned nbits = 32;		// Number of bits per dimension
unsigned pkbits = 6;		// Number of bits for parition
unsigned PFbufsz = 40000;	// ParticleFerry buffer size in particles


bool restart = false;		// Restart from a checkpoint
bool use_cwd = false;		// Use Node 0's current working directory on all nodes
int NICE = 0;			// Niceness level (default: 0)
int VERBOSE = 1;		// Chattiness for standard output
bool step_timing = false;	// Time parts of the step (set true by
				// DEFAULT>3)
bool initializing = false;	// Used by force methods to do "private things"
				// before the first step (e.g. run through
				// coefficient evaluations even when
				// self_consistent=false)

				// Total alloted runtime: 0 means ignore
double runtime = 0.0;

				// Files
string homedir  = "./";
string infile   = "restart.in";
string parmfile = "config";
string ratefile = "processor.rates";
string ldlibdir = ".";

double tnow     = 0.0;		// Per step variables
int psdump      = -1;
int this_step   = 0;
				// Global center of mass
double mtot;
double *gcom = new double [3];
double *gcov = new double [3];
bool global_cov = false;
bool eqmotion = true;
unsigned char stop_signal  = 0;
unsigned char dump_signal  = 0;
unsigned char quit_signal  = 0;
				// Multistep variables
unsigned        shiftlevl  = 0;
unsigned        maxlev     = 100;
int             mdrft      = 0;
int             mstep      = 0;
int             Mstep      = 1;
int             centerlevl = -1;
				// Timestep control
bool            DTold      = false;
double          dynfracS   = 1.00;
double          dynfracD   = 1.0e32;
double          dynfracV   = 0.01;
double          dynfracA   = 0.03;
double          dynfracP   = 0.05;

std::vector<int> mfirst, mintvl, stepL, stepN;
std::vector< vector<bool> > mactive;
std::vector< vector<int> > dstepL, dstepN;

#if HAVE_LIBCUDA==1
int cudaGlobalDevice;
#endif


				// Multithreading data structures for
				// incr_position and incr_velocity
/*struct thrd_pass_posvel 
{
  double dt;
  int mlevel;
  int id;
};*/

vector<thrd_pass_posvel> posvel_data;
vector<pthread_t> posvel_thrd;
int is_init=1;

				// List of host names and ranks
std::map<std::string, std::vector<int> > nameMap;

				// List of sibling ranks
std::vector<int> siblingList;


CylCoefHeader  coefheadercyl;
SphCoefHeader  coefheadersph;

ComponentContainer *comp = 0;
ExternalCollection *external = 0;
OutputContainer    *output = 0;
YAML::Node          parse;

map<string, maker_t *, less<string> > factory;
map<string, maker_t *, less<string> >::iterator fitr;

std::string lastPS, lastPSQ, lastPSR;
CheckpointTimer chktimer;
string restart_cmd;

unsigned int random_seed = 11;

MPI_Datatype MPI_EXP_KEYTYPE;

BarrierWrapper *barrier = 0;
bool barrier_check = false;
bool barrier_debug = false;
bool barrier_extra = false;
bool barrier_label = true;
bool barrier_light = true;
bool barrier_quiet = true;

bool cuda_prof     = false;
bool debug_wait    = false;
bool main_wait     = false;
bool mpi_wait      = false;
bool fpe_trap      = false;
bool fpe_trace     = false;
bool fpe_wait      = false;
bool traceback     = false;

bool ignore_info   = false;
bool all_couples   = true;

int  rlimit_val    = 0;
bool use_cuda      = false;
