#include <mpi.h>
#include <string>

#include <coef.H>
#include <ComponentContainer.H>
#include <ExternalForce.H>
#include <ExternalCollection.H>
#include <OutputContainer.H>
#include <ParamParseMPI.H>

				// Numerical parameters

int nbodmax = 20000;		// Maximum number of bodies; this is not
				// an intrinsic limitation just a sanity
				// value

int nsteps = 500;		// Number of steps to execute
int nscale = 20;		// Number of steps between rescaling
int nthrds = 2;			// Number of POSIX threads
int nbalance = 0;		// Steps between load balancing
double dbthresh = 0.05;		// Load balancing threshold (5% by default)
double dtime = 0.1;		// Default time step size

bool restart = false;		// Restart from a checkpoint
bool use_cwd = false;		// Use Node 0's current working directory on all nodes
int NICE = 10;			// Default niceness level
int VERBOSE = 1;		// Chattiness for standard output

				// Files
string homedir = "./";
string infile = "restart.in";
string parmfile = "PARAM";
string ratefile = "processor.rates";
string outdir = "";
string runtag = "newrun";
string ldlibdir = ".";

double tnow;			// Per step variables
int this_step;

				// Global center of mass
double mtot;
double *gcom = new double [3];
double *gcov = new double [3];
bool global_cov = false;
bool eqmotion = true;
unsigned char stop_signal = 0;
unsigned char dump_signal = 0;
				// Multistep variables
unsigned multistep = 0;
bool posnsync = true;
double dynfracV = 0.01;
double dynfracA = 0.03;
int Mstep = 0;
int mstep = 0;
vector<int> mfirst, mintvl, levpop, stepL, stepN;
vector< vector<bool> > mactive;


				// Multithreading data structures for
				// incr_position and incr_velocity
struct thrd_pass_posvel 
{
  double dt;
  int mlevel;
  int id;
};

vector<thrd_pass_posvel> posvel_data;
vector<pthread_t> posvel_thrd;


				// MPI variables
int is_init=1;
int numprocs, slaves, myid, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

MPI_Comm MPI_COMM_SLAVE;

char threading_on = 0;
pthread_mutex_t mem_lock;

CoefHeader coefheader;
CoefHeader2 coefheader2;

ComponentContainer comp;
ExternalCollection external;
OutputContainer output;
ParamParseMPI *parse;

map<string, maker_t *, less<string> > factory;
map<string, maker_t *, less<string> >::iterator fitr;

