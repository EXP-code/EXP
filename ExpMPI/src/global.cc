#include <mpi.h>
#include <string>

#include <coef.H>
#include <ComponentContainer.H>
#include <ExternalForce.H>
#include <ExternalCollection.H>
#include <OutputContainer.H>
#include <ParamParseMPI.H>

				// Numerical parameters (nearly unchanged)
int nbodmax = 20000;
int nsteps = 500;
int nscale = 20;
int nthrds = 2;			// Number of POSIX threads
double dtime = 0.1;
double rmax_tidal = 1.0e+04;

bool use_cwd = true;
bool restart = false;
int NICE = 10;

				// Files
string homedir = "./";
string infile = "restart.in";
string parmfile = "PARMS.FILE";
string ldlibdir = ".";

double tpos, tvel, tnow;	// Per step variables
int this_step;

				// Global center of mass
double mtot;
double *gcom = new double [3];
double *gcov = new double [3];
bool global_cov = false;
bool fixacc = false;

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

