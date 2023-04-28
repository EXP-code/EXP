/**
   Description:
   -----------
 
   This routine is an NBODY code whose force calculation is based 
   on a biorthogonal phase-space expansion; for now, the expansion
   set is based on the eigenfunctions of the spherical Laplacian:
 
              Y_{lm} * j_n
                ^       ^
                |       |
                |       spherical Bessel functions
                | spherical harmonics
 
 
   Call sequence:
   -------------
   exp [-f parameter file] [ [keyword=value] [keyword=value] ... ]
 
 
   Parameters:
   ----------
   The parameter file may be used instead of the keyword/value pairs.
   The available keywords are: 
 
   KEYWORD        EXPLANATION                            DEFAULT
   ----------     ----------------------------           --------
   lmax =         maximum harmonic order                 4
   nmax =         maxinum radial order                   10
   nlog =         interval between log output            10
   nlist =        interval between p-s output            50
   nsteps =       total number of steps to run           500
   fixacc =       add off set to acceleration to
                  conserve momentum                      1  (1=yes 0=no)
   selfgrav =     turn-on self-gravity                   1
   outbods =      output phase space                     1
   outengy =      output per particle energies           0
   diag =         diagnostic output                      0
   dtime =        time step                              0.1
   rmax =         maximum radial extent                  (computed from input)
   logfile =      logfile name                           LOG.FILE
   infile =       input parameter file                   IN.FILE
   outname =      output file prefix                     OUT
   parmfile =     parameter file name                    PARMS.FILE



  $Id$
*/

#ifndef _expand_H
#define _expand_H

#include <config_exp.h>

#define BSD_SOURCE

#include <mpi.h>

#include <pthread.h>

#include <stdexcept>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <limits>

#include <Eigen/Eigen>

using namespace std;

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

				// Turn on sequence checking on debug
#ifdef DEBUG
#define SEQCHECK
#endif

				// Function declarations
void init_velocity(void);
void begin_run(void);
void incr_position(double dt, int mlevel=0);
void incr_velocity(double dt, int mlevel=0);
void incr_com_position(double dt);
void incr_com_velocity(double dt);
void write_parm(void);
void initialize_multistep();
void sync_eval_multistep();
void adjust_multistep_level(bool all=true);
void initialize_cuda(void);



#ifndef MAX
#define MAX(A,B) (A>B ? A : B)
#endif
#ifndef MIN
#define MIN(A,B) (A<B ? A : B)
#endif

//*****Stuff for expand.c*****

double get_dens(double r, int l, double *coef);
double densi(double rr, int l, int n),potli(double rr, int l, int n);
void make_grid(double rmin, double rmax, int lmax, int nmax);
void get_potacc(double r, int l, double *coef, double *p, double *dp);
void set_global_com(void);
void create_tipsy(void);

//******Mathematical utilities******
double zbrent();
extern "C" double dgammln();
void locate();
double plgndr(int, int, double);
double dplgndr(int l, int m, double x);
double factrl(int n);

// Constants

static const double DSMALL = 1.0e-16;
static const int    mdtDim = 7;

// Parallel prototypes

void setup_distribution(void);
void redistribute_particles(void);
void test_mpi(void);
void parallel_gather_coefficients(void);
void parallel_distribute_coefficients(void);
void recompute_processor_rates(void);
void MPL_reset_timer(void);
void MPL_start_timer(void);
void MPL_stop_timer(void);
double MPL_read_timer(int reset);

// Utilities

void get_ultra(int nmax, double l, double x, Eigen::VectorXd& p);
string trimRight(const string);
string trimLeft(const string);
string trimComment(const string);
inline double atof(string& s) {return atof(s.c_str());}
inline int atoi(string& s) {return atoi(s.c_str());}
inline bool atol(string& s) {return atoi(s.c_str()) ? true : false;}


// Component structures

#include <ComponentContainer.H>

// External forces

#include <ExternalForce.H>

// Output list

#include <Output.H>

void do_output_init(void);

// Parameter database

#include <yaml-cpp/yaml.h>
#include <ComponentContainer.H>
#include <global.H>

// Exception handler

#include <EXPException.H>

#endif

