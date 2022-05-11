// EXP global library variables for standalone routines
//
#ifndef _LIBVARS_H
#define _LIBVARS_H

#include <random>

extern char threading_on;
extern pthread_mutex_t mem_lock;
extern std::string outdir, runtag;
extern int nthrds;
extern int this_step;
extern unsigned multistep;
extern unsigned maxlev;
extern int mstep;
extern int Mstep;
extern int myid;
extern std::mt19937 random_gen;

#endif	// END _LIBVARS_H
