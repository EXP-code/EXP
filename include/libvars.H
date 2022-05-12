// EXP global library variables for standalone routines
//
#ifndef _LIBVARS_H
#define _LIBVARS_H

#include <random>

#include <mpi.h>

namespace __EXP__
{
  //! @name Theading variables
  //@{
  //! signal that threading is in use
  extern char threading_on;

  //! mutex for memory
  extern pthread_mutex_t mem_lock;
  //@}

  //! Location for output directory
  extern std::string outdir;

  //! Run tag
  extern std::string runtag;

  //! Number of POSIX threads per process (e.g. one per processor)
  extern int nthrds;

  //! Current Nbody step id
  extern int this_step;

  //! Multistep levels (default: 0 means no multistepping)
  extern unsigned multistep;

  //! Multistep sublevel info used for nbody
  extern unsigned maxlev;

  //! Current point in multistepping for velocity step
  extern int mstep;

  //! Maximum number of multistep steps
  extern int Mstep;

  //! Random generator instance
  extern std::mt19937 random_gen;
};

#endif	// END _LIBVARS_H
