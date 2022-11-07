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

  //! Multistep levels (default: 0 means no multistepping)
  extern unsigned multistep;

  //! Random generator instance
  extern std::mt19937 random_gen;

  //! Control flag for source/line debug info
  extern bool sourceline;

};

#endif	// END _LIBVARS_H
