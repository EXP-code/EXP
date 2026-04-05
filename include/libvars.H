// EXP global library variables for standalone routines
//
#ifndef _LIBVARS_H
#define _LIBVARS_H

#include <random>

#include <mpi.h>

#include "config_exp.h"		// For version string
#include "EXPversion.H"		// For compile-time version parsing

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

  //! Sign convention element rank for eigenfunctions and eigenvectors
  extern int nevsign;

  //! Sanity tolerance for orthogonality
  extern double orthoTol;

  //! Compile-time version parsing
  static constexpr EXPversion exp_build = EXP_parse_version(VERSION);

};

#endif	// END _LIBVARS_H
