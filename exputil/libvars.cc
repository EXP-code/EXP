/**
   Globals used by EXP runtime environment supplied here for
   standalone utilities
 */

#include "libvars.H"

namespace __EXP__
{
  //@{
  //! POSIX thread support
  char             threading_on     = 0;
  pthread_mutex_t  mem_lock;
  //@}

  //! Location for output
  std::string      outdir           = "./";

  //! Run name for file labeling
  std::string      runtag           = "newrun";

  //! Number of POSIX threads (minimum: 1)
  int              nthrds           = 1;

  //! Multistep indices
  unsigned         multistep        = 0;

  //! Random number generator instance
  std::mt19937     random_gen;

  //! Source/line info exception strings
  bool             sourceline       = false;

  //! Sign convention element rank for eigenfunctions and eigenvectors
  int              nevsign          = 4;

  //! Sanity tolerance for orthogonality
  double           orthoTol         = 1.0e-2;
};


