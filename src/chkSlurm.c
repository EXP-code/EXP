/**
   Wrap slurm lib call to protect incompatible C structures from C++
 */

#include "config_exp.h"
#ifdef HAVE_LIBSLURM
#include <printf.h>
#include <slurm/slurm.h>

long rem_time(int jobid)
{
  return slurm_get_rem_time(jobid);
}

#else

long rem_time(int jobid)
{
  return 0;
}

#endif

