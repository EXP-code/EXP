#include <mpi.h>

static double MPL_accum_time=0.0;
static double MPL_last_time=0.0;

void MPL_reset_timer(void)
{
  MPL_accum_time=0.0;
}

void MPL_start_timer(void)
{
  MPL_last_time = MPI_Wtime();
}

void MPL_stop_timer(void)
{
  double curtime = MPI_Wtime();

  MPL_accum_time += curtime - MPL_last_time;
  MPL_last_time = curtime;
}

double MPL_read_timer(int reset)
{
  double save = MPL_accum_time;
  if (reset) MPL_accum_time = 0.0;
  return save;
}
