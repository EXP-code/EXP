/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine is an NBODY code whose force calculation is based 
 *  on a biorthogonal phase-space expansion; for now, the expansion
 *  set is based on the eigenfunctions of the spherical Laplacian:
 *
 *             Y_{lm} * j_n
 *               ^       ^
 *               |       |
 *               |       spherical Bessel functions
 *               | spherical harmonics
 *
 *
 *  Call sequence:
 *  -------------
 *  exp [-f parameter file] [ [keyword=value] [keyword=value] ... ]
 *
 *
 *  Parameters:
 *  ----------
 *  The parameter file may be used instead of the keyword/value pairs.
 *  The available keywords are: 
 *
 *  KEYWORD        EXPLANATION                            DEFAULT
 *  ----------     ----------------------------           --------
 *  nbodmax =      maximum # of bodies                    20000
 *  lmax =         maximum harmonic order                 4
 *  nmax =         maxinum radial order                   10
 *  nlog =         interval between log output            10
 *  nlist =        interval between p-s output            50
 *  nsteps =       total number of steps to run           500
 *  fixacc =       add off set to acceleration to
 *                 conserve momentum                      1  (1=yes 0=no)
 *  selfgrav =     turn-on self-gravity                   1
 *  outbods =      output phase space                     1
 *  outengy =      output per particle energies           0
 *  diag =         diagnostic output                      0
 *  dtime =        time step                              0.1
 *  rmax =         maximum radial extent                  (computed from input)
 *  logfile =      logfile name                           LOG.FILE
 *  infile =       input p-s file                         IN.FILE
 *  outname =      output file prefix                     OUT
 *  parmfile =     parameter file name                    PARMS.FILE
 *
 *
 *  Compile options:
 *  ---------------
 *
 *  FPETRAP        turn on floating point trapping
 *                 (for debugger only)
 *
 *  Notes:
 *  -----
 *
 *  Modifications:
 *  -------------
 *
 *	KL 3/20/92 - changed to allow fixed potential calculation.
 *		Expansion coefficients are computed during initialisation
 *		and then left unchanged. Set flag selfcons=0 to use
 *		fixed potential.
 *
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *
 ***************************************************************************/
#define MAIN 1
#include "expand.h"
#include "version.h"

void do_step(int);
void clean_up(void);

#include <sys/time.h>
#include <sys/resource.h>

/*
static char rcsid[] = "$Id$";
*/

void MPL_parse_args(int argc, char** argv);

int main(int argc, char** argv)
{

  int *nslaves, n;
  MPI_Group world_group, slave_group;

  /*===================*/
  /* MPI preliminaries */
  /*===================*/

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

				/* Make SLAVE group */
  slaves = numprocs - 1;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  nslaves = (int *)malloc(slaves*sizeof(int));
  if (!nslaves) {
    fprintf(stderr, "main: problem allocating <nslaves>\n");
    exit(-1);
  }
  for (n=1; n<numprocs; n++) nslaves[n-1] = n;
  MPI_Group_incl(world_group, slaves, nslaves, &slave_group);
  MPI_Comm_create(MPI_COMM_WORLD, slave_group, &MPI_COMM_SLAVE);
  free(nslaves);

    
  /* Debug id */
  MPI_Group_rank ( slave_group, &n );
  fprintf(stderr, "Process %d on %s    rank in SLAVE: %d\n", 
	   myid, processor_name, n);

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef MPE_PROFILE
  MPE_Init_log();

  if (myid == 0) {
    MPE_Describe_state(1, 2, "Distribute particles", "red:dimple3" );
    MPE_Describe_state(3, 4, "Gather particles", "green:dllines3" );
    MPE_Describe_state(5, 6, "Gather coefs", "cyan:hlines2" );
    MPE_Describe_state(7, 8, "Distribute coefs", "yellow:drlines4" );
    MPE_Describe_state(9, 10, "Compute coefs", "magenta:vlines3" );
    MPE_Describe_state(11, 12, "Compute forces", "orange3:gray" );
    MPE_Describe_state(13, 14, "Advance time", "purple:boxes" );
    MPE_Describe_state(15, 16, "Send energies", "blue:dllines4" );
  }
#endif


#ifdef DEBUG
  //  sleep(20);
#endif

  /*================*/
  /* Print welcome  */
  /*================*/

  if (myid==0) {
    fprintf(stdout, "\nThis is Exp %s\n\n", version_id);
  }


  /*============================*/
  /* Parse command line:        */
  /* broadcast to all processes */
  /*============================*/

  MPL_parse_args(argc,argv);


  /*================*/
  /* Nice process ? */
  /*================*/

  if (NICE>0) setpriority(PRIO_PROCESS, 0, NICE);


  /*==============================================*/
  /* Read in points and initialize expansion grid */
  /*==============================================*/

  begin_run();


  /*===========*/
  /* MAIN LOOP */
  /*===========*/

  for (this_step=1; this_step<=nsteps; this_step++) do_step(this_step);


  /*===========*/
  /* Finish up */
  /*===========*/

  clean_up();

  return 0;
}


