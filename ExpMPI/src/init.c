/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine initializes working vectors
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *
 ***************************************************************************/

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include "expand.h"

void set_global_com(void);

#define BUF 256

void do_bodies_init(void)
{
  int i;

  /* Initialize vectors */

				/* Center of mass storage for all */
  com1 = dvector(0, 2);
  com2 = dvector(0, 2);
  for (i=0; i<3; i++) com1[i] = com2[i] = 0.0;

				/* Body vectors */

  mass = dvector(1,nbodmax);
  x = dvector(1,nbodmax);
  y = dvector(1,nbodmax);
  z = dvector(1,nbodmax);
  rr = dvector(1,nbodmax);
  vx = dvector(1,nbodmax);
  vy = dvector(1,nbodmax);
  vz = dvector(1,nbodmax);
  ax = dvector(1,nbodmax);
  ay = dvector(1,nbodmax);
  az = dvector(1,nbodmax);
  pot = dvector(1,nbodmax);
  potext = dvector(1,nbodmax);

  if (relx) esave = dvector(1,nbodmax);
  if (scatter) mfp = dvector(1,nbodmax);
  
  component = (int *)malloc((unsigned)nbodmax*sizeof(int));
  if (!component) {
    fprintf(stderr,"Couldn't allocate component vector, myid=%d\n", myid);
    exit(-1);
  }
  component--;
  
#ifdef SEQCHECK
  sequence  = (int *)malloc((unsigned)nbodmax*sizeof(unsigned int));
  if (!sequence) {
    fprintf(stderr,"Couldn't allocate sequence vector, myid=%d\n", myid);
    exit(-1);
  }
  sequence--;
#endif

  return;

}
