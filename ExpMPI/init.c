/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine reads in phase-space and initializes working vectors
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

void read_bodies_and_init(void)
{
  int i,iret;
  FILE *fin;
  double max=0.0,r2,tmp;
  char buf[BUF];

  /* Initialize vectors */

				/* Center of mass storage for all */
  com1 = dvector(0, 2);
  com2 = dvector(0, 2);
  for (i=0; i<3; i++) com1[i] = com2[i] = 0.0;

				/* Body vectors */
  if (myid>0) {

    mass = dvector(1,nbodmax);
    initial_mass = dvector(1,nbodmax);
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
  
    component = (int *)malloc((unsigned)nbodmax*sizeof(int));
    if (!component) {
      fprintf(stderr,"Couldn't allocate component vector, myid=%d\n", myid);
      exit(-1);
    }
    component--;
    
    return;
				/* All slaves return */  
  }


  /* Read in phase space */

  if ( (fin=fopen(infile,"r")) == NULL) {
    fprintf(stderr,"Couldn't open %s . . . quitting\n",infile);
    exit(-1);
  }

  fscanf(fin,"%d %lf\n",&nbodies,&tnow);
  if (nbodies > nbodmax*slaves) {
    fprintf(stderr,"Not enough space on all processors to hold phase space\nnbodmax should be at least %d\n", (int)( (double)nbodies/slaves + 1) );
    exit(-1);
  }

  mass = dvector(1,nbodies);
  initial_mass = dvector(1,nbodies);
  x = dvector(1,nbodies);
  y = dvector(1,nbodies);
  z = dvector(1,nbodies);
  rr = dvector(1,nbodies);
  vx = dvector(1,nbodies);
  vy = dvector(1,nbodies);
  vz = dvector(1,nbodies);
  ax = dvector(1,nbodies);
  ay = dvector(1,nbodies);
  az = dvector(1,nbodies);
  pot = dvector(1,nbodies);
  potext = dvector(1,nbodies);

  component = (int *)malloc((unsigned)nbodies*sizeof(int));
  if (!component) {
    fprintf(stderr,"Couldn't allocate component vector, myid=%d\n", myid);
    exit(-1);
  }
  component--;

  if (fgets(buf,BUF,fin) == NULL) {
    fprintf(stderr,"Trouble reading file %s . . . quitting\n",infile);
    exit(-1);
  }

  iret = sscanf(buf,"%d %lf %lf %lf %lf %lf %lf %lf %lf",
		&component[1], &mass[1],
		&x[1],&y[1],&z[1],&vx[1],&vy[1],&vz[1],&tmp);

  if (iret==8) restart = 0;
  else if (iret==9) restart = 1;
  else {
    fprintf(stderr,"Trouble parsing string from %s . . . quitting\n",infile);
    exit(-1);
  }

  if (restart)
    printf("Restart . . .\n");
  else
    printf("New run . . .\n");
  fflush(stdout);


  pot[1] = potext[1] = 0.0;
  for (i=2; i<=nbodies; i++) {
    if (!restart) {
      fscanf(fin,"%d %lf %lf %lf %lf %lf %lf %lf",
	     &component[i],
	     &mass[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i]);
    }
    else {
      fscanf(fin,"%d %lf %lf %lf %lf %lf %lf %lf %lf",
	     &component[i],
	     &mass[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i],&tmp);
    }
    pot[i] = potext[i] = 0.0;
    r2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
    max = MAX(r2,max);
  }

				/* Default: set to max radius */

  if (rmax <= 0.0) rmax = sqrt(max);

				/* Initialize time */
  tpos = tvel = tnow;

				/* store initial masses */
				/* 
				    needed for mass_loss(), KL - 7/3/92
				*/
  for (i=1; i<=nbodies; i++)
  {
	initial_mass[i] = mass[i];
  }

				/* Set center of mass and com velocity
				   to zero for each component */
  if (zerocom || zerovel) set_global_com();

}
