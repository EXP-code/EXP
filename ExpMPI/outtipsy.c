#define DEBUG 1

/*
  spit out phase-space, one particle per line
*/

#include <stdlib.h>

#include "expand.h"
#include "tipsydefs.h"

int t_ngas = 0;
int t_ndark = 0;
int t_nstar = 0;

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void out_tipsy(int n)
{
  int i;
  char file[80];
  FILE *fout;

  int ngas ;
  int ndark ;
  int nstar ;

  /*
    struct gas_particle *gp;
    */
  struct dark_particle *dp;
  struct star_particle *sp;

  sprintf(file, "%s.bin", outname);
  if ( (fout=fopen(file,"a+")) == NULL) {
    fprintf(stderr,"Couldn't open %s . . . quitting\n", file);
    return;
  }

				/* Count up types */
  ngas = 0;
  ndark = 0;
  nstar = 0;

  for (i=1; i<=nbodies; i++) {
    switch(component[i]) {
    case -1:			/* Slab particles */
      nstar++;
      break;
    case 0:			/* Point particles */
      ndark++;
      break;
    case 1:			/* Dark particles */
      ndark++;
      break;
    case 2:			/* Disk particles */
      nstar++;
      break;
    } /* No gas, yet . . . */
  }
	
#ifdef DEBUG
  fprintf(stderr, "T=%f    ndark, nstar, ngas: %d, %d, %d\n", tnow, ndark, nstar, ngas);
#endif

				/* Hard wire for time begin */
  header.time    = tnow;
  header.ndim    = 3;
  header.nbodies = nbodies;

  header.nsph    = ngas;
  header.ndark   = ndark;
  header.nstar   = nstar;


  /*======================================================================*/
  /* Allocate gas particles                                               */
  /*======================================================================*/

  if (ngas != t_ngas) {

    if(gas_particles != NULL) free(gas_particles);

    if(ngas != 0) {
      gas_particles = 
	(struct gas_particle *) malloc(ngas*sizeof(*gas_particles));

      if(gas_particles == NULL) {
	fprintf(stderr,	"<sorry, no memory for gas particles, master>\n") ;
	return ;
      }

      t_ngas = ngas;

    }
    else
      gas_particles = NULL;
  }


  /*======================================================================*/
  /* Allocate dark particles                                              */
  /*======================================================================*/


  if (ndark != t_ndark) {

    if(dark_particles != NULL) free(dark_particles);

    if(ndark != 0) {
      dark_particles = 
	(struct dark_particle *) malloc(ndark*sizeof(*dark_particles));

      if(dark_particles == NULL) {
	fprintf(stderr,	"<sorry, no memory for dark particles, master>\n") ;
	return ;
      }

      t_ndark = ndark;

    }
    else
      dark_particles = NULL;
  }


  /*======================================================================*/
  /* Allocate star particles                                              */
  /*======================================================================*/
  
  if (nstar != t_nstar) {

    if(star_particles != NULL) free(star_particles);

    if(nstar != 0) {
      star_particles =
	(struct star_particle *)malloc(nstar*sizeof(*star_particles));

      if(star_particles == NULL) {
	fprintf(stderr,
		"<sorry, no memory for star particles, master>\n") ;
	return ;
      }

      t_nstar = nstar;

    }
    else
      star_particles = NULL;
  }

	
  /*======================================================================*/
  /* Populate tipsy structures                                            */
  /*======================================================================*/

  dp = dark_particles;
  sp = star_particles;

  for (i=1; i<=nbodies; i++) {

    switch(component[i]) {
    case 0:			/* Point */
    case 1:			/* Halo */
      
      dp->mass = mass[i];

      dp->pos[0] = x[i];
      dp->pos[1] = y[i];
      dp->pos[2] = z[i];

      dp->vel[0] = vx[i];
      dp->vel[1] = vy[i];
      dp->vel[2] = vz[i];

      dp->eps = 0;

      dp->phi = pot[i];

      dp++;

      break;
    case -1:			/* Slab */
    case 2:			/* Disk */
	
      sp->mass = mass[i];

      sp->pos[0] = x[i];
      sp->pos[1] = y[i];
      sp->pos[2] = z[i];

      sp->vel[0] = vx[i];
      sp->vel[1] = vy[i];
      sp->vel[2] = vz[i];

      sp->metals = 0;
      sp->tform = 0;
      sp->eps = 0;

      sp->phi = pot[i];

      sp++;

      break;
    }

  }

  /*======================================================================*/
  /* Write TIPSY binary file                                              */
  /*======================================================================*/

  fwrite((char *)&header, sizeof(header), 1, fout) ;
  fwrite((char *)gas_particles, sizeof(struct gas_particle),
	 ngas, fout) ;
  fwrite((char *)dark_particles, sizeof(struct dark_particle),
	 ndark, fout) ;
  fwrite((char *)star_particles, sizeof(struct star_particle),
	 nstar, fout) ;

  fclose(fout);
}


void out_chkpt(void)
{
  int i;
  char file[80];
  FILE *fout;

  int ngas ;
  int ndark ;
  int nstar ;

  /*
    struct gas_particle *gp;
    */
  struct dark_particle *dp;
  struct star_particle *sp;

  sprintf(file, "%s.chkpt", outname);
  if ( (fout=fopen(file,"w")) == NULL) {
    fprintf(stderr,"Couldn't open checkpoint file <%s> . . . going on\n", file);
    return;
  }

				/* Count up types */
  ngas = 0;
  ndark = 0;
  nstar = 0;

  for (i=1; i<=nbodies; i++) {
    switch(component[i]) {
    case -1:			/* Slab particles */
      nstar++;
      break;
    case 0:			/* Point particles */
    case 1:			/* Disk particles */
      ndark++;
      break;
    case 2:			/* Dark particles */
      nstar++;
      break;
    } /* No gas, yet . . . */
  }
	
#ifdef DEBUG
  fprintf(stderr, "T=%f    ndark, nstar, ngas: %d, %d, %d\n", tnow, ndark, nstar, ngas);
#endif

				/* Hard wire for time begin */
  header.time    = tnow;
  header.ndim    = 3;
  header.nbodies = nbodies;

  header.nsph    = ngas;
  header.ndark   = ndark;
  header.nstar   = nstar;


  /*======================================================================*/
  /* Allocate gas particles                                               */
  /*======================================================================*/

  if (ngas != t_ngas) {

    if(gas_particles != NULL) free(gas_particles);

    if(ngas != 0) {
      gas_particles = 
	(struct gas_particle *) malloc(ngas*sizeof(*gas_particles));

      if(gas_particles == NULL) {
	fprintf(stderr,	"<sorry, no memory for gas particles, master>\n") ;
	return ;
      }

      t_ngas = ngas;

    }
    else
      gas_particles = NULL;
  }


  /*======================================================================*/
  /* Allocate dark particles                                              */
  /*======================================================================*/


  if (ndark != t_ndark) {

    if(dark_particles != NULL) free(dark_particles);

    if(ndark != 0) {
      dark_particles = 
	(struct dark_particle *) malloc(ndark*sizeof(*dark_particles));

      if(dark_particles == NULL) {
	fprintf(stderr,	"<sorry, no memory for dark particles, master>\n") ;
	return ;
      }

      t_ndark = ndark;

    }
    else
      dark_particles = NULL;
  }


  /*======================================================================*/
  /* Allocate star particles                                              */
  /*======================================================================*/
  
  if (nstar != t_nstar) {

    if(star_particles != NULL) free(star_particles);

    if(nstar != 0) {
      star_particles =
	(struct star_particle *)malloc(nstar*sizeof(*star_particles));

      if(star_particles == NULL) {
	fprintf(stderr,
		"<sorry, no memory for star particles, master>\n") ;
	return ;
      }

      t_nstar = nstar;

    }
    else
      star_particles = NULL;
  }

	
  /*======================================================================*/
  /* Populate tipsy structures                                            */
  /*======================================================================*/

  dp = dark_particles;
  sp = star_particles;

  for (i=1; i<=nbodies; i++) {

    switch(component[i]) {
    case 0:			/* Point */
    case 1:			/* Halo */
      
      dp->mass = mass[i];

      dp->pos[0] = x[i];
      dp->pos[1] = y[i];
      dp->pos[2] = z[i];

      dp->vel[0] = vx[i];
      dp->vel[1] = vy[i];
      dp->vel[2] = vz[i];

      dp->eps = 0;

      dp->phi = pot[i];

      dp++;

      break;
    case -1:			/* Slab */
    case 2:			/* Disk */
	
      sp->mass = mass[i];

      sp->pos[0] = x[i];
      sp->pos[1] = y[i];
      sp->pos[2] = z[i];

      sp->vel[0] = vx[i];
      sp->vel[1] = vy[i];
      sp->vel[2] = vz[i];

      sp->metals = 0;
      sp->tform = 0;
      sp->eps = 0;

      sp->phi = pot[i];

      sp++;

      break;
    }

  }

  /*======================================================================*/
  /* Write TIPSY binary file                                              */
  /*======================================================================*/

  fwrite((char *)&header, sizeof(header), 1, fout) ;
  fwrite((char *)gas_particles, sizeof(struct gas_particle),
	 ngas, fout) ;
  fwrite((char *)dark_particles, sizeof(struct dark_particle),
	 ndark, fout) ;
  fwrite((char *)star_particles, sizeof(struct star_particle),
	 nstar, fout) ;

  fclose(fout);
}
