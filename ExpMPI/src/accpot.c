/*
  Compute accelerations, potential, and density.
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#undef TESTACC

extern void get_acceleration_and_potential_bes(void);
extern void get_acceleration_and_potential_CB(void);
extern void get_acceleration_and_potential_HERNQ(void);
extern void get_acceleration_and_potential_SLsph(void);
extern void get_acceleration_and_potential_CBDisk(void);
extern void get_acceleration_and_potential_cube(void);
extern void get_acceleration_and_potential_slab(void);
extern void get_acceleration_and_potential_slabSL(void);
extern void get_acceleration_and_potential_n2(void);
extern void get_acceleration_and_potential_Cyl(void);
extern void external_halo(void);
extern void external_shock(void);
extern void scatter_mfp(void);
extern void clouds(void);
extern void user_perturbation(void);

void compute_potential(void)
{
  int i;

#ifdef TESTACC
  FILE *out = NULL;
  
  if (component[1] == 0 && myid) {
    out = fopen("testacc.dat", "a");
    fprintf(out, "%3d %13.6e", myid, tnow);
  }
#endif

				/* Zero-out external potential */
  for (i=1; i<=nbodies; i++)
    potext[i] = 0.0;
  
				/* Zero-out potential and acceleration */

  for (i=1; i<=nbodies; i++) {
    ax[i] = 0.0;
    ay[i] = 0.0;
    az[i] = 0.0;
    pot[i] = 0.0;
  }

				/* Compute new accelerations and potential */
  if (selfgrav) {

    int done=0;
				/* Count particles contributing */
    used = 0;

				/* Component #0: point masses */
    if (pmnum0) {
      get_acceleration_and_potential_pointmass();
      done = 1;
    }

#ifdef TESTACC
    if (component[1] == 0 && myid) {
      fprintf(out, "  %13.6e  %13.6e  %13.6e", ax[1], ay[1], az[1]);
    }
#endif
				/* Component #1: halo/spherical particles */
    if (bessel_sph) { 
      get_acceleration_and_potential_bes();
      done = 1;
    }
    else if (c_brock) {
      get_acceleration_and_potential_CB();
      done = 1;
    }
    else if (c_brock_disk) {
      get_acceleration_and_potential_CBDisk();
      done = 1;
    }
    else if (hernq) {
      get_acceleration_and_potential_HERNQ();
      done = 1;
    }
    else if (sphereSL) {
      get_acceleration_and_potential_SLsph();
      done = 1;
    }
    else if (cube) {
      get_acceleration_and_potential_cube();
      done = 1;
    }
    else if (slab) {
      if (slabSL)
	get_acceleration_and_potential_slabSL();
      else
	get_acceleration_and_potential_slab();
      done = 1;
    }

#ifdef TESTACC
    if (component[1] == 0 && myid) {
      fprintf(out, "  %13.6e  %13.6e  %13.6e", ax[1], ay[1], az[1]);
    }
#endif

				/* Component #2: disk particles */
    if (cylinder) {
      get_acceleration_and_potential_Cyl();
      done = 1;
    }
    
#ifdef TESTACC
    if (component[1] == 0 && myid) {
      fprintf(out, "  %13.6e  %13.6e  %13.6e\n", ax[1], ay[1], az[1]);
    }
#endif

				/* No self-gravity whatsoever! */
    if (nulltest) {
      done = 1;
    }

    if (!done) {
      fprintf(stderr,"No potential selected!\n");
      exit(-1);
    }

  }

  if (myid == 0) return;



	/* include tidal perturbation (KL 5/27/92) */

  if (tides)
    {
      tidal_field();
    }
  

  if (shock)			/* (MDW 3/26/94) */
    {
      external_shock();
    }

				/* One would put in code to compute
				   an external potential here */

  if (halo) 			/* (MDW 3/10/95) */
    {
      external_halo();
    }


  if (scatter) 			/* (MDW 10/16/99) */
    {
      scatter_mfp();
    }

				/* To be supplied by user for testing */

  if (user)
    {
      user_perturbation();
    }


#ifdef TESTACC
  if (out) fclose(out);
#endif
}

