// This may look like C code, but it is really -*- C++ -*-

// #define DEBUG 1

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  the Cylindrical biorthogonal expansion
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
 *  Value
 *
 *  Notes:
 *  -----
 *
 *	KL 5/27/92   Modified to allow freezing of particles beyond some cutoff
 *    radius. This is needed when running systems in tidal fields. 
 *    The function freeze_particle() is called for each particle
 *    to decide whether to freeze. *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *      06/09/92 updated to use recursion relations rather than tables
 *
 ***************************************************************************/

#include <values.h>
#include <WghtOrth3.h>
#include "expand.h"

static char rcsid[] = "$Id$";

void determine_coefficients_Cyl(void);

extern "C" {
  void determine_acceleration_and_potential_Cyl(void);
  void determine_fields_at_point_Cyl(double r, double z, double phi,
				     double *tdens, double *tpotl, 
				     double *tpotr, double *tpotz, 
				     double *tpotp);
}

//======================================================================
//======================================================================
//======================================================================


static CylindricalSL *ortho;

extern "C" 
void get_acceleration_and_potential_Cyl(void)
{
  static int firstime=1;
  int l, m, n;

  if (firstime) {
    CylindricalSL::ZMAX = zmax;
    SLGridCyl::A = acyl;
    ortho = new CylindricalSL(nmax2, nfft, mmax2);
  }

  /*======================*/
  /* Compute coefficients */
  /*======================*/

  if (firstime || self_consistent) {
    firstime = 0;
    determine_coefficients_Cyl();
  }

  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  determine_acceleration_and_potential_Cyl();

  MPL_stop_timer();

}

void determine_coefficients_Cyl(void)
{
  int i, l;
  int use0, use1;
  double r, r2, phi;
  double xx, yy, zz;

#ifdef MPE_PROFILE
  MPE_Log_event(9, myid, "b_compute_coef");
#endif

  use0 = 0;
  use1 = 0;

  ortho->setup_accumulation();
    
  /*		Begin by finding positions */

  for (i=1; i<=nbodies; i++) {

    if (freeze_particle(i)) continue;		/* frozen particles don't
						   contribute to field.
						   KL 5/27/92 */
    xx = x[i] - com2[0];
    yy = y[i] - com2[1];
    zz = z[i] - com2[2];

    r2 = (xx*xx + yy*yy);
    r = rr[i] = sqrt(r2) + DSMALL;

    if (component[i] != 2) continue;

    if (r<=rmax && fabs(zz)<=zmax) {
      use1++;
      phi = atan2(yy,xx);
      
      ortho->accumulate(r, zz, phi, mass[i]);
    }
  }
  
  ortho->make_coefficients();

  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (myid==0) used += use0;

#ifdef MPE_PROFILE
  MPE_Log_event(10, myid, "e_compute_coef");
#endif

}

void check_force_values(double phi, double p, double fr, double fz, double fp)
{
  if (
      isinf(phi) || isnan(phi) ||
      isinf(p  ) || isnan(p  ) ||
      isinf(fr ) || isnan(fr ) ||
      isinf(fz ) || isnan(fz ) ||
      isinf(fp ) || isnan(fp ) ) 
    {
      cerr << "check_force_values: Illegal value\n";
    }
}


void determine_acceleration_and_potential_Cyl(void)
{
  int i, l;
  double r, r2, phi;
  double xx, yy, zz;
  double p, d, fr, fz, fp;

  
#ifdef MPE_PROFILE
  MPE_Log_event(11, myid, "b_compute_force");
#endif

  /* Determine potential and acceleration */

  for (i=1; i<=nbodies; i++) {

    if (!disk_on_halo && component[i] != 2) continue;

    if (freeze_particle(i)) continue;

    xx = x[i] - com2[0];
    yy = y[i] - com2[1];
    zz = z[i] - com2[2];

    r2 = xx*xx + yy*yy;
    r = rr[i] = sqrt(r2) + DSMALL;
    phi = atan2(yy, xx);

    ortho->accumulated_eval(r, zz, phi, p, fr, fz, fp);
    
#ifdef DEBUG
    check_force_values(phi, p, fr, fz, fp);
#endif

    pot[i] += p;

    ax[i] += fr*xx/r - fp*yy/r2;
    ay[i] += fr*yy/r + fp*xx/r2;
    az[i] += fz;
  }

#ifdef MPE_PROFILE
  MPE_Log_event(12, myid, "e_compute_force");
#endif

}


extern "C" 
void determine_fields_at_point_Cyl(double r, double z, double phi,
				   double *tdens, double *tpotl, 
				   double *tpotr, double *tpotz, double *tpotp)
{
  ortho->accumulated_eval(r, z, phi, *tpotl, *tpotr, *tpotz, *tpotp);
  *tdens = ortho->accumulated_dens_eval(r, z, phi);
}

				/* Dump coefficients to a file */

extern "C" void dump_coefs_Cyl(FILE *fout)
{
  /*
  fprintf(fout, "Time=%f\n", tnow);
  ortho->dump_coefs(fout);
  */
  ortho->dump_coefs_binary(fout, tnow);
}
