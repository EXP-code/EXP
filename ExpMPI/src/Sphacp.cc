// This may look like C code, but it is really -*- C++ -*-

// #define DEBUG 1

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  the spherical Sturm-Liouville basis
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

#define DENSITY
// #define SELECTOR

#include <values.h>

#include "expand.h"

#include <gaussQ.h>
#include <SphericalSL.h>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void determine_coefficients_Cyl(void);
void dump_mzero(char *, int);


//======================================================================
//======================================================================
//======================================================================

static SphericalSL *orthoS;
static int eof;

extern "C" {
  void determine_acceleration_and_potential_SLsph(void);

  void determine_fields_at_point_SLsphere(double r, double z, double phi,
					  double *tdens, double *tpotl, 
					  double *tpotr, double *tpotz, 
					  double *tpotp) {
    orthoS->determine_fields_at_point_SLsph(r, z, phi,
					    tdens, tpotl, tpotr, tpotz, tpotp);
  }

  void get_potl_dens_SLsph(int l, int n, double r) {
    orthoS->get_potl_dens_SLsph(l, n, r);
  }

  void get_dens_coefs_SLsph(int l, Vector& coef, double *p) {
    orthoS->get_dens_coefs_SLsph(l, coef, p);
  }
    
  void get_pot_coefs_SLsph(int l, Vector& coef, double *p, double *dp) {
    get_pot_coefs_SLsph(l, coef, p, dp);
  }

  void dump_coefs_SLsph(FILE *file) {
    orthoS->dump_coefs_SLsph(file);
  }
}



extern "C" 
void get_acceleration_and_potential_SLsph(void)
{
  static int firstime=1;

  if (firstime) {

    SphericalSL::RMAX = rsphSL;

    orthoS = new SphericalSL(lmax, nmax);

  }

  /*======================*/
  /* Compute coefficients */
  /*======================*/

  if (firstime || self_consistent) {
    
    orthoS->determine_coefficients_SLsph();

    firstime = 0;
  }


  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  orthoS->determine_acceleration_and_potential_SLsph();

  MPL_stop_timer();
}

