/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Update coefficient file with coefficients at current time
 *
 *
 *  Notes:
 *  -----
 *  Output is a binary file consisting of doubles.  At each call,
 *  the current time is recorded followed by the entire expcoef array
 *  [dimensions: lmax*(lmax+2) by nmax] with the n dimension changing
 *  most rapidly
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *
 ***************************************************************************/

/*
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void dump_coefs_CB(FILE *);
void dump_coefs_CBDisk(FILE *);
void dump_coefs_HERNQ(FILE *);
void dump_coefs_bes(FILE *);
void dump_coefs_SLsph(FILE *);
void dump_coefs_Cyl(FILE *);
void dump_coefs_slab(FILE *);
void dump_coefs_slabSL(FILE *);

void out_coef(int n)
{
  FILE *fcoef;
				/* Open coefficient file for append */

  if ( (fcoef=fopen(coeffile,"a+")) == NULL) {
    fprintf(stderr,"Couldn't open %s . . . quitting\n",coeffile);
    exit(-1);
  }

  if (c_brock)
    dump_coefs_CB(fcoef);
  else if (c_brock_disk)
    dump_coefs_CBDisk(fcoef);
  else if (hernq)
    dump_coefs_HERNQ(fcoef);
  else if (bessel_sph)
    dump_coefs_bes(fcoef);
  else if (sphereSL)
    dump_coefs_SLsph(fcoef);
  else if (slab) {
    if (slabSL)
      dump_coefs_slabSL(fcoef);
    else
      dump_coefs_slab(fcoef);
  }

  fclose(fcoef);
}

void out_coef_cyl(int n)
{
  FILE *fcoef;
				/* Open coefficient file for append */

  if ( (fcoef=fopen(coeffilecyl,"a+")) == NULL) {
    fprintf(stderr,"Couldn't open %s . . . quitting\n",coeffilecyl);
    exit(-1);
  }

  dump_coefs_Cyl(fcoef);

  fclose(fcoef);
}

