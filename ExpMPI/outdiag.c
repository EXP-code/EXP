/*
  diagnostic output: compute density, potential and d(pot)/dr
                     along the direction (PHI,COSTH)
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#define PHI 0.0			/* "Line of Sight" */
#define THETA (0.5*M_PI)
#define TOLR (1.0e-3)

				/* Function defs */

void determine_fields_at_point_bes(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);

void determine_fields_at_point_CB(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);

void determine_fields_at_point_CBDisk(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);

void determine_fields_at_point_HERNQ(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);

void determine_fields_at_point_Cyl(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);



#define NUM 100
void out_diag(int n)
{

  int i,ndigits;
  char fstring[20],file[80];
  static int nsnap=0;
  FILE *fout;
  double r,dr,dens;
  double potl,potr,pott,potp;
  double phi=PHI,theta=THETA;
    
  if (!bessel_sph && !c_brock && !c_brock_disk && !hernq && !cylinder) return;

  if (nsnap) 
    ndigits = (int)log10((double)nsnap)+1; 
  else
    ndigits = 1;

  sprintf(fstring,"%%s.%%s.%%%dd",ndigits);

  sprintf(file,fstring,outname,"diag",nsnap);
  if ( (fout=fopen(file,"w")) == NULL) {
    fprintf(stderr,"Couldn't open %s . . . quitting\n",file);
    exit(-1);
  }
  
  fprintf(fout,"! 1) Radius; 2) Rho; 3) Pot; 4) d(Pot)/dr; 5) d(Pot)/d cos(theta); 6) d(Pot)/d phi\n");

    /* Determine potential and acceleration */

  dr = (rdiag-TOLR)/(double)NUM;
  for (i=0; i<=NUM; i++) {
    r = TOLR + dr*i;

    if (bessel_sph)
      determine_fields_at_point_bes(r, theta, phi, 
				    &dens, &potl, &potr, &pott, &potp);
    else if (c_brock)
      determine_fields_at_point_CB(r, theta, phi, 
				    &dens, &potl, &potr, &pott, &potp);
    else if (c_brock_disk)
      determine_fields_at_point_CBDisk(r, theta, phi, 
				    &dens, &potl, &potr, &pott, &potp);
    else if (hernq)
      determine_fields_at_point_HERNQ(r, theta, phi, 
				    &dens, &potl, &potr, &pott, &potp);
    else
      dens = potl = potr = pott = potp = 0.0;

    fprintf(fout," %e %e %e %e %e %e",r,dens,potl,potr,pott,potp);

    if (cylinder) {
      determine_fields_at_point_Cyl(r, 0.0, phi, 
				    &dens, &potl, &potr, &pott, &potp);
      fprintf(fout," %e %e %e %e %e %e",r,dens,potl,potr,pott,potp);
    }

    fprintf(fout,"\n");
  }

  fclose(fout);

  nsnap=nsnap+1;

}

