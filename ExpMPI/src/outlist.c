/*
  spit out phase-space, one particle per line
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void out_list(int n)
{

  int i,ndigits;
  double lx,ly,lz,ek,ep,energy;
  char fstring[20],file[80];
  static int nsnap=0;
  FILE *fout1=NULL, *fout2=NULL;

  
  if (nsnap) 
    ndigits = (int)log10((double)nsnap)+1; 
  else
    ndigits = 1;

  sprintf(fstring,"%%s.%%s.%%%dd",ndigits);

  if (outbods) {
    sprintf(file,fstring,outname,"bods",nsnap);
    if ( (fout1=fopen(file,"w")) == NULL) {
      fprintf(stderr,"Couldn't open %s . . . quitting\n",file);
      exit(-1);
    }
  }

  if (outengy) {
    sprintf(file,fstring,outname,"engy",nsnap);
    if ( (fout2=fopen(file,"w")) == NULL) {
      fprintf(stderr,"Couldn't open %s . . . quitting\n",file);
      exit(-1);
    }
  }

  
  if (outbods) fprintf(fout1,"%d %e\n",nbodies,tnow);

  for (i=1; i<=nbodies; i++) {
    if (outbods) fprintf(fout1,"%d %e %e %e %e %e %e %e %e",
	    component[i],
	    mass[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],pot[i]+potext[i]);

    if (outbods) fprintf(fout1,"\n");

    if (outengy) {
      lx=mass[i]*(y[i]*vz[i]-z[i]*vy[i]);
      ly=mass[i]*(z[i]*vx[i]-x[i]*vz[i]);
      lz=mass[i]*(x[i]*vy[i]-y[i]*vx[i]);
      ek=0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
      ep=mass[i]*pot[i]+mass[i]*potext[i];
      energy=ek+ep;
      fprintf(fout2,"%e %e %e %e %e\n",tnow,lx,ly,lz,energy);
    }
  }

  if (outbods) fclose(fout1);
  if (outengy) fclose(fout2);

  nsnap=nsnap+1;

}
