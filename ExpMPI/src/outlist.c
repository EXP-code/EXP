/*
  spit out phase-space, one particle per line
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void out_list(int n)
{
  struct Partstruct
  {
    int component;
    double mass;
    double pos[3];
    double vel[3];
    double pot;
    double potext;
    double esave;
    double mfp;
  } *p;
  int number=-1;

  int i,ndigits;
  double lx,ly,lz,ek,ep,energy;
  char fstring[20],file[80];
  static int nsnap=0;
  FILE *fout1=NULL, *fout2=NULL;


  if (myid==0) {

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

  }

  while (1) {

    get_particles(&number, p);

    if (number == 0) break;

    if (myid==0) {

      for (i=0; i<number; i++) {

	if (outbods) 
	  fprintf(fout1,"%d %e %e %e %e %e %e %e %e",
		  p[i].component, p[i].mass, 
		  p[i].pos[0], p[i].pos[1], p[i].pos[2],
		  p[i].vel[0], p[i].vel[1], p[i].vel[2],
		  p[i].pot, p[i].potext);
	
	if (outbods) fprintf(fout1,"\n");

	if (outengy) {
	  lx=p[i].mass*(p[i].pos[1]*p[i].vel[2]-p[i].pos[2]*p[i].vel[1]);
	  ly=p[i].mass*(p[i].pos[2]*p[i].vel[0]-p[i].pos[0]*p[i].vel[2]);
	  lz=p[i].mass*(p[i].pos[0]*p[i].vel[1]-p[i].pos[1]*p[i].vel[0]);
	  ek=0.5*p[i].mass*(p[i].vel[0]*p[i].vel[0] + 
			    p[i].vel[1]*p[i].vel[1] +
			    p[i].vel[2]*p[i].vel[2]);
	  ep=p[i].mass*p[i].pot+p[i].mass*p[i].potext;
	  energy=ek+ep;
	  fprintf(fout2,"%e %e %e %e %e\n",tnow,lx,ly,lz,energy);
	}
      }
    }

  }

  if (myid==0) {
    
    if (outbods) fclose(fout1);
    if (outengy) fclose(fout2);
    
    nsnap=nsnap+1;

  }

}
