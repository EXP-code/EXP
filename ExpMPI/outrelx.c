/*
  diagnostic output: compute relaxation diagnostics
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void out_relx(int n)
{

  int i;
  static char file[80];
  static int firstime=1;
  FILE *fout;
  double delta_e, e_average=0.0, e_absolute=0.0, variance=0.0;
  double e_average1=0.0, e_absolute1=0.0, variance1=0.0;
  int used1 = 0, used0 = 0;
  

#ifdef MPE_PROFILE
  MPE_Log_event(15, myid, "b_relx");
#endif

  if (firstime) {

    distribute_relx();

    if (myid==0) {
				/* Initialize output file */
      sprintf(file,"%s.%s", "relx", outname);
      if ( (fout=fopen(file,"w")) == NULL) {
	fprintf(stderr,"Couldn't open %s . . . quitting\n", file);
	exit(-1);
    }
      fprintf(fout,"! 1) time 2) step 3) Delta E; 4) Root variance; 5) |Delta E|\n");
      fclose(fout);
    }
    
    firstime = 0;

    return;
  }

  if (myid != 0) {

    for (i=1; i<=nbodies; i++) {
      if (freeze_particle(i)) continue;

      delta_e = 0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]) + 
	mass[i]*(pot[i] + potext[i]) - esave[i];
	
      e_average1 += delta_e/esave[i];
      e_absolute1 += fabs(delta_e/esave[i]);
      variance1 += (delta_e*delta_e/(esave[i]*esave[i]));
      used1++;
    }

  }

				/* Collect values */

  MPI_Reduce(&used1, &used0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&e_average1, &e_average, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&e_absolute1, &e_absolute, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&variance1, &variance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid == 0) {

    e_average /= used0;
    e_absolute /= used0;
    variance = (variance - e_average*e_average)/(used0 - 1);

    if ( (fout=fopen(file,"a+")) == NULL) {
      fprintf(stderr,"Couldn't reopen %s\n",file);
      return;
    }
    fprintf(fout, "%13.6e  %5d  %13.6e  %13.6e  %13.6e\n", tnow, n, e_average,
	    sqrt(variance), e_absolute);
    fclose(fout);
  }

#ifdef MPE_PROFILE
  MPE_Log_event(16, myid, "e_relx");
#endif

}


