/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Update log file with current global values:
 *  10 lines per call
 *
 *  Output format:
 *  -------------
 *  LINE #   VALUES                      EXPLANATION
 *  ------   ------                      -----------
 *  1        tnow                        current time
 *  2        mtot,nbodies,ninteract,nmerge
 *                         current total mass, # bodies, # interact, # merge
 *  3        xcm,ycm,zcm                 center of mass
 *  4        vxcm,vycm,vzcm              center of mass velocity
 *  5        lxtot,lytot,lztot           total angular momentum
 *  6        ektot,eptot,epselfg,etot    KE, total PE, self PE only, total E
 *  7a       m2tw,clausius,m2claus       -2T/W, Virial of Clausius (VC), -2T/VC
 *  7b       m2tw,clausius,m2claus,m2c1,m2c2
 *                                       -2T/W, Virial of Clausius (VC), -2T/VC, -2T/VC|_1, -2T/VC|_2
 *  8        cpu,used                    mean cputime/step, # points on grid
 *  9        --blank--
 *  10       --blank--
 *
 *
 *  Notes:
 *  -----
 *  People have claimed that the "Virial of Clausius" is a better indicated
 *  of the true potential when computed by n-body scheme.  Suit yourself.
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

#ifdef MPI_LINUX
#include <time.h>
#endif

static char rcsid[] = "$Id$";

void return_stats(double *, double *, double *, double *, double *, double *,
		  double *, double *, double *, double *);

void out_log(int n)
{
  int i;
  double lxtot,lytot,lztot,mtot,vxcm,vycm,vzcm,etot,ektot,eptot;
  double m2tw,clausius,m2claus,xcm,ycm,zcm,epselfg,wtime;
  double lxtot1,lytot1,lztot1,mtot1,vxcm1,vycm1,vzcm1,etot1,ektot1,eptot1;
  double ektotc1, ektotc11, epself11, epself1, clausiusc1, clausiusc11;
  double ektotc2, ektotc21, epself21, epself2, clausiusc2, clausiusc21;
  double clausius1,xcm1,ycm1,zcm1,epselfg1;
  /*
  double vrel,vrel2,r2stars,r2clds,v2stars,v2clds,masstar,mcld,mcld1,mcld2;
  */
  double Mxx, Myy, Mzz, Mxy, Mxz, Myz;
  double Mxx1, Myy1, Mzz1, Mxy1, Mxz1, Myz1;
  static double curwtime, lastwtime;
  static double ektotxy=0.0;
  double ektotxy1=0.0;
  static int laststep,firstime=1;
  FILE *flog;

				/* Use clock() to time step */
  wtime = 0.0;
  if (firstime) {
    lastwtime = MPI_Wtime();
    laststep = n;
  }
  else if (n>laststep) {
    curwtime = MPI_Wtime();
    wtime = (curwtime-lastwtime)/(n-laststep);
    lastwtime = curwtime;
    laststep = n;
  }

				/* Compute phase-space diagnostics */
  mtot = 0.0;          mtot1 = 0.0;    
  xcm = 0.0;           xcm1 = 0.0;     
  ycm = 0.0;           ycm1 = 0.0;     
  zcm = 0.0;           zcm1 = 0.0;     
  vxcm = 0.0;          vxcm1 = 0.0;    
  vycm = 0.0;          vycm1 = 0.0;    
  vzcm = 0.0;          vzcm1 = 0.0;    
  lxtot = 0.0;         lxtot1 = 0.0;   
  lytot = 0.0;         lytot1 = 0.0;   
  lztot = 0.0;         lztot1 = 0.0;   
  etot = 0.0;          etot1 = 0.0;    
  ektot = 0.0;         ektot1 = 0.0;   
  eptot = 0.0;         eptot1 = 0.0;   
  ektotc1 = 0.0;       ektotc11 = 0.0;   
  ektotc2 = 0.0;       ektotc21 = 0.0;   
  epself1 = 0.0;       epself11 = 0.0;   
  epself2 = 0.0;       epself21 = 0.0;   
  epselfg = 0.0;       epselfg1 = 0.0; 
  clausius = 0.0;      clausius1 = 0.0;
  clausiusc1 = 0.0;    clausiusc11 = 0.0;
  clausiusc2 = 0.0;    clausiusc21 = 0.0;


			/* inertia tensor elements -KL 6/23/92*/
  Mxx = 0.0;           Mxx1 = 0.0;
  Myy = 0.0;	       Myy1 = 0.0;
  Mzz = 0.0;	       Mzz1 = 0.0;
  Mxy = 0.0;	       Mxy1 = 0.0;
  Mxz = 0.0;	       Mxz1 = 0.0;
  Myz = 0.0;	       Myz1 = 0.0;
  
				/* Collect info */
  if (myid > 0) {

    for (i=1; i<=nbodies; i++) {
				/* ignore frozen particles -KL 6/23/92*/
      if (freeze_particle(i)) continue;
      mtot1 = mtot1 + mass[i];
      xcm1 = xcm1 + mass[i]*x[i];
      ycm1 = ycm1 + mass[i]*y[i];
      zcm1 = zcm1 + mass[i]*z[i];
      vxcm1 = vxcm1 + mass[i]*vx[i];
      vycm1 = vycm1 + mass[i]*vy[i];
      vzcm1 = vzcm1 + mass[i]*vz[i];
      eptot1 = eptot1 + 0.5*mass[i]*pot[i] + mass[i]*potext[i];
      epselfg1 = epselfg1 + 0.5*mass[i]*pot[i];
      ektot1 = ektot1 + 0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
      if (slab && firstime)
	ektotxy1 = ektotxy1 + 0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]);
      clausius1 = clausius1 + mass[i]*(x[i]*ax[i]+y[i]*ay[i]+z[i]*az[i]);
      if (component[i]==1) {
	ektotc11 +=  0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
	clausiusc11 += mass[i]*(x[i]*ax[i]+y[i]*ay[i]+z[i]*az[i]);
	epself11 += 0.5*mass[i]*pot[i];
      }
      if (component[i]==2) {
	ektotc21 += 0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
	clausiusc21 += mass[i]*(x[i]*ax[i]+y[i]*ay[i]+z[i]*az[i]);
	epself21 += 0.5*mass[i]*pot[i];
      }
				/* add only bound particles to L and M */
      if (0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]) + pot[i] <= 0.0)
	{
	  lxtot1 = lxtot1 + mass[i]*(y[i]*vz[i]-z[i]*vy[i]);
	  lytot1 = lytot1 + mass[i]*(z[i]*vx[i]-x[i]*vz[i]);
	  lztot1 = lztot1 + mass[i]*(x[i]*vy[i]-y[i]*vx[i]);
	  Mxx1 = Mxx1 + mass[i]*x[i]*x[i];
	  Myy1 = Myy1 + mass[i]*y[i]*y[i];
	  Mzz1 = Mzz1 + mass[i]*z[i]*z[i];
	  Mxy1 = Mxy1 + mass[i]*x[i]*y[i];
	  Mxz1 = Mxz1 + mass[i]*x[i]*z[i];
	  Myz1 = Myz1 + mass[i]*y[i]*z[i];
	}
    }

  }
  
				/* Send back to Process 0 */

  MPI_Reduce(&mtot1, &mtot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&xcm1, &xcm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ycm1, &ycm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&zcm1, &zcm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&vxcm1, &vxcm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&vycm1, &vycm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&vzcm1, &vzcm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&eptot1, &eptot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&epselfg1, &epselfg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ektot1, &ektot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (slab && firstime)
    MPI_Reduce(&ektotxy1, &ektotxy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&clausius1, &clausius, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ektotc11, &ektotc1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ektotc21, &ektotc2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&epself11, &epself1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&epself21, &epself2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&clausiusc11, &clausiusc1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&clausiusc21, &clausiusc2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(&lxtot1, &lxtot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&lytot1, &lytot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&lztot1, &lztot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(&Mxx1, &Mxx, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Myy1, &Myy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Mzz1, &Mzz, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Mxy1, &Mxy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Mxz1, &Mxz, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Myz1, &Myz, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if (myid == 0) {
				/* Open logfile for append */

    if ( (flog=fopen(logfile,"a+")) == NULL) {
      fprintf(stderr,"Couldn't open %s . . . quitting\n",logfile);
    exit(-1);
    }

    xcm = xcm/mtot;
    ycm = ycm/mtot;
    zcm = zcm/mtot;
    vxcm = vxcm/mtot;
    vycm = vycm/mtot;
    vzcm = vzcm/mtot;
  
    if (slab) ektot -= ektotxy;
    
    etot = ektot+eptot;
    m2tw =  -2.0*ektot/eptot;
    m2claus =  -2.0*ektot/clausius;
    
    fprintf(flog,"%e\n", tnow);
    fprintf(flog,"%e %d %d %d\n", mtot,nbodies,ninteract,nmerge);
    fprintf(flog,"%e %e %e\n", xcm,ycm,zcm);
    fprintf(flog,"%e %e %e\n", vxcm,vycm,vzcm);
    fprintf(flog,"%e %e %e\n", lxtot,lytot,lztot);
    fprintf(flog,"%e %e %e %e", ektot,eptot,epselfg,etot);
    if (epself1<0.0 && epself2<0.0)
      fprintf(flog," %e %e\n", epself1, epself2);
    else
      fprintf(flog,"\n");
    fprintf(flog,"%e %e %e", m2tw,clausius,m2claus);
    if (clausiusc1<0.0 && clausiusc2<0.0) {
      double m2c1 = -2.0*ektotc1/clausiusc1;
      double m2c2 = -2.0*ektotc2/clausiusc2;
      
      /* test version */
      fprintf(flog," %e %e %e %e %e %e\n", 
	      ektotc1, clausiusc1, ektotc2, clausiusc2, m2c1, m2c2);
      /*
	fprintf(flog," %e %e\n", m2c1, m2c2); 
	*/
    }
    else
      fprintf(flog,"\n");
    fprintf(flog,"%e %d\n", wtime, used);
    if (inertia) {
      fprintf(flog, "%e %e %e\n%e %e %e\n%e %e %e\n", Mxx, Mxy, Mxz,
	      Mxy, Myy, Myz, Mxz, Myz, Mzz);
    }
    else
      fprintf(flog,"\n");
    fprintf(flog,"\n");

    fclose(flog);
  }

  firstime = 0;
}
