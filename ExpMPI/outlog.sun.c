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
 *  7        m2tw,clausius,m2claus       -2T/W, Virial of Clausius (VC), -2T/VC
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

static char rcsid[] = "$Id$";

void return_stats(double *, double *, double *, double *, double *, double *,
		  double *, double *, double *, double *);

void out_log(int n)
{
  int i;
  double lxtot,lytot,lztot,mtot,vxcm,vycm,vzcm,etot,ektot,eptot;
  double m2tw,t1,clausius,m2claus,xcm,ycm,zcm,epselfg,cpu;
  double vrel,vrel2,r2stars,r2clds,v2stars,v2clds,masstar,mcld,mcld1,mcld2;
  double Mxx, Myy, Mzz, Mxy, Mxz, Myz;
  long curcpu,clock(void);
  static long lastcpu;
  static int laststep,firstime=1;
  FILE *flog;

				/* Use clock() to time step */
  cpu = 0.0;
  if (firstime) {
    lastcpu = clock();
    laststep = n;
    firstime = 0;
  }
  else if (n>laststep) {
    curcpu = clock();
    cpu = (curcpu-lastcpu)*1.0e-6/(n-laststep);
    lastcpu = curcpu;
    laststep = n;
  }

				/* Open logfile for append */

  if ( (flog=fopen(logfile,"a+")) == NULL) {
    fprintf(stderr,"Couldn't open %s . . . quitting\n",logfile);
    exit(-1);
  }

				/* Compute phase-space diagnostics */
  mtot = 0.0;
  xcm = 0.0;
  ycm = 0.0;
  zcm = 0.0;
  vxcm = 0.0;
  vycm = 0.0;
  vzcm = 0.0;
  lxtot = 0.0;
  lytot = 0.0;
  lztot = 0.0;
  etot = 0.0;
  ektot = 0.0;
  eptot = 0.0;
  epselfg = 0.0;
  clausius = 0.0;

  Mxx = 0.0;			/* inertia tensor elements -KL 6/23/92*/
  Myy = 0.0;
  Mzz = 0.0;
  Mxy = 0.0;
  Mxz = 0.0;
  Myz = 0.0;
  
  for (i=1; i<=nbodies; i++) {
    if (freeze_particle(i)) continue;	/* ignore frozen particles -KL 6/23/92*/
    mtot = mtot + mass[i];
    xcm = xcm + mass[i]*x[i];
    ycm = ycm + mass[i]*y[i];
    zcm = zcm + mass[i]*z[i];
    vxcm = vxcm + mass[i]*vx[i];
    vycm = vycm + mass[i]*vy[i];
    vzcm = vzcm + mass[i]*vz[i];
    eptot = eptot + 0.5*mass[i]*pot[i] + mass[i]*potext[i];
    epselfg = epselfg + 0.5*mass[i]*pot[i];
    ektot = ektot + 0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
    clausius = clausius + mass[i]*(x[i]*ax[i]+y[i]*ay[i]+z[i]*az[i]);

    /* add only bound particles to L and M */
    if (0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]) + pot[i] <= 0.0)
    {
       lxtot = lxtot + mass[i]*(y[i]*vz[i]-z[i]*vy[i]);
       lytot = lytot + mass[i]*(z[i]*vx[i]-x[i]*vz[i]);
       lztot = lztot + mass[i]*(x[i]*vy[i]-y[i]*vx[i]);
       Mxx = Mxx + mass[i]*x[i]*x[i];
       Myy = Myy + mass[i]*y[i]*y[i];
       Mzz = Mzz + mass[i]*z[i]*z[i];
       Mxy = Mxy + mass[i]*x[i]*y[i];
       Mxz = Mxz + mass[i]*x[i]*z[i];
       Myz = Myz + mass[i]*y[i]*z[i];
    }
  }

  xcm = xcm/mtot;
  ycm = ycm/mtot;
  zcm = zcm/mtot;
  vxcm = vxcm/mtot;
  vycm = vycm/mtot;
  vzcm = vzcm/mtot;
  
  etot = ektot+eptot;
  m2tw =  -2.0*ektot/eptot;
  m2claus =  -2.0*ektot/clausius;

  fprintf(flog,"%e\n", tnow);
  fprintf(flog,"%e %d %d %d\n", mtot,nbodies,ninteract,nmerge);
  fprintf(flog,"%e %e %e\n", xcm,ycm,zcm);
  fprintf(flog,"%e %e %e\n", vxcm,vycm,vzcm);
  fprintf(flog,"%e %e %e\n", lxtot,lytot,lztot);
  fprintf(flog,"%e %e %e %e\n", ektot,eptot,epselfg,etot);
  fprintf(flog,"%e %e %e\n", m2tw,clausius,m2claus);
  if (inertia) {
    fprintf(flog, "%e %e %e\n%e %e %e\n%e %e %e\n", Mxx, Mxy, Mxz,
	Mxy, Myy, Myz, Mxz, Myz, Mzz);
  }
  fprintf(flog,"%e %d\n", cpu,used);
  if (cloud) {
    return_stats(&vrel,&vrel2,&r2stars,&r2clds,&v2stars,&v2clds,
		 &masstar,&mcld,&mcld1,&mcld2);

    fprintf(flog,"%e %e %e %e %e %e %e %e %e %e\n",
	    vrel,vrel2,r2stars,r2clds,v2stars,v2clds,masstar,mcld,mcld1,mcld2);
  }
  else
    fprintf(flog,"\n");

  fprintf(flog,"\n");

  fclose(flog);
}
