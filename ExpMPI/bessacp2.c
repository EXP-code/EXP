#include "expand.h"
#include <exp_thread.h>

static char rcsid[] = "$Id$";

/******Definitions for routines <biorth_bes.c>******/

/*   Set up a structure for roots of a given order for biorthonormal 
     functions based on bessels */
struct ROOTS {
  int l;
  int n;
  double *a;
} *get_ln(int l, int n);

void rid_ln(struct ROOTS *p);
void set_radius(double R);
double dens(double r, int n, struct ROOTS *p),potl(double r, int n, struct ROOTS *p);

void get_dpotl_bes(int lmax, int nmax, double r, double **p, double **dp);
void get_potl_bes(int lmax, int nmax, double r, double **p);
void get_dens_bes(int lmax, int nmax, double r, double **p);
void get_potl_dens_bes(int lmax, int nmax, double r, double **p, double **d);
void legendre(int lmax, double x, double ***pp);
void dlegendre(int lmax, double x, double ***pp, double ***dpp);
void sinecosine(int mmax, double phi, double **cc, double **ss);
void get_pot_coefs_bes(int, double *, double *, double *);
void get_pot_coefs_bes_safe(int, double *, double *, double *, 
			    double **, double **);
void get_dens_coefs_bes(int, double *, double *);
void determine_coefficients_bes(void);
void determine_acceleration_and_potential_bes(void);
void determine_fields_at_point_bes(double r, double theta, double phi,
   double *tdens, double *tpotl, double *tpotr, double *tpott, double *tpotp);

extern double **krnl, *u, *du, **dend, ***potd, ***dpot;
extern double **factorial,**cosm,**sinm,***legs,***dlegs;

extern pthread_mutex_t expc_lock, cc_lock;
extern int *use, compute;

void determine_coefficients_bes_thread(void * arg)
{
  int l, i, loffset, moffset, m, n, nn, use1;
  double r, r2, rs, fac1, fac2, costh, phi;
  double facs1, facs2, fac0=4.0*M_PI;
  double xx, yy, zz;

  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  use[id] = 0;

  for (i=nbeg; i<=nend; i++) {

    if (freeze_particle(i)) continue;		/* frozen particles don't
						   contribute to field.
						   KL 5/27/92 */
    xx = x[i] - com1[0];
    yy = y[i] - com1[1];
    zz = z[i] - com1[2];

    r2 = (xx*xx + yy*yy + zz*zz);
    r = rr[i] = sqrt(r2) + DSMALL;
      
    if (component[i] != 1) continue;

    if (r<=rmax) {
      use[id]++;
      costh = zz/r;
      phi = atan2(yy,xx);
      rs = r/scale;
	
      
      legendre(lmax, costh, &legs[id]);
      sinecosine(lmax, phi, &cosm[id], &sinm[id]);

      get_potl_bes(lmax, nmax, rs, potd[id]);

      /*		l loop */
      for (l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {
	/*		m loop */
	for (m=0, moffset=0; m<=l; m++) {
	  if (m==0) {
	    if (selector && compute)
	      facs1 = legs[id][l][m]*legs[id][l][m]*mass[i];
	    for (n=1; n<=nmax; n++) {
	      pthread_mutex_lock(&expc_lock);
	      expcoef1[loffset+moffset][n] += potd[id][l][n]*legs[id][l][m]*mass[i]*
		fac0;
	      pthread_mutex_unlock(&expc_lock);

	      if (selector && compute) {
		pthread_mutex_lock(&cc_lock);
		for (nn=n; nn<=nmax; nn++)
		  cc1[loffset+moffset][n][nn] += potd[id][l][n]*potd[id][l][nn]*facs1;
		pthread_mutex_unlock(&cc_lock);
	      }
	    }
	    moffset++;
	  }
	  else {
	    fac1 = legs[id][l][m]*cosm[id][m];
	    fac2 = legs[id][l][m]*sinm[id][m];
	    if (selector && compute) {
	      facs1 = fac1*fac1*mass[i];
	      facs2 = fac2*fac2*mass[i];
	    }
	    for (n=1; n<=nmax; n++) {
	      pthread_mutex_lock(&expc_lock);
	      expcoef1[loffset+moffset][n] += potd[id][l][n]*fac1*mass[i]*fac0;
	      expcoef1[loffset+moffset+1][n] += potd[id][l][n]*fac2*mass[i]*fac0;
	      pthread_mutex_unlock(&expc_lock);

	      if (selector && compute) {
		pthread_mutex_lock(&cc_lock);
		for (nn=n; nn<=nmax; nn++) {
		  cc1[loffset+moffset][n][nn] += 
		    potd[id][l][n]*potd[id][l][nn]*facs1;
		  cc1[loffset+moffset+1][n][nn] +=
		    potd[id][l][n]*potd[id][l][nn]*facs2;
		}
		pthread_mutex_unlock(&cc_lock);
	      }
	      
	    }
	    moffset+=2;
	  }
	}
      }
    }
  }

}

void determine_acceleration_and_potential_bes_thread(void * arg)
{
  int l,i,loffset,moffset,m,ioff;
  double r,rs,r0,fac,fac1,fac2,fac3,fac4,costh,phi,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps,facp,facdp;
  double dfac=0.25/M_PI;
  double xx, yy, zz;

  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  for (i=nbeg; i<=nend; i++) {

    if (freeze_particle(i)) 	/* frozen particles do not respond */
      continue;			/* KL 5/27/92 */

    fac1 = dfac;

    xx = x[i] - com1[0];
    yy = y[i] - com1[1];
    zz = z[i] - com1[2];

    if (!self_consistent) {
      r = rr[i] = sqrt(xx*xx + yy*yy + zz*zz) + DSMALL;
    }
    else
      r = rr[i];


    costh = zz/r;
    rs = r/scale;
    phi = atan2(yy, xx);

    dlegendre(lmax, costh, &legs[id], &dlegs[id]);
    sinecosine(lmax, phi, &cosm[id], &sinm[id]);

    if (r>rmax) {
      ioff = 1;
      r0 = r;
      r = rmax;
      rs = r/scale;
    }
    else
      ioff = 0;


    get_dpotl_bes(lmax, nmax, rs, potd[id], dpot[id]);
    get_pot_coefs_bes_safe(0, expcoef[0], &p, &dp, potd[id], dpot[id]);
    if (ioff) {
      p *= rmax/r0;
      dp = -p/r0;
    }
    potl = fac1*p;
    potr = fac1*dp;
    pott = potp = 0.0;
      
    /*		l loop */
    
    for (l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {

      /*		m loop */
      for (m=0, moffset=0; m<=l; m++) {
	fac1 = (2.0*l+1.0)/(4.0*M_PI);
	if (m==0) {
				/* Suppress L=1 terms? */
	  if ( !NO_L1 || !(l==1) ) {
	    fac2 = fac1*legs[id][l][m];
	    get_pot_coefs_bes_safe(l, expcoef[loffset+moffset], &p, &dp,
			      potd[id], dpot[id]);
	    if (ioff) {
	      p *= pow(rmax/r0,(double)(l+1));
	      dp = -p/r0 * (l+1);
	    }
	    potl += fac2*p;
	    potr += fac2*dp;
	    pott += fac1*dlegs[id][l][m]*p;
	  }
	  moffset++;
	}
	else {
				/* Suppress L=1 terms? */
	  if ( !NO_L1 || !(l==1) ) {
	    fac2 = 2.0 * fac1 * factorial[l][m];
	    fac3 = fac2 * legs[id][l][m];
	    fac4 = fac2 * dlegs[id][l][m];
	    get_pot_coefs_bes_safe(l, expcoef[loffset+moffset], &pc, &dpc,
				   potd[id], dpot[id]);
	    get_pot_coefs_bes_safe(l, expcoef[loffset+moffset+1] ,&ps, &dps,
				   potd[id], dpot[id]);
	    if (ioff) {
	      facp = pow(rmax/r0,(double)(l+1));
	      facdp = -1.0/r0 * (l+1);
	      pc *= facp;
	      ps *= facp;
	      dpc = pc*facdp;
	      dps = ps*facdp;
	    }
	    potl += fac3*(pc*cosm[id][m] + ps*sinm[id][m]);
	    potr += fac3*(dpc*cosm[id][m] + dps*sinm[id][m]);
	    pott += fac4*(pc*cosm[id][m] + ps*sinm[id][m]);
	    potp += fac3*(-pc*sinm[id][m] + ps*cosm[id][m])*m;
	  }
	  moffset +=2;
	}
      }
    }

    fac = xx*xx + yy*yy;

    potr /= scale*scale;
    potl /= scale;
    pott /= scale;
    potp /= scale;

    ax[i] += -(potr*xx/r - pott*xx*zz/(r*r*r) );
    ay[i] += -(potr*yy/r - pott*yy*zz/(r*r*r) );
    az[i] += -(potr*zz/r + pott*fac/(r*r*r));
    if (fac > DSMALL2) {
      ax[i] +=  potp*yy/fac;
      ay[i] += -potp*xx/fac;
    }
    pot[i] += potl;

  }

}
