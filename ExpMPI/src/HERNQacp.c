/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  the Clutton-Brock ultraspherical biorthogonal expansion for the Hernquist
 *  model
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

#include "expand.h"
#include "exp_thread.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void get_ultra(int nmax, double l, double x, double *u);
void get_dpotl_HERNQ(int lmax, int nmax, double r, double **p, double **dp);
void get_dpotl_HERNQ_safe(int lmax, int nmax, double r, double **p, double **dp,
			  double *uu, double *duu);
void get_potl_HERNQ(int lmax, int nmax, double r, double **p);
void get_potl_HERNQ_safe(int lmax, int nmax, double r, double **p, double *uu);
void get_dens_HERNQ(int lmax, int nmax, double r, double **p);
void get_potl_dens_HERNQ(int lmax, int nmax, double r, double **p, double **d);
double normHERNQ(int,int);
double knlHERNQ(int,int);
void legendre(int lmax, double x, double ***pp);
void dlegendre(int lmax, double x, double ***pp, double ***dpp);
void sinecosine(int mmax, double phi, double **cc, double **ss);
void legendre_R(int lmax, double x, double **p);
void dlegendre_R(int lmax, double x, double **p, double **dp);
void sinecosine_R(int mmax, double phi, double *c, double *s);
void get_pot_coefs_HERNQ(int, double *, double *, double *);
void get_pot_coefs_HERNQ_safe(int, double *, double *, double *,
			      double **, double **);
void get_dens_coefs_HERNQ(int, double *, double *);
void determine_coefficients_HERNQ(void);
void determine_acceleration_and_potential_HERNQ(void);
void determine_fields_at_point_HERNQ(double r, double theta, double phi,
   double *tdens, double *tpotl, double *tpotr, double *tpott, double *tpotp);


static double **krnl, **u, **du, **dend, ***potd, ***dpot;
static double **factorial,**cosm,**sinm,***legs,***dlegs;
static double ***expcoef0;
static int Lmax;

void test_eigen(char *);
void pca_hall(int compute);

void get_acceleration_and_potential_HERNQ(void)
{
  static int firstime=1;
  int l, m, n, i;

  if (firstime) {
    
    Lmax = lmax<1 ? 1 : lmax;

    if (nthrds<1) nthrds=1;


    /*	Allocate coefficient matrix                                     */
    expcoef = dmatrix(0,Lmax*(Lmax+2),1,nmax);
    expcoef1 = dmatrix(0,Lmax*(Lmax+2),1,nmax);

    expcoef0  = (double ***)malloc(nthrds*sizeof(double**));
    if (!expcoef0) {
      fprintf(stderr, "HERNQacp: problem allocating <expcoef0>\n");
      exit(-1);
    }
    for (i=0; i<nthrds; i++)
      expcoef0[i] = dmatrix(0, Lmax*(Lmax+2), 1, nmax);


    if (selector) {
      cc = (double ***)malloc((unsigned) (Lmax*(Lmax+2)+1)*sizeof(double **));
      if (!cc) {
	fprintf(stderr, "HERNQacp: problem allocating <cc>\n");
	exit(-1);
      }
/*      cc -= 1; */
      for (l=0; l<=Lmax*(Lmax+2); l++)
	cc[l] = dmatrix(1, nmax, 1, nmax);
  
      cc1 = (double ***)malloc((unsigned) (Lmax*(Lmax+2)+1)*sizeof(double **));
      if (!cc1) {
	fprintf(stderr, "HERNQacp: problem allocating <cc1>\n");
	exit(-1);
      }
/*      cc1 -= 1; */
      for (l=0; l<=Lmax*(Lmax+2); l++)
	cc1[l] = dmatrix(1, nmax, 1, nmax);
  
    }

    /* Allocate and compute normalization matrix */

    normM = dmatrix(0,Lmax,1,nmax);
    krnl  = dmatrix(0,Lmax,1,nmax);
    dend  = dmatrix(0,Lmax,1,nmax);


    potd  = (double ***)malloc(nthrds*sizeof(double**));
    if (!potd) {
      fprintf(stderr, "HERNQacp: problem allocating <potd>\n");
      exit(-1);
    }

    dpot  = (double ***)malloc(nthrds*sizeof(double**));
    if (!dpot) {
      fprintf(stderr, "HERNQacp: problem allocating <dpot>\n");
      exit(-1);
    }

    for (i=0; i<nthrds; i++) {
      potd[i]  = dmatrix(0,Lmax,1,nmax);
      dpot[i]  = dmatrix(0,Lmax,1,nmax);
    }

    /* Sin, cos, legendre */

    cosm = (double **)malloc(nthrds*sizeof(double*));
    if (!cosm) {
      fprintf(stderr, "HERNQacp: problem allocating <cosm>\n");
      exit(-1);
    }
    sinm = (double **)malloc(nthrds*sizeof(double*));
    if (!sinm) {
      fprintf(stderr, "HERNQacp: problem allocating <sinm>\n");
      exit(-1);
    }
    legs = (double ***)malloc(nthrds*sizeof(double**));
    if (!legs) {
      fprintf(stderr, "HERNQacp: problem allocating <legs>\n");
      exit(-1);
    }
    dlegs = (double ***)malloc(nthrds*sizeof(double**));
    if (!dlegs) {
      fprintf(stderr, "HERNQacp: problem allocating <dlegs>\n");
      exit(-1);
    }

    for (i=0; i<nthrds; i++) {
      cosm[i] = dvector(0,Lmax);
      sinm[i] = dvector(0,Lmax);
      legs[i] = dmatrix(0,Lmax,0,Lmax);
      dlegs[i] = dmatrix(0,Lmax,0,Lmax);
    }

				/* Work vectors */
    u = (double **)malloc(nthrds*sizeof(double*));
    du = (double **)malloc(nthrds*sizeof(double*));
    if (!u) {
      fprintf(stderr, "HERNQacp: problem allocating <u>\n");
      exit(-1);
    }
    if (!du) {
      fprintf(stderr, "HERNQacp: problem allocating <du>\n");
      exit(-1);
    }

    for (i=0; i<nthrds; i++) {
      u[i]  = dvector(0,nmax);
      du[i] = dvector(0,nmax);
    }

    for (l=0; l<=lmax; l++) {
      for (n=1; n<=nmax; n++) {
	normM[l][n] = normHERNQ(n-1,l);
	krnl[l][n] = knlHERNQ(n-1,l);
      }
    }
    
    /* Factorial matrix */

    factorial = dmatrix(0,Lmax,0,Lmax);

    for (l=0; l<=lmax; l++) {
      for (m=0; m<=l; m++) 
	factorial[l][m] = factrl(l-m)/factrl(l+m);
    }

  }


  /*======================*/
  /* Compute coefficients */
  /*======================*/

  if (firstime || self_consistent) {
    firstime = 0;
    determine_coefficients_HERNQ();
  }

  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  determine_acceleration_and_potential_HERNQ();

  MPL_stop_timer();

}

/* static pthread_mutex_t expc_lock, cc_lock; */
static pthread_mutex_t cc_lock;
static int *use, compute;

void * determine_coefficients_HERNQ_thread(void * arg)
{
  int l, i, loffset, moffset, m, n, nn;
  double r, r2, rs, fac1, fac2, costh, phi;
  double facs1=0.0, facs2=0.0, fac0=4.0*M_PI;
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
	
      
      legendre_R(lmax, costh, legs[id]);
      sinecosine_R(lmax, phi, cosm[id], sinm[id]);
      /*
	legendre(lmax, costh, &legs[id]);
	sinecosine(lmax, phi, &cosm[id], &sinm[id]);
      */

      get_potl_HERNQ_safe(lmax, nmax, rs, potd[id], u[id]);
      
      /*		l loop */
      for (l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {
	/*		m loop */
	for (m=0, moffset=0; m<=l; m++) {
	  if (m==0) {
	    if (selector && compute)
	      facs1 = legs[id][l][m]*legs[id][l][m]*mass[i];
	    for (n=1; n<=nmax; n++) {
	      expcoef0[id][loffset+moffset][n] += potd[id][l][n]*legs[id][l][m]*mass[i]*
		fac0/normM[l][n];

	      if (selector && compute) {
		pthread_mutex_lock(&cc_lock);
		for (nn=n; nn<=nmax; nn++)
		  cc1[loffset+moffset][n][nn] += potd[id][l][n]*potd[id][l][nn]*
		    facs1/(normM[l][n]*normM[l][nn]);
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
	      expcoef0[id][loffset+moffset][n] += potd[id][l][n]*fac1*mass[i]*
		fac0/normM[l][n];
	      expcoef0[id][loffset+moffset+1][n] += potd[id][l][n]*fac2*mass[i]*
		fac0/normM[l][n];

	      if (selector && compute) {
		pthread_mutex_lock(&cc_lock);
		for (nn=n; nn<=nmax; nn++) {
		  cc1[loffset+moffset][n][nn] += 
		    potd[id][l][n]*potd[id][l][nn]*facs1/(normM[l][n]*normM[l][nn]);
		  cc1[loffset+moffset+1][n][nn] +=
		    potd[id][l][n]*potd[id][l][nn]*facs2/(normM[l][n]*normM[l][nn]);
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

  return (NULL);
}

void determine_coefficients_HERNQ(void)
{
  static int firstime=1;
  static char routine[] = "determine_coefficients_HERNQ";
  int l, i, loffset, moffset, m, n, nn, use0, use1;

  if (selector) compute = !(this_step%npca) || firstime;

  /*		Clean */
  for (n=1; n<=nmax; n++) {
      for (l=0; l<=lmax*(lmax+2); l++) {
	expcoef[l][n] = 0.0;
	expcoef1[l][n] = 0.0;
	for (i=0; i<nthrds; i++) expcoef0[i][l][n] = 0.0;
	if (selector && compute) 
	  for (nn=n; nn<=nmax; nn++) cc1[l][n][nn] = 0.0;
      }
    }

  use0 = 0;
  use1 = 0;

  use = (int *) malloc(nthrds*sizeof(int));
  if (!use) {
    fprintf(stderr, "HERNQacp: problem allocating <use>\n");
    exit(-1);
  }

				/* Initialize locks */
  /*    make_mutex(&cc_lock, routine, "cc_lock"); */

  exp_thread_fork(determine_coefficients_HERNQ_thread, routine);

  /*    kill_mutex(&cc_lock, routine, "cc_lock"); */

  for (i=0; i<nthrds; i++) {
    use1 += use[i];
    for (l=0; l<= Lmax*(Lmax+2); l++)
      for (n=1; n<=nmax; n++)
	expcoef1[l][n] += expcoef0[i][l][n];
  }
    
  free(use);

  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (myid==0) used += use0;


  if (!selector) {

    for (l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {
	/*		m loop */
	for (m=0, moffset=0; m<=l; m++) {
	  if (m==0) {

	    MPI_Allreduce ( &expcoef1[loffset+moffset][1],
			   &expcoef[loffset+moffset][1],
			   nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	    moffset++;
	  }
	  else {

	    MPI_Allreduce ( &expcoef1[loffset+moffset][1],
			   &expcoef[loffset+moffset][1],
			   nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	    MPI_Allreduce ( &expcoef1[loffset+moffset+1][1],
			   &expcoef[loffset+moffset+1][1],
			   nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    moffset+=2;
	  }
	}
      }
  }
 
  if (selector) {

    parallel_gather_coefficients();

    if (myid == 0) pca_hall(compute);

    parallel_distribute_coefficients();

    firstime = 0;
  }

}

void  * determine_acceleration_and_potential_HERNQ_thread(void * arg)
{
  int l,i,loffset,moffset,m;
  double r,rs,fac,fac1,fac2,fac3,fac4,costh,phi,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps;
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
    phi = atan2(yy,xx);

    dlegendre_R(lmax, costh, legs[id], dlegs[id]);
    sinecosine_R(lmax, phi, cosm[id], sinm[id]);
    /*
      dlegendre(lmax, costh, &legs[id], &dlegs[id]);
      sinecosine(lmax, phi, &cosm[id], &sinm[id]);
      */
    get_dpotl_HERNQ_safe(lmax, nmax, rs, potd[id], dpot[id], u[id], du[id]);


    get_pot_coefs_HERNQ_safe(0, expcoef[0], &p, &dp, potd[id], dpot[id]);
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
	    get_pot_coefs_HERNQ_safe(l, expcoef[loffset+moffset], &p, &dp,
				  potd[id], dpot[id]);
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
	    get_pot_coefs_HERNQ_safe(l,expcoef[loffset+moffset],&pc,&dpc,
				  potd[id],dpot[id]);
	    get_pot_coefs_HERNQ_safe(l,expcoef[loffset+moffset+1],&ps,&dps,
				  potd[id],dpot[id]);
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

  return (NULL);
}


void determine_acceleration_and_potential_HERNQ(void)
{
  static char routine[] = "determine_acceleration_and_potential_HERNQ";

#ifdef MPE_PROFILE
  MPE_Log_event(11, myid, "b_compute_force");
#endif

  exp_thread_fork(determine_acceleration_and_potential_HERNQ_thread, routine);

#ifdef MPE_PROFILE
  MPE_Log_event(12, myid, "e_compute_force");
#endif

}


void determine_fields_at_point_HERNQ(double r, double theta, double phi,
   double *tdens, double *tpotl, double *tpotr, double *tpott, double *tpotp)
{
  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps,dens;
  double dfac=0.25/M_PI;


  rs = r/scale;
  costh = cos(theta);

  fac1 = dfac;

  dlegendre_R(lmax, costh, legs[0], dlegs[0]);
  sinecosine_R(lmax, phi, cosm[0], sinm[0]);
  get_dens_HERNQ(lmax, nmax, rs, dend);
  get_dpotl_HERNQ(lmax, nmax, rs, potd[0], dpot[0]);

  get_dens_coefs_HERNQ(0,expcoef[0],&dens);
  dens *= dfac*dfac;

  get_pot_coefs_HERNQ(0,expcoef[0],&p,&dp);
  potl = fac1*p;
  potr = fac1*dp;
  pott = potp = 0.0;
  
      
  /*		l loop */
    
  for (l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
    
    /*		m loop */
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
	fac2 = fac1*legs[0][l][m];
	get_dens_coefs_HERNQ(l,expcoef[loffset+moffset],&p);
	dens += dfac*fac2*p;
	get_pot_coefs_HERNQ(l,expcoef[loffset+moffset],&p,&dp);
	potl += fac2*p;
	potr += fac2*dp;
	pott += fac1*dlegs[0][l][m]*p;
	moffset++;
      }
      else {
	fac2 = 2.0 * fac1 * factorial[l][m];
	fac3 = fac2 * legs[0][l][m];
	fac4 = fac2 * dlegs[0][l][m];
	
	get_dens_coefs_HERNQ(l,expcoef[loffset+moffset],&pc);
	get_dens_coefs_HERNQ(l,expcoef[loffset+moffset+1],&ps);
	dens += dfac*fac3*(pc*cosm[0][m] + ps*sinm[0][m]);
	
	get_pot_coefs_HERNQ(l,expcoef[loffset+moffset],&pc,&dpc);
	get_pot_coefs_HERNQ(l,expcoef[loffset+moffset+1],&ps,&dps);
	potl += fac3*(pc*cosm[0][m] + ps*sinm[0][m]);
	potr += fac3*(dpc*cosm[0][m] + dps*sinm[0][m]);
	pott += fac4*(pc*cosm[0][m] + ps*sinm[0][m]);
	potp += fac3*(-pc*sinm[0][m] + ps*cosm[0][m])*m;
	moffset +=2;
      }
    }
  }

  *tdens = dens/(scale*scale*scale);
  *tpotl = potl/scale;
  *tpotr = potr/(scale*scale);
  *tpott = pott/scale;
  *tpotp = potp/scale;
  
}


void get_pot_coefs_HERNQ(int l, double *coef, double *p, double *dp)
{
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (i=1; i<=nmax; i++) {
    pp  += potd[0][l][i] * coef[i];
    dpp += dpot[0][l][i] * coef[i];
  }

  *p = -pp;
  *dp = -dpp;
}


void get_pot_coefs_HERNQ_safe(int l, double *coef, double *p, double *dp,
			      double **potd1, double **dpot1)
{
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (i=1; i<=nmax; i++) {
    pp  += potd1[l][i] * coef[i];
    dpp += dpot1[l][i] * coef[i];
  }

  *p = -pp;
  *dp = -dpp;
}


void get_dens_coefs_HERNQ(int l, double *coef, double *p)
{
  double pp;
  int i;

  pp = 0.0;

  for (i=1; i<=nmax; i++)
    pp  += dend[l][i] * coef[i];

  *p = pp;
}


				/* Get potential functions by recursion */

void get_dpotl_HERNQ(int lmax, int nmax, double r, double **p, double **dp)
{
  double x, dx, fac, rfac, drfac, dfac1, dfac2;
  double r_to_rq(double), d_r_to_rq(double);
  int l, n;

  x = r_to_rq(r);
  dx = d_r_to_rq(r);

  fac = 0.25*(1.0 - x*x);
  rfac = 0.5*(1.0 - x);
  drfac = -1.0/(1.0 - x*x);
  
  for (l=0; l<=lmax; l++) {
    dfac1 = 1.0 + x + 2.0*x*l;
    dfac2 = 4.0*l + 3.0;

    get_ultra(nmax-1, 2.0*l+0.5, x, u[0]);
    get_ultra(nmax-1, 2.0*l+1.5, x, du[0]);

    for (n=0; n<nmax; n++) p[l][n+1] = rfac*u[0][n];
    dp[l][1] = dx*drfac*dfac1*rfac*u[0][0];
    for (n=1; n<nmax; n++) dp[l][n+1] = dx*rfac*(drfac*dfac1*u[0][n] + dfac2*du[0][n-1]);

    rfac *= fac;
  }

/* check */
/*  {
    double tmp, dtmp, del, potl_HERNQ(int, int, double);

    del = 1.0e-4;
    for (l=0; l<=lmax; l++) {
      for (n=1; n<nmax; n++) {
	tmp = potl_HERNQ(n, l, x);
	dtmp = dx*(potl_HERNQ(n, l, x+del) - potl_HERNQ(n, l, x-del))/(2.0*del);
	printf("%d %d> %e %e %e %e\n",l, n, tmp, p[l][n], dtmp, dp[l][n]);
      }
    }
  } */
}

void get_dpotl_HERNQ_safe(int lmax, int nmax, double r, double **p, double **dp,
			  double *uu, double *duu)
{
  double x, dx, fac, rfac, drfac, dfac1, dfac2;
  double r_to_rq(double), d_r_to_rq(double);
  int l, n;

  x = r_to_rq(r);
  dx = d_r_to_rq(r);

  fac = 0.25*(1.0 - x*x);
  rfac = 0.5*(1.0 - x);
  drfac = -1.0/(1.0 - x*x);
  
  for (l=0; l<=lmax; l++) {
    dfac1 = 1.0 + x + 2.0*x*l;
    dfac2 = 4.0*l + 3.0;

    get_ultra(nmax-1, 2.0*l+0.5, x, uu);
    get_ultra(nmax-1, 2.0*l+1.5, x, duu);

    for (n=0; n<nmax; n++) p[l][n+1] = rfac*uu[n];
    dp[l][1] = dx*drfac*dfac1*rfac*uu[0];
    for (n=1; n<nmax; n++) dp[l][n+1] = dx*rfac*(drfac*dfac1*uu[n] + dfac2*duu[n-1]);

    rfac *= fac;
  }

/* check */
/*  {
    double tmp, dtmp, del, potl_HERNQ(int, int, double);

    del = 1.0e-4;
    for (l=0; l<=lmax; l++) {
      for (n=1; n<nmax; n++) {
	tmp = potl_HERNQ(n, l, x);
	dtmp = dx*(potl_HERNQ(n, l, x+del) - potl_HERNQ(n, l, x-del))/(2.0*del);
	printf("%d %d> %e %e %e %e\n",l, n, tmp, p[l][n], dtmp, dp[l][n]);
      }
    }
  } */
}

void get_potl_HERNQ(int lmax, int nmax, double r, double **p)
{
  double x, fac, rfac;
  double r_to_rq(double);
  int l, n;

  x = r_to_rq(r);

  fac = 0.25*(1.0 - x*x);
  rfac = 0.5*(1.0 - x);
  
  for (l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, u[0]);
    for (n=0; n<nmax; n++) p[l][n+1] = rfac*u[0][n];
    rfac *= fac;
  }

/* check */
/*  {
    double tmp, potl_HERNQ(int, int, double);

    printf("\nx=%e\n",x);
    for (l=0; l<=lmax; l++) {
      for (n=1; n<nmax; n++) {
	tmp = potl_HERNQ(n, l, x);
	printf("%d %d> %e %e\n",l, n, tmp, p[l][n]);
      }
    }
  } */

}

void get_potl_HERNQ_safe(int lmax, int nmax, double r, double **p, double *uu)
{
  double x, fac, rfac;
  double r_to_rq(double);
  int l, n;

  x = r_to_rq(r);

  fac = 0.25*(1.0 - x*x);
  rfac = 0.5*(1.0 - x);
  
  for (l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, uu);
    for (n=0; n<nmax; n++) p[l][n+1] = rfac*uu[n];
    rfac *= fac;
  }
}

void get_dens_HERNQ(int lmax, int nmax, double r, double **p)
{
  double x, fac, rfac;
  double r_to_rq(double);
  int l, n;

  x = r_to_rq(r);
  fac = 0.25*(1.0 - x*x);
  rfac = 0.25*pow(1.0 - x, 5.0)/(1.0 - x*x);

  for (l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, u[0]);
    for (n=0; n<nmax; n++) p[l][n+1] = krnl[l][n+1]*rfac*u[0][n];
    rfac *= fac;
  }

}

void get_potl_dens_HERNQ(int lmax, int nmax, double r, double **p, double **d)
{
  double x, fac, rfacp, rfacd;
  double r_to_rq(double);
  int l, n;

  x = r_to_rq(r);

  fac = 0.25*(1.0 - x*x);
  rfacp = 0.5*(1.0 - x);
  rfacd = 0.25*pow(1.0 - x, 5.0)/(1.0 - x*x);


  for (l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, u[0]);
    for (n=0; n<nmax; n++) p[l][n+1] = rfacp*u[0][n];
    for (n=0; n<nmax; n++) d[l][n+1] = krnl[l][n+1]*rfacd*u[0][n];
    rfacp *= fac;
    rfacd *= fac;
  }
}

				/* Dump coefficients to a file */

void dump_coefs_HERNQ(FILE *fout)
{
  int ir, l;

  fwrite(&tnow, sizeof(double), 1, fout);

  for (ir=1; ir<=nmax; ir++) {
    for (l=0; l<=lmax*(lmax+2); l++)
      fwrite(&expcoef[l][ir], sizeof(double), 1, fout);
  }

}

