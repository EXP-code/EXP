/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  the Clutton-Brock flat disk expansion
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

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void get_dpotl_CBDisk(int lmax, int nmax, double r, double **p, double **dp);
void get_potl_CBDisk(int lmax, int nmax, double r, double **p);
void get_dens_CBDisk(int lmax, int nmax, double r, double **p);
void get_dens_potl_CBDisk(int lmax, int nmax, double r, double **p, double **d);
double normCBDisk(int,int);
void sinecosine(int mmax, double phi, double **cc, double **ss);
void get_pot_coefs_CBDisk(int, double *, double *, double *);
void get_dens_coefs_CBDisk(int, double *, double *);
void determine_coefficients_CBDisk(void);
void determine_acceleration_and_potential_CBDisk(void);
void determine_fields_at_point_CBDisk(double r, double theta, double phi,
   double *tdens, double *tpotl, double *tpotr, double *tpott, double *tpotp);
void pca_hall(int);


double *u, *du, **dend, **potd, **dpot, **work;
double *cosm, *sinm;

void get_acceleration_and_potential_CBDisk(void)
{
  static int firstime=1;
  int l, n;

  if (firstime) {
    
    dof = 2;			/* Two degrees of freedom */
				/* Three is default */

    /*	Allocate coefficient matrix                                     */
    /*    First row will always contain zeros for ease in interpolating */


    expcoef = dmatrix(0,2*lmax+1,1,nmax);
    expcoef1 = dmatrix(0,2*lmax+1,1,nmax);

    if (selector) {
      cc = (double ***)malloc((unsigned) (2*lmax+1)*sizeof(double **));
      if (!cc) {
	fprintf(stderr, "CBdiskacp: problem allocating <cc>\n");
	exit(-1);
      }
      for (l=0; l<=2*lmax; l++)
	cc[l] = dmatrix(1, nmax, 1, nmax);
  
      cc1 = (double ***)malloc((unsigned) (2*lmax+1)*sizeof(double **));
      if (!cc1) {
	fprintf(stderr, "CBdiskacp: problem allocating <cc1>\n");
	exit(-1);
      }
      for (l=0; l<=2*lmax; l++)
	cc1[l] = dmatrix(1, nmax, 1, nmax);
  
    }
  
    /* Allocate and compute normalization matrix */

    normM = dmatrix(0,lmax,1,nmax);
    dend = dmatrix(0,lmax,1,nmax);
    potd = dmatrix(0,lmax+1,1,nmax);
    dpot = dmatrix(0,lmax,1,nmax);
    work = dmatrix(0,lmax+1,1,nmax);

				/* Work vectors */
    u  = dvector(0,nmax);
    du = dvector(0,nmax);

    for (l=0; l<=lmax; l++) {
      for (n=1; n<=nmax; n++) {
	normM[l][n] = normCBDisk(n-1,l);
      }
    }
    
    /* Sin, cos */

    cosm = dvector(0,lmax);
    sinm = dvector(0,lmax);

  }


  /*======================*/
  /* Compute coefficients */
  /*======================*/

  if (firstime || self_consistent) {
    firstime = 0;
    determine_coefficients_CBDisk();
  }

  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  if (myid>0) determine_acceleration_and_potential_CBDisk();

  MPL_stop_timer();

}

void determine_coefficients_CBDisk(void)
{
  static int firstime=1, compute;
  int l, i, n, nn, use0, use1;
  double r, r2, rs, fac1, fac2, phi;

  if (selector) compute = !(this_step%npca) || firstime;

  /*		Clean */
  for (n=1; n<=nmax; n++) {
    for (l=0; l<=2*lmax; l++) {
      expcoef[l][n] = 0.0;
      expcoef1[l][n] = 0.0;
      if (selector && compute) {
	for (nn=n; nn<=nmax; nn++) cc1[l][n][nn] = 0.0;
      }
    }
  }

  use0 = 0;
  use1 = 0;

  if (myid>0) {

    /*		Begin by finding positions */
    for (i=1; i<=nbodies; i++) {

      if (freeze_particle(i)) continue;		/* frozen particles don't
						   contribute to field.
						   KL 5/27/92 */
      r2 = (x[i]*x[i] + y[i]*y[i]);
      r = rr[i] = sqrt(r2) + DSMALL;

      if (component[i] != 1) continue;

      if (r<=rmax) {
	use1++;
	phi = atan2(y[i],x[i]);
	rs = r/scale;
	
      
	sinecosine(lmax, phi, &cosm, &sinm);
	get_potl_CBDisk(lmax, nmax, rs, potd);

	/*		l loop */

	for (n=1; n<=nmax; n++) {
	  expcoef1[0][n] += potd[0][n]*mass[i]/normM[0][n];
	  if (selector && compute) {
	    for (nn=n; nn<=nmax; nn++)
	      cc1[0][n][nn] += potd[0][n]*potd[0][nn]*mass[i]/
		(normM[0][n]*normM[0][nn]);
	  }
	}
	
	for (l=1;l<=lmax; l++) {

	  fac1 = cosm[l];
	  fac2 = sinm[l];

	  for (n=1; n<=nmax; n++) {
	    expcoef1[2*l - 1][n] +=  potd[l][n]*fac1*mass[i]/normM[l][n];
	    expcoef1[2*l    ][n] +=  potd[l][n]*fac2*mass[i]/normM[l][n];
	    if (selector && compute) {
	      for (nn=n; nn<=nmax; nn++) {
		cc1[2*l - 1][n][nn] += potd[l][n]*potd[l][nn]*fac1*fac1*
		  mass[i]/(normM[l][n]*normM[l][nn]);
		cc1[2*l    ][n][nn] += potd[l][n]*potd[l][nn]*fac2*fac2*
		  mass[i]/(normM[l][n]*normM[l][nn]);
	      }
	    }
	  }
	}
      }
    }
  }


  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (myid==0) used += use0;


  if (!selector) {

    MPI_Allreduce ( &expcoef1[0][1],
		    &expcoef[0][1],
		    nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	
    for (l=1;l<=lmax; l++) {

      MPI_Allreduce ( &expcoef1[2*l - 1][1],
		      &expcoef[2*l - 1][1],
		      nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce ( &expcoef1[2*l    ][1],
		      &expcoef[2*l    ][1],
		      nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
  }

  if (selector) {

    parallel_gather_coefficients();

    if (myid == 0) pca_hall(compute);

    parallel_distribute_coefficients();

    firstime = 0;
  }

}

void determine_acceleration_and_potential_CBDisk(void)
{
  int l,i;
  double r,rs,fac,phi,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps;

  /* Determine potential and acceleration */

  for (i=1; i<=nbodies; i++) {

    if (freeze_particle(i)) 	/* frozen particles do not respond */
      continue;			/* KL 5/27/92 */

    if (!self_consistent) {
      r = rr[i] = sqrt(x[i]*x[i] + y[i]*y[i]) + DSMALL;
    }
    else
      r = rr[i];

    rs = r/scale;
    phi = atan2(y[i],x[i]);

    get_dpotl_CBDisk(lmax, nmax, rs, potd, dpot);


    get_pot_coefs_CBDisk(0,expcoef[0],&p,&dp);
    potl = p;
    potr = dp;
    pott = potp = 0.0;
      
    /*		l loop */
    
    for (l=1; l<=lmax; l++) {

      get_pot_coefs_CBDisk(l,expcoef[2*l - 1],&pc,&dpc);
      get_pot_coefs_CBDisk(l,expcoef[2*l    ],&ps,&dps);
      potl += pc*cosm[l] + ps*sinm[l];
      potr += dpc*cosm[l] + dps*sinm[l];
      potp += (-pc*sinm[l] + ps*cosm[l])*l;
    }

    fac = x[i]*x[i] + y[i]*y[i];
    
    potr /= scale*scale;
    potl /= scale;
    potp /= scale;

    ax[i] += -potr*x[i]/r;
    ay[i] += -potr*y[i]/r;
    az[i] += -potr*z[i]/r;
    if (fac > DSMALL2) {
      ax[i] +=  potp*y[i]/fac;
      ay[i] += -potp*x[i]/fac;
    }
    pot[i] += potl;

  }


}


void determine_fields_at_point_CBDisk(double r, double theta, double phi,
   double *tdens, double *tpotl, double *tpotr, double *tpott, double *tpotp)
{
  int l;
  double rs,costh,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps,dens;

  rs = r/scale;
  costh = cos(theta);

  get_dens_CBDisk(lmax, nmax, rs, dend);
  get_dpotl_CBDisk(lmax, nmax, rs, potd, dpot);

  get_dens_coefs_CBDisk(0,expcoef[0],&dens);

  get_pot_coefs_CBDisk(0,expcoef[0],&p,&dp);
  potl = p;
  potr = dp;
  pott = potp = 0.0;
  
      
  /*		l loop */
    
  for (l=1; l<=lmax; l++) {
    
    get_dens_coefs_CBDisk(l,expcoef[2*l - 1],&pc);
    get_dens_coefs_CBDisk(l,expcoef[2*l    ],&ps);
    dens += pc*cosm[l] + ps*sinm[l];
    
    get_pot_coefs_CBDisk(l,expcoef[2*l - 1],&pc,&dpc);
    get_pot_coefs_CBDisk(l,expcoef[2*l    ],&ps,&dps);
    potl += pc*cosm[l] + ps*sinm[l];
    potr += dpc*cosm[l] + dps*sinm[l];
    potp += (-pc*sinm[l] + ps*cosm[l])*l;
  }

  *tdens = dens/(scale*scale*scale);
  *tpotl = potl/scale;
  *tpotr = potr/(scale*scale);
  *tpotp = potp/scale;
  
}


void get_pot_coefs_CBDisk(int l, double *coef, double *p, double *dp)
{
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (i=1; i<=nmax; i++) {
    pp  += potd[l][i] * coef[i];
    dpp += dpot[l][i] * coef[i];
  }

  *p = -pp;
  *dp = -dpp;
}


void get_dens_coefs_CBDisk(int l, double *coef, double *p)
{
  double pp;
  int i;

  pp = 0.0;

  for (i=1; i<=nmax; i++)
    pp  += dend[l][i] * coef[i];

  *p = pp;
}


				/* Get potential functions by recursion */

void get_dpotl_CBDisk(int lmax, int nmax, double r, double **p, double **dp)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac), rcum = 1.0;
  int l, nn;

  for (l=0; l<=lmax+1; l++) {
    cur = cur0;

    p[l][1] = cur*rcum;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      p[l][nn+1] = cur*rcum;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcum *= r;
  }


  for (l=0; l<=lmax; l++) {

    dp[l][1] = p[l][1]*l/r - p[l+1][1];
    if (nmax<1) break;
    dp[l][2] = p[l][2]*l/r - (p[l+1][2] - 2*p[l+1][1]);
    if (nmax<2) break;

    for (nn=2; nn<nmax; nn++) {
      dp[l][nn+1] = p[l][nn+1]*l/r - 
	( p[l+1][nn+1] - 2*p[l+1][nn] + p[l+1][nn-1] );
    }
  }

}

void get_potl_CBDisk(int lmax, int nmax, double r, double **p)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac), rcum = 1.0;
  int l, nn;

  for (l=0; l<=lmax+1; l++) {
    cur = cur0;

    work[l][1]  = cur;
    p[l][1] = cur*rcum;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      p[l][nn+1] = cur*rcum;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcum *= r;
  }

}


void get_dens_CBDisk(int lmax, int nmax, double r, double **d)
{

  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac);
  int l, nn;
  double rcum = 0.5/M_PI;

  for (l=0; l<=lmax+1; l++) {
    cur = cur0;

    work[l][1] = cur;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      work[l][nn+1] = cur;
    }
    cur0 *= fac*(2*(l+1) - 1);
  }

  for (l=0; l<=lmax; l++) {
    d[l][1] = work[l+1][1]*rcum;
    d[l][2] = work[l+1][2]*rcum;
    for (nn=2; nn<nmax; nn++)
      d[l][nn+1] = (work[l+1][nn+1] - work[l+1][nn-1])*rcum;

    rcum *= r;
  }

}

void get_potl_dens_CBDisk(int lmax, int nmax, double r, double **p, double **d)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac);
  int l, nn;
  double  rcump = 1.0, rcumd = 0.5/M_PI;

  for (l=0; l<=lmax+1; l++) {
    cur = cur0;

    work[l][1] = cur;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      p[l][nn+1] = cur*rcump;
      work[l][nn+1] = cur;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcump *= r;
  }

  for (l=0; l<=lmax; l++) {
    d[l][1] = work[l+1][1]*rcumd;
    d[l][2] = work[l+1][2]*rcumd;
    for (nn=2; nn<nmax; nn++)
      d[l][nn+1] = (work[l+1][nn+1] - work[l+1][nn-1])*rcumd;

    rcumd *= r;
  }

}

double normCBDisk(int n, int m)
{
  double ans = 1.0;
  int i;
 
  for (i=n+1; i<=n+2*m; i++)
    ans *= i;

  return pow(0.5, 2*m+1)*ans;
}



				/* Dump coefficients to a file */

void dump_coefs_CBDisk(FILE *fout)
{
  int ir, l;

  fwrite(&tnow, sizeof(double), 1, fout);
  fwrite(&scale, sizeof(double), 1, fout);
  fwrite(&nmax, sizeof(int), 1, fout);
  fwrite(&lmax, sizeof(int), 1, fout);

  for (ir=1; ir<=nmax; ir++) {
    for (l=0; l<=2*lmax; l++)
      fwrite(&expcoef[l][ir], sizeof(double), 1, fout);
  }

}

