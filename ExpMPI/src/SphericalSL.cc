/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  the Sturm Liouville direct solution
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

#include <stdlib.h>
#include <iostream>
#include <stdexcept>

#include "expand.h"

#include <SphericalSL.h>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void legendre_R(int lmax, double x, Matrix& p);
void dlegendre_R(int lmax, double x, Matrix &p, Matrix &dp);
void sinecosine_R(int mmax, double phi, Vector& c, Vector& s);
extern "C" {
  void get_potl_dens_SLsph(int l, int n, double r);
  void get_dens_coefs_SLsph(int l, Vector& coef, double *p);
  void get_pot_coefs_SLsph(int l, Vector& coef, double *p, double *dp);
}

void pca_hall(int compute, Matrix& expcoef, Matrix*& cc, Matrix& normM);
void parallel_gather_coefficients(Matrix& expcoef, Matrix& expcoef1,
				  Matrix*& cc, Matrix*& cc1,
				  int lmax);
void parallel_distribute_coefficients(Matrix& expcoef, int lmax);


double SphericalSL::RMIN = 0.001;
double SphericalSL::RMAX = 100.0;
int SphericalSL::NUMR = 1000;

SphericalSL::SphericalSL(void)
{
}

SphericalSL::SphericalSL(int lmax, int nmax)
{
  reset(lmax, nmax);
}

void SphericalSL::reset(int lmax, int nmax)
{
  NMAX = nmax;
  LMAX = lmax<1 ? 1 : lmax;
  
  if (nthrds<1) nthrds=1;
  
				//  Allocate coefficient matrix
  expcoef  = Matrix(0, LMAX*(LMAX+2), 1, NMAX);
  expcoef1 = Matrix(0, LMAX*(LMAX+2), 1, NMAX);
  
  try {
    expcoef0  = new Matrix [nthrds];
  }
  catch(exception const & msg)
    {
      cerr << "SphericalSL: " << msg.what() << endl;
      return;
    }
  
  for (int j=0; j<nthrds; j++)
    expcoef0[j] = Matrix(0, LMAX*(LMAX+2), 1, NMAX);
  
  if (selector) {
#ifdef DEBUG
    cerr << "Process " << myid << ": Selector is 1\n";
#endif

    try {
      cc = new Matrix [LMAX*(LMAX+2) + 1];
      cc1 = new Matrix [LMAX*(LMAX+2) + 1];
    }
    catch(exception const & msg)
      {
	cerr << "SphericalSL: " << msg.what() << endl;
	return;
      }
    
    for (int l=0; l<=LMAX*(LMAX+2); l++)
      cc[l] = Matrix(1, NMAX, 1, NMAX);
    
    for (int l=0; l<=LMAX*(LMAX+2); l++)
      cc1[l] = Matrix(1, NMAX, 1, NMAX);
  }
  
  // Allocate and compute normalization matrix
  
  normM = Matrix(0, LMAX, 1, NMAX);
  krnl  = Matrix(0, LMAX, 1, NMAX);
  dend  = Matrix(0, LMAX, 1, NMAX);
  
  try {
    potd  = new Matrix [nthrds];
    dpot  = new Matrix [nthrds];
    legs  = new Matrix [nthrds];
    dlegs = new Matrix [nthrds];
    cosm  = new Vector [nthrds];
    sinm  = new Vector [nthrds];
  }
  catch(exception const & msg)
    {
      cerr << "SphericalSL: " << msg.what() << endl;
      return;
    }
  
				// Potential
  for (int i=0; i<nthrds; i++) {
    potd[i]  = Matrix(0, LMAX, 1, NMAX);
    dpot[i]  = Matrix(0, LMAX, 1, NMAX);
  }
  
				// Sin, cos, legendre
  
  for (int i=0; i<nthrds; i++) {
    cosm[i] = Vector(0, LMAX);
    sinm[i] = Vector(0, LMAX);
    legs[i] = Matrix(0, LMAX, 0, LMAX);
    dlegs[i] = Matrix(0, LMAX, 0, LMAX);
  }
  
  for (int l=0; l<=LMAX; l++) {
    for (int n=1; n<=NMAX; n++) {
      normM[l][n] = 1.0;
      krnl[l][n] = 1.0;
    }
  }
  
				// Factorial matrix
  
  factorial = Matrix(0, LMAX, 0, LMAX);
  
  for (int l=0; l<=LMAX; l++) {
    for (int m=0; m<=l; m++) 
      factorial[l][m] = factrl(l-m)/factrl(l+m);
  }
  
  use = new int [nthrds];

				// Generate Sturm-Liouville grid
  SLGridSph::mpi = 1;		// Turn on MPI
  ortho = new SLGridSph(LMAX, NMAX, NUMR, RMIN, RMAX);

  if (selector) pthread_mutex_init(&cc_lock, NULL);

}


SphericalSL::~SphericalSL(void)
{
  delete [] expcoef0;
  if (selector) {
    delete [] cc;
    delete [] cc1;
  }
  delete [] potd;
  delete [] dpot;
  delete [] cosm;
  delete [] sinm;
  delete [] legs;
  delete [] dlegs;
  delete [] use;

  delete ortho;
  if (selector) pthread_mutex_destroy(&cc_lock);

}

void SphericalSL::get_acceleration_and_potential_SLsph(void)
{
  static bool firstime=true;
  int l, m, n, i;


  /*======================*/
  /* Compute coefficients */
  /*======================*/

  if (firstime || self_consistent) {
    firstime = false;
    determine_coefficients_SLsph();
  }

  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  cout << "Process " << myid << ": about to call determine_acceration\n"
       << flush;
  determine_acceleration_and_potential_SLsph();

}

void SphericalSL::thread_call(void *arg)
{
  if (do_accel) {
    determine_acceleration_and_potential_SLsph_thread(arg);
  }
  else {
    determine_coefficients_SLsph_thread(arg);
  }
}

void SphericalSL::determine_coefficients_SLsph_thread(void * arg)
{
  int l, i, loffset, moffset, m, n, nn;
  double r, r2, fac1, fac2, costh, phi;
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

    if (r<=rmax && r<=RMAX) {
      use[id]++;
      costh = zz/r;
      phi = atan2(yy,xx);
      
      legendre_R(LMAX, costh, legs[id]);
      sinecosine_R(LMAX, phi, cosm[id], sinm[id]);

      // get_potl_SLsph_safe(LMAX, NMAX, r, potd[id], u[id]);
      
      ortho->get_pot(potd[id], r);

      /*		l loop */
      for (l=0, loffset=0; l<=LMAX; loffset+=(2*l+1), l++) {
	/*		m loop */
	for (m=0, moffset=0; m<=l; m++) {
	  if (m==0) {
	    if (selector && compute)
	      facs1 = legs[id][l][m]*legs[id][l][m]*mass[i];
	    for (n=1; n<=NMAX; n++) {
	      expcoef0[id][loffset+moffset][n] += potd[id][l][n]*legs[id][l][m]*mass[i]*
		fac0/normM[l][n];

	      if (selector && compute) {
		pthread_mutex_lock(&cc_lock);
		for (nn=n; nn<=NMAX; nn++)
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
	    for (n=1; n<=NMAX; n++) {
	      expcoef0[id][loffset+moffset][n] += potd[id][l][n]*fac1*mass[i]*
		fac0/normM[l][n];

	      expcoef0[id][loffset+moffset+1][n] += potd[id][l][n]*fac2*mass[i]*
		fac0/normM[l][n];

	      if (selector && compute) {
		pthread_mutex_lock(&cc_lock);
		for (nn=n; nn<=NMAX; nn++) {
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

}

void SphericalSL::determine_coefficients_SLsph(void)
{
  static int firstime=1;
  int l, i, loffset, moffset, m, n, nn, use0, use1;

  if (selector) compute = !(this_step%npca) || firstime;

  /*		Clean */
  for (n=1; n<=NMAX; n++) {
      for (l=0; l<=LMAX*(LMAX+2); l++) {
	expcoef[l][n] = 0.0;
	expcoef1[l][n] = 0.0;
	for (i=0; i<nthrds; i++) expcoef0[i][l][n] = 0.0;
	if (selector && compute) 
	  for (nn=n; nn<=NMAX; nn++) cc1[l][n][nn] = 0.0;
      }
    }

  use0 = 0;
  use1 = 0;
  do_accel = false;

  start_threads(nthrds);

  for (i=0; i<nthrds; i++) {
    use1 += use[i];
    for (l=0; l<= LMAX*(LMAX+2); l++)
      for (n=1; n<=NMAX; n++)
	expcoef1[l][n] += expcoef0[i][l][n];
  }
    

  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (myid==0) used += use0;


  if (!selector) {

    for (l=0, loffset=0; l<=LMAX; loffset+=(2*l+1), l++) {
	/*		m loop */
	for (m=0, moffset=0; m<=l; m++) {
	  if (m==0) {

	    MPI_Allreduce ( &expcoef1[loffset+moffset][1],
			   &expcoef[loffset+moffset][1],
			   NMAX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	    moffset++;
	  }
	  else {

	    MPI_Allreduce ( &expcoef1[loffset+moffset][1],
			   &expcoef[loffset+moffset][1],
			   NMAX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	    MPI_Allreduce ( &expcoef1[loffset+moffset+1][1],
			   &expcoef[loffset+moffset+1][1],
			   NMAX, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    moffset+=2;
	  }
	}
      }
  }
 
  if (selector) {


    parallel_gather_coefficients(expcoef, expcoef1, cc, cc1, lmax);

    if (myid == 0) pca_hall(compute, expcoef, cc, normM);

    parallel_distribute_coefficients(expcoef, lmax);

    firstime = 0;
  }

#ifdef DEBUG
				// Check coefficients
  int iflg = 0;

  for (n=1; n<=NMAX; n++) {
    for (l=0; l<=LMAX*(LMAX+2); l++) {
      if (isnan(expcoef[l][n])) {
	cerr.form("expcoef[%d][%d] is NaN\n", l, n);
	iflg++;
      }
    }
  }
  if (iflg) {
    cerr << iflg << " NaNs\n";
    MPI_Finalize();
    exit(-11);
  }
#endif
}

void SphericalSL::determine_acceleration_and_potential_SLsph_thread(void * arg)
{
  int l,i,loffset,moffset,m;
  double r,fac,fac1,fac2,fac3,fac4,costh,phi,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps;
  double dfac=0.25/M_PI, r1, pfext1, pfext2;
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
    phi = atan2(yy,xx);

    dlegendre_R(LMAX, costh, legs[id], dlegs[id]);
    sinecosine_R(LMAX, phi, cosm[id], sinm[id]);
    // get_dpotl_SLsph_safe(LMAX, NMAX, r, potd[id], dpot[id], u[id], du[id]);

				// For exterior solution
    pfext1 = 1.0;
    pfext2 = 1.0;
    r1 = r;
    if (r>RMAX) {
      pfext1 = RMAX/r;
      pfext2 = pfext1;
      r1 = RMAX;
    }
      
    ortho->get_pot(potd[id], r1);
    ortho->get_force(dpot[id], r1);

    get_pot_coefs_SLsph_safe(0, expcoef[0], &p, &dp, potd[id], dpot[id]);
    potl = fac1*p * pfext2;
    potr = fac1*dp * pfext2*pfext1;
    pott = potp = 0.0;
      
    /*		l loop */
    
    for (l=1, loffset=1; l<=LMAX; loffset+=(2*l+1), l++) {

      pfext2 *= pfext1;

      /*		m loop */
      for (m=0, moffset=0; m<=l; m++) {
	fac1 = (2.0*l+1.0)/(4.0*M_PI);
	if (m==0) {
				/* Suppress L=1 terms? */
	  if ( !NO_L1 || !(l==1) ) {
	    fac2 = fac1*legs[id][l][m];
	    get_pot_coefs_SLsph_safe(l, expcoef[loffset+moffset], &p, &dp,
				  potd[id], dpot[id]);

				// Exterior solution
	    p *= pfext2;
	    dp *= pfext2*pfext1;

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
	    get_pot_coefs_SLsph_safe(l,expcoef[loffset+moffset],&pc,&dpc,
				  potd[id],dpot[id]);
	    get_pot_coefs_SLsph_safe(l,expcoef[loffset+moffset+1],&ps,&dps,
				  potd[id],dpot[id]);

				// Exterior solution
	    pc *= pfext2;
	    dpc *= pfext2*pfext1;
	    ps *= pfext2;
	    dps *= pfext2*pfext1;

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

    ax[i] += -(potr*xx/r - pott*xx*zz/(r*r*r) );
    ay[i] += -(potr*yy/r - pott*yy*zz/(r*r*r) );
    az[i] += -(potr*zz/r + pott*fac/(r*r*r));
    if (fac > DSMALL2) {
      ax[i] +=  potp*yy/fac;
      ay[i] += -potp*xx/fac;
    }
    pot[i] += potl;

  }

  return;
}


void SphericalSL::determine_acceleration_and_potential_SLsph(void)
{
#ifdef MPE_PROFILE
  MPE_Log_event(11, myid, "b_compute_force");
#endif

  do_accel = true;
  start_threads(nthrds);

#ifdef MPE_PROFILE
  MPE_Log_event(12, myid, "e_compute_force");
#endif
}


void SphericalSL::determine_fields_at_point_SLsph(
   double r, double theta, double phi,
   double *tdens, double *tpotl, double *tpotr, double *tpott, double *tpotp)
{
  int l,loffset,moffset,m;
  double fac1,fac2,fac3,fac4,costh,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps,dens;
  double dfac=0.25/M_PI, pfext1, pfext2, r1;


  costh = cos(theta);

  fac1 = dfac;

  dlegendre_R(LMAX, costh, legs[0], dlegs[0]);
  sinecosine_R(LMAX, phi, cosm[0], sinm[0]);
  // get_dens_SLsph(LMAX, NMAX, r, dend);
  // get_dpotl_SLsph(LMAX, NMAX, r, potd[0], dpot[0]);

				// For exterior solution
  pfext1 = 1.0;
  pfext2 = 1.0;
  r1 = r;
  if (r>RMAX) {
    pfext1 = RMAX/r;
    pfext2 = pfext1;
    r1 = RMAX;
  }

  ortho->get_dens(dend, r);
  ortho->get_pot(potd[0], r1);
  ortho->get_force(dpot[0], r1);

  get_dens_coefs_SLsph(0,expcoef[0],&dens);
  dens *= dfac*dfac;

  get_pot_coefs_SLsph(0,expcoef[0],&p,&dp);
  potl = fac1*p * pfext2;
  potr = fac1*dp * pfext2*pfext1;
  pott = potp = 0.0;
  
      
  /*		l loop */
    
  for (l=1, loffset=1; l<=LMAX; loffset+=(2*l+1), l++) {
    
    /*		m loop */
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
	fac2 = fac1*legs[0][l][m];
	get_dens_coefs_SLsph(l,expcoef[loffset+moffset],&p);
	dens += dfac*fac2*p;
	get_pot_coefs_SLsph(l,expcoef[loffset+moffset],&p,&dp);

				// External solution
	p *= pfext2;
	dp *= pfext2*pfext1;

	potl += fac2*p;
	potr += fac2*dp;
	pott += fac1*dlegs[0][l][m]*p;
	moffset++;
      }
      else {
	fac2 = 2.0 * fac1 * factorial[l][m];
	fac3 = fac2 * legs[0][l][m];
	fac4 = fac2 * dlegs[0][l][m];
	
	get_dens_coefs_SLsph(l,expcoef[loffset+moffset],&pc);
	get_dens_coefs_SLsph(l,expcoef[loffset+moffset+1],&ps);
	dens += dfac*fac3*(pc*cosm[0][m] + ps*sinm[0][m]);
	
	get_pot_coefs_SLsph(l,expcoef[loffset+moffset],&pc,&dpc);
	get_pot_coefs_SLsph(l,expcoef[loffset+moffset+1],&ps,&dps);

				// External solution
	pc *= pfext2;
	dpc *= pfext2*pfext1;
	ps *= pfext2;
	dps *= pfext2*pfext1;

	potl += fac3*(pc*cosm[0][m] + ps*sinm[0][m]);
	potr += fac3*(dpc*cosm[0][m] + dps*sinm[0][m]);
	pott += fac4*(pc*cosm[0][m] + ps*sinm[0][m]);
	potp += fac3*(-pc*sinm[0][m] + ps*cosm[0][m])*m;
	moffset +=2;
      }
    }
  }

  *tdens = dens;
  *tpotl = potl;
  *tpotr = potr;
  *tpott = pott;
  *tpotp = potp;
  
}


void SphericalSL::get_pot_coefs_SLsph(int l, Vector& coef, 
				      double *p, double *dp)
{
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (i=1; i<=NMAX; i++) {
    pp  += potd[0][l][i] * coef[i];
    dpp += dpot[0][l][i] * coef[i];
  }

  *p = -pp;
  *dp = -dpp;
}


void SphericalSL::get_pot_coefs_SLsph_safe(int l, Vector& coef, 
					   double *p, double *dp,
					   Matrix& potd1, Matrix& dpot1)
{
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (i=1; i<=NMAX; i++) {
    pp  += potd1[l][i] * coef[i];
    dpp += dpot1[l][i] * coef[i];
  }

  *p = -pp;
  *dp = -dpp;
}


void SphericalSL::get_dens_coefs_SLsph(int l, Vector& coef, double *p)
{
  double pp;
  int i;

  pp = 0.0;

  for (i=1; i<=NMAX; i++)
    pp  += dend[l][i] * coef[i];

  *p = pp;
}
				/* Dump coefficients to a file */

void SphericalSL::dump_coefs_SLsph(FILE *fout)
{
  int ir, l;

  fwrite(&tnow, sizeof(double), 1, fout);

  for (ir=1; ir<=NMAX; ir++) {
    for (l=0; l<=LMAX*(LMAX+2); l++)
      fwrite(&expcoef[l][ir], sizeof(double), 1, fout);
  }

}

