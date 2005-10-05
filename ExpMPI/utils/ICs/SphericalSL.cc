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
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *      06/09/92 updated to use recursion relations rather than tables
 *
 ***************************************************************************/

#include <stdlib.h>
#include <values.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include <SphericalSL.h>

#ifdef RCSID
static char rcsid[] = 
"$Id$";
#endif

void legendre_R(int lmax, double x, Matrix& p);
void dlegendre_R(int lmax, double x, Matrix &p, Matrix &dp);
void sinecosine_R(int mmax, double phi, Vector& c, Vector& s);
double factrl(int n);

double SphericalSL::RMIN = 0.001;
double SphericalSL::RMAX = 100.0;
int SphericalSL::NUMR = 1000;

int SphericalSL::selector = 0;
int SphericalSL::tk_type = 2;
double SphericalSL::tksmooth = 1.0;
double SphericalSL::tkcum = 0.95;

SphericalSL::SphericalSL(void)
{
  compute = 0;
}

SphericalSL::SphericalSL(int lmax, int nmax, double SCALE)
{
  reset(lmax, nmax, SCALE);
  compute = 0;
}

void SphericalSL::reset(int lmax, int nmax, double SCALE)
{
  NMAX = nmax;
  LMAX = lmax<1 ? 1 : lmax;
  
				//  Allocate coefficient matrix
  expcoef  = Matrix(0, LMAX*(LMAX+2), 1, NMAX);
  expcoef1 = Matrix(0, LMAX*(LMAX+2), 1, NMAX);
  
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
  
				// Potential
  potd = Matrix(0, LMAX, 1, NMAX);
  dpot = Matrix(0, LMAX, 1, NMAX);
  
				// Sin, cos, legendre
  
  cosm = Vector(0, LMAX);
  sinm = Vector(0, LMAX);
  legs = Matrix(0, LMAX, 0, LMAX);
  dlegs = Matrix(0, LMAX, 0, LMAX);
  
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
  
				// Generate Sturm-Liouville grid
  SLGridSph::mpi = 1;		// Turn on MPI
  ortho = new SLGridSph(LMAX, NMAX, NUMR, RMIN, RMAX, 1, SCALE);

}


SphericalSL::~SphericalSL(void)
{
  if (selector) {
    delete [] cc;
    delete [] cc1;
  }

  delete ortho;
}


void SphericalSL::compute_coefficients(vector<Particle> &part)
{
  int l, loffset, moffset, m, n, nn;
  double r, r2, fac1, fac2, costh, phi;
  double facs1=0.0, facs2=0.0, fac0=4.0*M_PI;
  double xx, yy, zz, mass;

  use = 0;

  vector<Particle>::iterator p;
  for (p=part.begin(); p!=part.end(); p++) {
    
    xx = p->pos[0];
    yy = p->pos[1];
    zz = p->pos[2];
    mass = p->mass;

    r2 = (xx*xx + yy*yy + zz*zz);
    r = sqrt(r2) + MINDOUBLE;

    if (r<=RMAX) {
      use++;
      costh = zz/r;
      phi = atan2(yy,xx);
      
      legendre_R(LMAX, costh, legs);
      sinecosine_R(LMAX, phi, cosm, sinm);

      // get_potl_safe(LMAX, NMAX, r, potd, u);
      
      ortho->get_pot(potd, r);

      /*		l loop */
      for (l=0, loffset=0; l<=LMAX; loffset+=(2*l+1), l++) {
	/*		m loop */
	for (m=0, moffset=0; m<=l; m++) {
	  if (m==0) {
	    if (selector && compute)
	      facs1 = legs[l][m]*legs[l][m]*mass;
	    for (n=1; n<=NMAX; n++) {
	      expcoef1[loffset+moffset][n] += potd[l][n]*legs[l][m]*mass*
		fac0/normM[l][n];

	      if (selector && compute) {
		for (nn=n; nn<=NMAX; nn++)
		  cc1[loffset+moffset][n][nn] += potd[l][n]*potd[l][nn]*
		    facs1/(normM[l][n]*normM[l][nn]);
	      }
	    }
	    moffset++;
	  }
	  else {
	    fac1 = legs[l][m]*cosm[m];
	    fac2 = legs[l][m]*sinm[m];
	    if (selector && compute) {
	      facs1 = fac1*fac1*mass;
	      facs2 = fac2*fac2*mass;
	    }
	    for (n=1; n<=NMAX; n++) {
	      expcoef1[loffset+moffset][n] += potd[l][n]*fac1*mass*
		fac0/normM[l][n];

	      expcoef1[loffset+moffset+1][n] += potd[l][n]*fac2*mass*
		fac0/normM[l][n];

	      if (selector && compute) {
		for (nn=n; nn<=NMAX; nn++) {
		  cc1[loffset+moffset][n][nn] += 
		    potd[l][n]*potd[l][nn]*facs1/(normM[l][n]*normM[l][nn]);
		  cc1[loffset+moffset+1][n][nn] +=
		    potd[l][n]*potd[l][nn]*facs2/(normM[l][n]*normM[l][nn]);
		}
	      }
		
	    }
	    moffset+=2;
	  }
	}
      }
    }
  }

				// Check coefficients
  int iflg = 0;

  for (n=1; n<=NMAX; n++) {
    for (l=0; l<=LMAX*(LMAX+2); l++) {
      if (isnan(expcoef[l][n])) {
	cerr << "expcoef[" << l << "][" << n << "] is NaN" << endl;
	iflg++;
      }
    }
  }
  if (iflg) {
    cerr << iflg << " NaNs\n";
    MPI_Finalize();
    exit(-11);
  }
}


void SphericalSL::accumulate(vector<Particle> &part)
{
  static int firstime=1;
  int l, loffset, moffset, m, n, nn, use0, use1;
  
  used = 0;

  if (selector) compute = firstime;

  /*		Clean */
  for (n=1; n<=NMAX; n++) {
      for (l=0; l<=LMAX*(LMAX+2); l++) {
	expcoef[l][n] = 0.0;
	expcoef1[l][n] = 0.0;
	if (selector && compute) 
	  for (nn=n; nn<=NMAX; nn++) cc1[l][n][nn] = 0.0;
      }
    }

  use0 = 0;
  use1 = 0;

  compute_coefficients(part);

  use1 += use;
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

    parallel_gather_coefficients(expcoef, expcoef1, cc, cc1, LMAX);

    if (myid == 0) pca_hall(compute, expcoef, cc, normM);

    parallel_distribute_coefficients(expcoef, LMAX);

    if (firstime) {
      compute = 0;
      firstime = 0;
    }
  }

#ifdef DEBUG
				// Check coefficients
  int iflg = 0;

  for (n=1; n<=NMAX; n++) {
    for (l=0; l<=LMAX*(LMAX+2); l++) {
      if (isnan(expcoef[l][n])) {
	cerr << "expcoef[" << l << "][" << n << "] is NaN" << endl;
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



void SphericalSL::determine_fields_at_point
(
 double r, double theta, double phi,
 double *tdens, double *tpotl, double *tpotr, double *tpott, double *tpotp
 )
{
  int l,loffset,moffset,m;
  double fac1,fac2,fac3,fac4,costh,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps,dens;
  double dfac=0.25/M_PI, pfext1, pfext2, r1;


  costh = cos(theta);

  fac1 = dfac;

  dlegendre_R(LMAX, costh, legs, dlegs);
  sinecosine_R(LMAX, phi, cosm, sinm);

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
  ortho->get_pot(potd, r1);
  ortho->get_force(dpot, r1);

  get_dens_coefs(0,expcoef[0],&dens);
  dens *= dfac*dfac;

  get_pot_coefs(0,expcoef[0],&p,&dp);
  potl = fac1*p * pfext2;
  potr = fac1*dp * pfext2*pfext1;
  pott = potp = 0.0;
  
      
  /*		l loop */
    
  for (l=1, loffset=1; l<=LMAX; loffset+=(2*l+1), l++) {
    
    /*		m loop */
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
	fac2 = fac1*legs[l][m];
	get_dens_coefs(l,expcoef[loffset+moffset],&p);
	dens += dfac*fac2*p;
	get_pot_coefs(l,expcoef[loffset+moffset],&p,&dp);

				// External solution
	p *= pfext2;
	dp *= pfext2*pfext1;

	potl += fac2*p;
	potr += fac2*dp;
	pott += fac1*dlegs[l][m]*p;
	moffset++;
      }
      else {
	fac2 = 2.0 * fac1 * factorial[l][m];
	fac3 = fac2 * legs[l][m];
	fac4 = fac2 * dlegs[l][m];
	
	get_dens_coefs(l,expcoef[loffset+moffset],&pc);
	get_dens_coefs(l,expcoef[loffset+moffset+1],&ps);
	dens += dfac*fac3*(pc*cosm[m] + ps*sinm[m]);
	
	get_pot_coefs(l,expcoef[loffset+moffset],&pc,&dpc);
	get_pot_coefs(l,expcoef[loffset+moffset+1],&ps,&dps);

				// External solution
	pc *= pfext2;
	dpc *= pfext2*pfext1;
	ps *= pfext2;
	dps *= pfext2*pfext1;

	potl += fac3*(pc*cosm[m] + ps*sinm[m]);
	potr += fac3*(dpc*cosm[m] + dps*sinm[m]);
	pott += fac4*(pc*cosm[m] + ps*sinm[m]);
	potp += fac3*(-pc*sinm[m] + ps*cosm[m])*m;
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


void SphericalSL::get_pot_coefs(int l, Vector& coef, 
				      double *p, double *dp)
{
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (i=1; i<=NMAX; i++) {
    pp  += potd[l][i] * coef[i];
    dpp += dpot[l][i] * coef[i];
  }

  *p = -pp;
  *dp = -dpp;
}


void SphericalSL::get_pot_coefs_safe(int l, Vector& coef, 
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


void SphericalSL::get_dens_coefs(int l, Vector& coef, double *p)
{
  double pp;
  int i;

  pp = 0.0;

  for (i=1; i<=NMAX; i++)
    pp  += dend[l][i] * coef[i];

  *p = pp;
}
				/* Dump coefficients to a file */

void SphericalSL::dump_coefs(ofstream& out)
{
  double tnow = 0.0;
  int ir, l;

  out.write((char *)&tnow, sizeof(double));

  for (ir=1; ir<=NMAX; ir++) {
    for (l=0; l<=LMAX*(LMAX+2); l++)
      out.write((char *)&expcoef[l][ir], sizeof(double));
  }

}

void SphericalSL::parallel_gather_coefficients
(
 Matrix& expcoef, Matrix& expcoef1,
 Matrix*& cc, Matrix*& cc1,
 int lmax)
{
  int Ldim, L0, loffset, moffset, l, m, n, nn;

  Ldim = lmax*(lmax + 2) + 1;
  L0 = 0;
  
  if (myid == 0) {

    for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

      for (m=0, moffset=0; m<=l; m++) {

	if (m==0) {
	  for (n=1; n<=NMAX; ++n) {
	    expcoef[loffset+moffset][n] = 0.0;

	    for (nn=n; nn<=NMAX; nn++)
	      cc[loffset+moffset][n][nn] = 0.0;
	  }
	  moffset++;
	}
	else {
	  for (n=1; n<=NMAX; ++n) {
	    expcoef[loffset+moffset][n] = 0.0;
	    expcoef[loffset+moffset+1][n] = 0.0;

	    for (nn=n; nn<=NMAX; nn++) {
	      cc[loffset+moffset][n][nn] = 0.0;
	      cc[loffset+moffset+1][n][nn] = 0.0;
	    }
	  }
	  moffset+=2;
	}
      }
    }
  }


  for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

    for (m=0, moffset=0; m<=l; m++) {

      if (m==0) {
	MPI_Reduce(&expcoef1[loffset+moffset][1], 
		   &expcoef[loffset+moffset][1], NMAX, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	for (n=1; n<=NMAX; n++)
	  MPI_Reduce(&cc1[loffset+moffset][n][n],
		     &cc[loffset+moffset][n][n], NMAX-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	moffset++;
      }
      else {
	MPI_Reduce(&expcoef1[loffset+moffset][1], 
		   &expcoef[loffset+moffset][1], NMAX, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&expcoef1[loffset+moffset+1][1],
		   &expcoef[loffset+moffset+1][1], NMAX, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	for (n=1; n<=NMAX; n++) {
	  MPI_Reduce(&cc1[loffset+moffset][n][n],
		     &cc[loffset+moffset][n][n], NMAX-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&cc1[loffset+moffset+1][n][n],
		     &cc[loffset+moffset+1][n][n], NMAX-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	moffset+=2;
      }
    }
  }

}

void SphericalSL::parallel_distribute_coefficients(Matrix& expcoef, int lmax)
{
  int Ldim, L0, loffset, moffset, l, m;

  Ldim = lmax*(lmax + 2) + 1;
  L0 = 0;

  for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

      for (m=0, moffset=0; m<=l; m++) {

	if (m==0) {
	  MPI_Bcast(&expcoef[loffset+moffset][1], NMAX, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
	  moffset++;
	}
	else {
	  MPI_Bcast(&expcoef[loffset+moffset][1], NMAX, MPI_DOUBLE,
		     0, MPI_COMM_WORLD);
	  MPI_Bcast(&expcoef[loffset+moffset+1][1], NMAX, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
	  moffset+=2;
	}
      }
  }

}

void SphericalSL::pca_hall
(int compute, Matrix& expcoef, Matrix*& cc, Matrix& normM)
{
  int l, Ldim;
  double dd, fac, fac02, var, b;

  static Vector smth;
  static Vector *weight;
  static Vector *b_Hall;
  static Vector inv;
  static Vector eval;
  static Vector cuml;
  static Matrix *evec;
  static Matrix Tevec;
  static Matrix sqnorm;

  static Matrix covar;

  static int setup = 0;
  
  if (!setup) {

    Ldim = LMAX*(LMAX + 2) + 1;
    
    weight = new Vector [Ldim];
    b_Hall = new Vector [Ldim];
    evec = new Matrix [Ldim];
    
    for (l=0; l<Ldim; l++) {
      weight[l].setsize(1, NMAX);
      b_Hall[l].setsize(1, NMAX);
      evec[l].setsize(1, NMAX, 1, NMAX);
    }

    smth.setsize(1, NMAX);
    inv.setsize(1, NMAX);
    eval.setsize(1, NMAX);
    cuml.setsize(1, NMAX);
    Tevec.setsize(1, NMAX, 1, NMAX);
    covar.setsize(1, NMAX, 1, NMAX);
    sqnorm.setsize(0, LMAX, 1, NMAX);
      
    for (l=0; l<=LMAX; l++)
      for (int n=1; n<=NMAX; n++) sqnorm[l][n] = sqrt(normM[l][n]);

    setup = 1;
  }


  int m, loffset, moffset, n, nn, indx, L0, lm;

  L0 = 0;
  fac02 = 16.0*M_PI*M_PI;

  for (l=L0, loffset=0; l<=LMAX; loffset+=(2*l+1), l++) {

    for (m=0, moffset=0; m<=l; m++) {

      lm = l;
      indx = loffset+moffset;

      if (m==0) {

	if (compute) {

	  for(n=1; n<=NMAX; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(n=1; n<=NMAX; n++) {
	    for(nn=n; nn<=NMAX; nn++) {
	      fac = sqnorm[lm][n]*sqnorm[lm][nn];
	      covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	      if (n!=nn)
		covar[nn][n] = covar[n][nn];
	    }    
	  }

				/* Diagonalize variance */

#ifdef GHQL
	  eval = covar.Symmetric_Eigenvalues_GHQL(evec[indx]);
#else
	  eval = covar.Symmetric_Eigenvalues(evec[indx]);
#endif
	  Tevec = evec[indx].Transpose();

	  if (tk_type == 2) {
	    cuml = eval;
	    for (n=2; n<=NMAX; n++) cuml[n] += cuml[n-1];
	    var = cuml[NMAX];
	    for (n=1; n<=NMAX; n++) cuml[n] /= var;
	  }

	  for (n=1; n<=NMAX; n++) {

	    for (dd=0.0, nn=1; nn<=NMAX; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn]*sqnorm[lm][nn];

	    var = eval[n]/used - dd*dd;

	    if (tk_type == 1) {

	      if (tksmooth*var > dd*dd)
		weight[indx][n] = 0.0;
	      else
		weight[indx][n] = 1.0;

	    }
	    else if (tk_type == 2) {
	      
	      if (n==1 || cuml[n] <= tkcum)
		weight[indx][n] = 1.0;
	      else
		weight[indx][n] = 0.0;
		
	    }
	    else if (tk_type == 3) {
	      
	      weight[indx][n] = 1.0/(1.0 + var/(dd*dd + 1.0e-14));
		
	    }
	    else
		weight[indx][n] = 1.0;

	    smth[n] = dd * weight[indx][n];
	  }

	}
	else {
	  Tevec = evec[indx].Transpose();
	  for (n=1; n<=NMAX; n++) {
	    for (dd=0.0, nn=1; nn<=NMAX; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}
	    
	inv = evec[indx]*smth;
	for (n=1; n<=NMAX; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	moffset++;
      }
      else {

	if (compute) {

	  for(n=1; n<=NMAX; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(n=1; n<=NMAX; n++) {
	    for(nn=n; nn<=NMAX; nn++) {
	      fac = sqnorm[lm][n] * sqnorm[lm][nn];
	      covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	      if (n!=nn)
		covar[nn][n] = covar[n][nn];
	    }
	  }  

				/* Diagonalize variance */

#ifdef GHQL
	  eval = covar.Symmetric_Eigenvalues_GHQL(evec[indx]);
#else
	  eval = covar.Symmetric_Eigenvalues(evec[indx]);
#endif
	  Tevec = evec[indx].Transpose();

	  if (tk_type == 2) {
	    cuml = eval;
	    for (n=2; n<=NMAX; n++) cuml[n] += cuml[n-1];
	    var = cuml[NMAX];
	    for (n=1; n<=NMAX; n++) cuml[n] /= var;
	  }

	  for (n=1; n<=NMAX; n++) {

	    for (dd=0.0, nn=1; nn<=NMAX; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn]*sqnorm[lm][nn];

	    var = eval[n]/used - dd*dd;

	    if (tk_type == 1) {

	      if (tksmooth*var > dd*dd)
		weight[indx][n] = 0.0;
	      else
		weight[indx][n] = 1.0;

	    }
	    else if (tk_type == 2) {
	      
	      if (n==1 || cuml[n] <= tkcum)
		weight[indx][n] = 1.0;
	      else
		weight[indx][n] = 0.0;
		
	    }
	    else if (tk_type == 3) {
	      
	      weight[indx][n] = 1.0/(1.0 + var/(dd*dd + 1.0e-14));
		
	    }
	    else
		weight[indx][n] = 1.0;

	    smth[n] = dd * weight[indx][n];
	  }
	}
	else {
	  Tevec = evec[indx].Transpose();
	  for (n=1; n<=NMAX; n++) {
	    for (dd=0.0, nn=1; nn<=NMAX; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}
	    
	inv = evec[indx]*smth;
	for (n=1; n<=NMAX; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	indx++;

	if (compute) {

	  for(n=1; n<=NMAX; n++) {
	    b = (cc[indx][n][n]*fac02 - expcoef[indx][n]*expcoef[indx][n]) /
	      (expcoef[indx][n]*expcoef[indx][n]*used);
	    b_Hall[indx][n] = 1.0/(1.0 + b);
	  }
    
	  for(n=1; n<=NMAX; n++) {
	    for(nn=n; nn<=NMAX; nn++) {
	      fac = sqnorm[lm][n] * sqnorm[lm][nn];
	      covar[n][nn] = fac * expcoef[indx][n]*expcoef[indx][nn];
	      if (n!=nn)
		covar[nn][n] = covar[n][nn];
	    }    
	  }

				/* Diagonalize variance */

#ifdef GHQL
	  eval = covar.Symmetric_Eigenvalues_GHQL(evec[indx]);
#else
	  eval = covar.Symmetric_Eigenvalues(evec[indx]);
#endif
	  Tevec = evec[indx].Transpose();

	  if (tk_type == 2) {
	    cuml = eval;
	    for (n=2; n<=NMAX; n++) cuml[n] += cuml[n-1];
	    var = cuml[NMAX];
	    for (n=1; n<=NMAX; n++) cuml[n] /= var;
	  }

	  for (n=1; n<=NMAX; n++) {

	    for (dd=0.0, nn=1; nn<=NMAX; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn]*sqnorm[lm][nn];

	    var = eval[n]/used - dd*dd;

	    if (tk_type == 1) {

	      if (tksmooth*var > dd*dd)
		weight[indx][n] = 0.0;
	      else
		weight[indx][n] = 1.0;

	    }
	    else if (tk_type == 2) {
	      
	      if (n==1 || cuml[n] <= tkcum)
		weight[indx][n] = 1.0;
	      else
		weight[indx][n] = 0.0;
		
	    }
	    else if (tk_type == 3) {
	      
	      weight[indx][n] = 1.0/(1.0 + var/(dd*dd + 1.0e-14));
		
	    }
	    else
		weight[indx][n] = 1.0;

	    smth[n] = dd * weight[indx][n];
	  }
	}
	else {
	  Tevec = evec[indx].Transpose();
	  for (n=1; n<=NMAX; n++) {
	    for (dd=0.0, nn=1; nn<=NMAX; nn++) 
	      dd += Tevec[n][nn]*expcoef[indx][nn] * sqnorm[lm][nn];
	    smth[n] = dd * weight[indx][n];
	  }
	}

	inv = evec[indx]*smth;
	for (n=1; n<=NMAX; n++) {
	  expcoef[indx][n] = inv[n]/sqnorm[lm][n];
	  if (tk_type == 0) expcoef[indx][n] *= b_Hall[indx][n];
	}
	
	moffset += 2;
      }
    }
  }

}


void SphericalSL::dump_basis(string& dumpname)
{
  static string labels ="pot.";
  
  double rmax = 0.33*RMAX;
  int numr = 400;
  double r, dr = rmax/numr;

  for (int L=0; L<=LMAX; L++) {
    
    ostringstream outs;
    outs << "sphbasis." << L << "." << dumpname.c_str() << '\0';

    ofstream out(outs.str().c_str());
    out.precision(3);
    out.setf(ios::scientific);

    for (int i=0; i<numr; i++) {
      r = dr*(0.5+i);

      out << setw(12) << r;
      for (int n=1; n<=min<int>(NMAX, 3); n++) {
	out
	  << setw(12) << ortho->get_pot(r, L, n, 1)
	  << setw(12) << ortho->get_dens(r, L, n, 1)
	  << setw(12) << ortho->get_force(r, L, n, 1);
      }
      out << endl;
    }
    out.close();
  }

}
    
