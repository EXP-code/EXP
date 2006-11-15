// This may look like C code, but it is really -*- C++ -*-

// #define TEST

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Test QLD program (for 2-dimensional function)
 *  Distribution function for Toomre Disk
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
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/20/91
 *  updated 6/10/94
 *
 ***************************************************************************/

using namespace std;

#include <iostream>
#include <fstream>
#include <string>

#include <malloc.h>
#include <math.h>

#include <numerical.h>
#include <Vector.h>
#include <gaussQ.h>
#include <interp.h>
#include <logic.h>

#include <massmodel.h>

#include <QPDistF.h>

bool QPDistF::MassEGrid = true;
int QPDistF::ITERMAX = 1000;
double QPDistF::ITERTOL = 1.0e-6;
double QPDistF::FSIGE = 1.2;
double QPDistF::FSIGK = 2.0;

extern "C" int ql0001_(int *m,int *me,int *mmax,int *n,int *nmax,int *mnn,
            double *c,double *d,double *a,double *b,double *xl,
            double *xu,double *x,double *u,int *iout,int *ifail,
            int *iprint,double *war,int *lwar,int *iwar,int *liwar,
            double *eps1);

//======================================================================
// Utility functions
//======================================================================

double QPDistF::kernel(double x, double y, 
		       double x0, double y0, 
		       double sigmax, double sigmay)
{
  return exp( 
	     -0.5*(x-x0)*(x-x0)/(sigmax*sigmax) 
	     -0.5*(y-y0)*(y-y0)/(sigmay*sigmay) 
	     )/(2.0*M_PI*sigmax*sigmay);
}

double QPDistF::kernel_x(double x, double y, 
		       double x0, double y0, 
		       double sigmax, double sigmay)
{
  return -(x-x0)/(sigmax*sigmax) * 
    exp( 
	-0.5*(x-x0)*(x-x0)/(sigmax*sigmax) 
	-0.5*(y-y0)*(y-y0)/(sigmay*sigmay) 
	)/(2.0*M_PI*sigmax*sigmay);
}

double QPDistF::kernel_y(double x, double y, 
		       double x0, double y0, 
		       double sigmax, double sigmay)
{
  return -(y-y0)/(sigmay*sigmay) * 
    exp( 
	-0.5*(x-x0)*(x-x0)/(sigmax*sigmax) 
	-0.5*(y-y0)*(y-y0)/(sigmay*sigmay) 
	)/(2.0*M_PI*sigmax*sigmay);
}



static AxiSymModel *m;
static double find_rofm_mass;

static double find_rofm(double r)
{
  return find_rofm_mass - m->get_mass(r);
}

QPDistF::QPDistF(AxiSymModel *T, double rmmax, double remax, 
		 int egrid, int kgrid, int mgrid,
		 double sigma, double lambda, double alpha, double beta,
		 double gama, double roff, double eoff, double koff, 
		 double kmin, double kmax,
		 int nint, int numt)
{
				// Model
  t = T;
				// Parameters
  RMMAX = rmmax;
  REMAX = remax;
  EGRID = egrid;
  KGRID = kgrid;
  MGRID = mgrid;
  SIGMA = sigma;
  LAMBDA = lambda;
  ALPHA = alpha;
  BETA = beta;
  GAMA = gama;
  ROFF = roff;
  EOFF = eoff;
  KOFF = koff;
  KMIN = kmin;
  KMAX = kmax;
  NINT = nint;
  NUMT = numt;

  NJMAX=100;
  VERBOSE=0;

  df_computed = false;
}

QPDistF::QPDistF(AxiSymModel *T, string file)
{
				// Model
  t = T;

				// Read in computed model
  read_state(file);
				// Initialize spherical orbit
  orb = new SphericalOrbit(t);

  df_computed = true;
				// Check registered model with ID

				// Would like to add ID checking here
				// but haven't yet
}

QPDistF::~QPDistF(void)
{
  if (df_computed) delete orb;
}

void QPDistF::set_verbose(void)
{
  if (VERBOSE) 
    VERBOSE=0;
  else
    VERBOSE=1;
}

void QPDistF::compute_distribution(void)
{
  df_computed = true;

				// Indices
  int i, j, k;

  //
  // Set-up model
  //

  m = t;
				// Set up energy and kappa grid
  Emin = t->get_pot(t->get_min_radius());
  Emax = t->get_pot(REMAX);
  double dE = (Emax-Emin)/EGRID;
  double dK = (KMAX-KMIN)/KGRID;
  Egrid.setsize(1, EGRID);
  Kgrid.setsize(1, KGRID);
  sigma_E.setsize(1, EGRID);
  sigma_K.setsize(1, KGRID);

  if (MassEGrid) {
    double rmin, rmax;
    double Mmin = t->get_mass((rmin=t->get_min_radius()));
    double Mmax = t->get_mass((rmax=RMMAX));
    //    double dM = (Mmax - Mmin)/EGRID;
    double dM = (log(Mmax) - log(Mmax/300.0))/EGRID;

#ifdef TEST
    cerr << "QPDistF: [Mmin, Mmax] = [" << Mmin << ", " << Mmax << "]\n";
#endif

    double rcur, rlo, rhi, mcur, mlo, mhi;

    for (i=1; i<=EGRID; i++) {
      double target = Mmax/300.0 * exp(dM * i);

      rlo = rmin;
      rhi = rmax;
      for (int it=0; it<ITERMAX; it++) {
	rcur = 0.5*(rlo + rhi);
	mcur = t->get_mass(rcur);
	if (mcur<target) {
	  rlo = rcur;
	  mlo = mcur;
	}
	else {
	  rhi = rcur;
	  mhi = mcur;
	}

	if (rhi-rlo < ITERTOL) break;
      }

      rcur = 0.5*(rlo + rhi);
      Egrid[i] = t->get_pot(rcur);
      rmin = rcur;

#ifdef TEST
      cerr << "QPDistF: Egrid[" << target << "[" << rcur << "]] = " << Egrid[i] << ",  eps = " << rhi-rlo << endl;
#endif

      if (i==1)
	dE = 2.0*(Egrid[i] - Emin);
      else
	dE = Egrid[i] - Egrid[i-1];

      sigma_E[i] = SIGMA * dE * FSIGE;
    }
  }
  else {
    for (i=1; i<=EGRID; i++) {
      double fac = ((double)i - EOFF)/EGRID;
      double fac2 = pow(fac, GAMA);
      Egrid[i] = Emin + (Emax-Emin)*fac2;
      sigma_E[i] = SIGMA * dE * FSIGE * GAMA*fac2/fac;
    }
  }

  for (i=1; i<=KGRID; i++) {
    Kgrid[i] = KMIN + dK*((double)i - KOFF);
    sigma_K[i] = SIGMA * dK * FSIGK;
  }

				// Set up radial grid
  double Rmin = t->get_min_radius();
  double Mmin = t->get_mass(Rmin);
  double Mmax = t->get_mass(RMMAX);
  double Mtotal = Mmax - Mmin;
  double dM = Mtotal/MGRID;
  Vector Rgrid(1, MGRID);
  Vector Dgrid(1, MGRID);
  for (i=1; i<=MGRID; i++) {
    find_rofm_mass = Mtotal*pow(dM/Mtotal*((double)i-ROFF), BETA) + Mmin;
    Rgrid[i] = zbrent(find_rofm, Rmin, 4.0*RMMAX, 1.0e-8);
    Dgrid[i] = t->get_density(Rgrid[i]);
  }

  double vrmax, vv, th, pot, E, K, dt=0.5*M_PI/NUMT;
  LegeQuad wk(NINT);
  orb = new SphericalOrbit(t);
  Matrix *basis = new Matrix [MGRID] - 1;
  for (k=1; k<=MGRID; k++) {

    basis[k].setsize(1, EGRID, 1, KGRID);
    basis[k].zero();

    pot = t->get_pot(Rgrid[k]);
    vrmax = sqrt(2.0*(Emax - pot));

/*
    for (int iR=1; iR<=NINT; iR++) {
      vr = vrmax * wk.knot(iR);

      vtmax = sqrt(vrmax*vrmax - vr*vr);
      for (int iT=1; iT<=NINT; iT++) {
	vt = vtmax * wk.knot(iT);
	
	E = 0.5*(vr*vr + vt*vt) + pot;
	orb->new_orbit(E, 0.5);
	K = Rgrid[k]*vt/orb->Jmax();

	for (i=1; i<=EGRID; i++) {
	  for (j=1; j<=KGRID; j++) {
	    basis[k][i][j] += 4.0*wk.weight(iR)*wk.weight(iT)*vrmax*vtmax *
	      kernel(E, K, Egrid[i], Kgrid[j], sigma_E[i], sigma_K[j]);
	  }
	}
      }
    }
*/
#ifdef DEBUG
    double ktmax = 0.0;
#endif

    for (int iR=1; iR<=NINT; iR++) {
      vv = vrmax * sqrt( wk.knot(iR) );
      E = 0.5*vv*vv + pot;

      for (int iT=1; iT<=NUMT; iT++) {
	th = dt * ((double)iT - 0.5);
	orb->new_orbit(E, 0.5);
	K = Rgrid[k]*vv*sin(th)/orb->Jmax();

#ifdef DEBUG
	if (K>ktmax) ktmax = K;
#endif
	for (i=1; i<=EGRID; i++) {
	  for (j=1; j<=KGRID; j++) {
	    basis[k][i][j] += 2.0*dt*wk.weight(iR)*vrmax*vrmax *
	      kernel(E, K, Egrid[i], Kgrid[j], sigma_E[i], sigma_K[j]);
	  }
	}
      }
#ifdef TEST
      cout << "QPDistF::compute_distirbution: Kmax = " << ktmax << endl;
#endif
    }
  }
  

  //======================================================================
  // Set up for QLD
  //======================================================================

  int M=0;			// Number of constraints
  int ME=0;			// Number of equality constraints
  int MMAX=1;			// Row dimension of A
  int N=EGRID*KGRID;		// Number of variables
  int NMAX=N;			// Row dimension of C
  int MNN=M+2*N;

				// Objective function
  Matrix C(1, N, 1, N);
				// Constant vector
  Vector D(1, N);

				// Data for linear constraints
  Matrix A(1, MMAX, 1, N);
				// Constant vector for linear constraints
  Vector B(1, MMAX);

				
  Vector XL(1, N);		// lower bounds
  Vector XU(1, N);		// upper bounds
  X.setsize(1, N);		// On return, X contains the optimal
				//      solution vector

  Vector U(1, MNN);		// Lagrange multipliers on return
  int IOUT = 6;			// Fortran unit stream
//  int IFAIL;			// return flag
  int IPRINT = 1;		// print flag (>0 -> output)

				// work array
  int LWAR = 3*NMAX*NMAX/2 + 10*NMAX + 2*MMAX;
  double* WAR = new double [LWAR];

				// work array
  int LIWAR = N;
  int* IWAR = new int [LIWAR];
  double EPS = 1.0e-14;		// machine precision

  IWAR[0] = 1;			// Standard Hessian matrix


//======================================================================
// Set up for fitting problem
//======================================================================

  
				// Objective function
  C.zero();
  int ix, iy, jx, jy;
  
  for (i=1, ix=1; ix<=EGRID; ix++) {
    for (iy=1; iy<=KGRID; iy++) {

      for (j=1, jx=1; jx<=EGRID; jx++) {
	for (jy=1; jy<=KGRID; jy++) {
	  
	  for (k=1; k<=MGRID; k++)
	    C[i][j] += basis[k][ix][iy]*basis[k][jx][jy];

	  j++;
	}
      }
      i++;
    }
  }

  Matrix C0(C);

  if (LAMBDA>1.0e-20) {

    for (i=1, ix=1; ix<=EGRID; ix++) {
      for (iy=1; iy<=KGRID; iy++) {

	for (j=1, jx=1; jx<=EGRID; jx++) {
	  for (jy=1; jy<=KGRID; jy++) {
	    C[i][j] += LAMBDA * pow(Kgrid[iy]*Kgrid[jy], ALPHA);
	    j++;
	  }
	}
	i++;
      }
    }

  }

				// And constant vector
  D.zero();
  for (i=1, ix=1; ix<=EGRID; ix++) {
    for (iy=1; iy<=KGRID; iy++) {
      for (k=1; k<=MGRID; k++)
	D[i] -= Dgrid[k] *  basis[k][ix][iy];
      i++;
    }
  }
  
				// Constant
  double constant=0.0;
  for (k=1; k<=MGRID; k++)
    constant += Dgrid[k] * Dgrid[k];
     
				// Limits
  for (i=1, ix=1; ix<=EGRID; ix++) {
    for (iy=1; iy<=KGRID; iy++) {
      XL[i] = 0.0;
      XU[i] = 1.0e8;
      i++;
    }
  }
  
  delete [] (basis + 1);


//======================================================================
// Convert matrices and vectors to arrays for passing to footron
//======================================================================

  double* c = convert(C);
  double* d = convert(D);
  double* a = convert(A, MMAX);
  double* b = convert(B, MMAX);
  double* xl = convert(XL);
  double* xu = convert(XU);
  double* x = convert(X);
  double* u = convert(U);

  ql0001_(&M, &ME, &MMAX, &N, &NMAX, &MNN, c, d, a, b, xl, xu, x, u, &IOUT,
	  &IFAIL, &IPRINT, WAR, &LWAR, IWAR, &LIWAR, &EPS);

				// Copy output
  X = Vector(1, N, x-1);
  U = Vector(1, MNN, u-1);

				// Delete temp arrays
  delete [] WAR;
  delete [] IWAR;

  free( c );
  free( d );
  free( a );
  free( b );
  free( xl );
  free( xu );
  free( x );
  free( u );

//======================================================================
//======================================================================

				// Diagnostic output

  obj0 = 0.5*constant;
  obj = obj0;
  for (i=1; i<=N; i++) {
    for (j=1; j<=N; j++) {
      obj0 += 0.5 * X[i] * C0[i][j] * X[j];
      obj  += 0.5 * X[i] * C[i][j] *  X[j];
    }
    obj0 += D[i]*X[i];
    obj  += D[i]*X[i];
  }

  if (VERBOSE) {
    cout << "Solution vector: " << endl;
    for (i=1; i<=N; i++)
      if (fabs(X[i])>1.0e-10) cout << i << "> " << X[i] << endl;
    cerr << "Objective0: " << obj0 << endl;
    cerr << "Objective:  " << obj  << endl;
    cerr << "Condition:  " << IFAIL << endl;
  }

  if (IFAIL) bomb("distribution function solution did not converge");

//======================================================================
// Set up JMax grid for distribution function
//======================================================================
  
  TOLE=0.005*(Emax-Emin);
  dE = (Emax-TOLE-Emin)/(NJMAX-1);
  JMAXE.setsize(1, NJMAX);
  JMAX.setsize(1, NJMAX);
  JMAX2.setsize(1, NJMAX);
  for (i=1; i<=NJMAX; i++) {
    JMAXE[i] = Emin + TOLE + dE*(i-1);
    orb->new_orbit(JMAXE[i], 0.5);
    JMAX[i] = orb->Jmax();
  }
  Spline(JMAXE, JMAX, 1.0e30, 1.0e30, JMAX2);

}


double* QPDistF::convert(Matrix& a, int mm, int nn)
{
   int m = a.getnrows();
   int n = a.getncols();
   int rlow = a.getrlow();
   int clow = a.getclow();

   mm = m>mm ? m : mm;
   nn = n>nn ? n : nn;

   double* temp = (double *)calloc(mm*nn, sizeof(double));

   for (int i=0; i<n; i++)	// loop thru columns
     for (int j=0; j<m; j++)	// loop thru rows
       temp[mm*i+j] = a[j+rlow][i+clow];

   return temp;
}

double* QPDistF::convert(Vector& a, int nn)
{
  int low = a.getlow();
  int n = a.getlength();
  
  nn = n>nn ? n : nn;

  double* temp = (double *)calloc(nn, sizeof(double));

  for (int i=0; i<n; i++)
    temp[i] = a[i+low];

   return temp;
}

double QPDistF::distf(double E, double L)
{
  if (!df_computed) compute_distribution();

  double jmax;
  if (E < JMAXE[JMAXE.gethigh()] && E > JMAXE[JMAXE.getlow()])
    Splint1(JMAXE, JMAX, JMAX2, E, jmax, 1);
  else {
    orb->new_orbit(E, 0.5);
    jmax = orb->Jmax();
  }
  return distf_EK(E, L/jmax);
}

double QPDistF::dfdE(double E, double L)
{
  if (!df_computed) compute_distribution();

  double dfde, dfdk, jmax, djmax;

  Splint2(JMAXE, JMAX, JMAX2, E, jmax, djmax, 1);

  if (E > Emin+TOLE) {
    orb->new_orbit(E, 0.5);
    jmax = orb->Jmax();
  }
  double K = L/jmax;

  dfde = dfdE_EK(E, K);
  dfdk = dfdK_EK(E, K);
  return dfde - dfdk*djmax/jmax*K;
}


double QPDistF::dfdL(double E, double L)
{
  if (!df_computed) compute_distribution();

  double jmax;

  if (E < JMAXE[JMAXE.gethigh()] && E > JMAXE[JMAXE.getlow()])
    Splint1(JMAXE, JMAX, JMAX2, E, jmax, 1);
  else {
    orb->new_orbit(E, 0.5);
    jmax = orb->Jmax();
  }
  double K = L/jmax;

  return dfdK_EK(E, K)/jmax;
}


double QPDistF::distf_EK(double E, double K)
{
  if (!df_computed) compute_distribution();

  double ans = 0.0;
  int i=1;

  for (int ix=1; ix<=EGRID; ix++) {
    for (int iy=1; iy<=KGRID; iy++) {
      if (X[i] > 1.0e-10)
	ans += X[i] * 
	  kernel(E, K, Egrid[ix], Kgrid[iy], sigma_E[ix], sigma_K[iy]);
      i++;
    }
  }
  
  return ans;
}

double QPDistF::dfdE_EK(double E, double K)
{
  if (!df_computed) compute_distribution();

  double ans = 0.0;
  int i=1;

  for (int ix=1; ix<=EGRID; ix++) {
    for (int iy=1; iy<=KGRID; iy++) {
      if (X[i] > 1.0e-10)
	ans += X[i] * 
	  kernel_x(E, K, Egrid[ix], Kgrid[iy], sigma_E[ix], sigma_K[iy]);
      i++;
    }
  }
  
  return ans;
}

double QPDistF::dfdK_EK(double E, double K)
{
  if (!df_computed) compute_distribution();

  double ans = 0.0;
  int i=1;

  for (int ix=1; ix<=EGRID; ix++) {
    for (int iy=1; iy<=KGRID; iy++) {
      if (X[i] > 1.0e-10)
	ans += X[i] * 
	  kernel_y(E, K, Egrid[ix], Kgrid[iy], sigma_E[ix], sigma_K[iy]);
      i++;
    }
  }
  
  return ans;
}

// Write out all the necessary information to recreate the DF from a
// file

void QPDistF::write_state(string& name)
{
  int i;

  ofstream out(name.c_str());
  if (!out) {
    cerr << "Couldn't open <" << name << "> to save state!\n";
    return;
  }

  out.write((char *)&RMMAX, sizeof(double));
  out.write((char *)&REMAX, sizeof(double));
  out.write((char *)&EGRID, sizeof(int));
  out.write((char *)&KGRID, sizeof(int));
  out.write((char *)&MGRID, sizeof(int));
  out.write((char *)&SIGMA, sizeof(double));
  out.write((char *)&LAMBDA, sizeof(double));
  out.write((char *)&ALPHA, sizeof(double));
  out.write((char *)&BETA, sizeof(double));
  out.write((char *)&GAMA, sizeof(double));
  out.write((char *)&ROFF, sizeof(double));
  out.write((char *)&EOFF, sizeof(double));
  out.write((char *)&KOFF, sizeof(double));
  out.write((char *)&NINT, sizeof(double));
  out.write((char *)&NUMT, sizeof(double));

  for (i=1; i<=EGRID; i++) out.write((char *)&Egrid[i], sizeof(double));
  for (i=1; i<=EGRID; i++) out.write((char *)&sigma_E[i], sizeof(double));
  for (i=1; i<=KGRID; i++) out.write((char *)&Kgrid[i], sizeof(double));
  for (i=1; i<=KGRID; i++) out.write((char *)&sigma_K[i], sizeof(double));
  for (i=1; i<=EGRID*KGRID; i++) out.write((char *)&X[i], sizeof(double));
  out.write((char *)&obj0, sizeof(double));
  out.write((char *)&obj , sizeof(double));

  out.write((char *)&NJMAX , sizeof(int));
  out.write((char *)&Emin, sizeof(double));
  out.write((char *)&Emax, sizeof(double));
  out.write((char *)&TOLE, sizeof(double));
  for (i=1; i<=NJMAX; i++) out.write((char *)&JMAXE[i], sizeof(double));
  for (i=1; i<=NJMAX; i++) out.write((char *)&JMAX[i],  sizeof(double));
  for (i=1; i<=NJMAX; i++) out.write((char *)&JMAX2[i], sizeof(double));

}

// Reinitialize the DF from a saved state

void QPDistF::read_state(string& name)
{
  int i;

  ifstream in(name.c_str());
  if (!in) {
    cerr << "Couldn't open <" << name << "> to read state!\n";
    exit(-1);
  }

  in.read((char *)&RMMAX, sizeof(double));
  in.read((char *)&REMAX, sizeof(double));
  in.read((char *)&EGRID, sizeof(int));
  in.read((char *)&KGRID, sizeof(int));
  in.read((char *)&MGRID, sizeof(int));
  in.read((char *)&SIGMA, sizeof(double));
  in.read((char *)&LAMBDA, sizeof(double));
  in.read((char *)&ALPHA, sizeof(double));
  in.read((char *)&BETA, sizeof(double));
  in.read((char *)&GAMA, sizeof(double));
  in.read((char *)&ROFF, sizeof(double));
  in.read((char *)&EOFF, sizeof(double));
  in.read((char *)&KOFF, sizeof(double));
  in.read((char *)&NINT, sizeof(double));
  in.read((char *)&NUMT, sizeof(double));

  Egrid.setsize(1, EGRID);
  for (i=1; i<=EGRID; i++) in.read((char *)&Egrid[i], sizeof(double));
  sigma_E.setsize(1, EGRID);
  for (i=1; i<=EGRID; i++) in.read((char *)&sigma_E[i], sizeof(double));
  Kgrid.setsize(1, KGRID);
  for (i=1; i<=KGRID; i++) in.read((char *)&Kgrid[i], sizeof(double));
  sigma_K.setsize(1, KGRID);
  for (i=1; i<=KGRID; i++) in.read((char *)&sigma_K[i], sizeof(double));
  X.setsize(1, EGRID*KGRID);
  for (i=1; i<=EGRID*KGRID; i++) in.read((char *)&X[i], sizeof(double));
  in.read((char *)&obj0, sizeof(double));
  in.read((char *)&obj , sizeof(double));

  in.read((char *)&NJMAX , sizeof(int));
  in.read((char *)&Emin, sizeof(double));
  in.read((char *)&Emax, sizeof(double));
  in.read((char *)&TOLE, sizeof(double));
  JMAXE.setsize(1, NJMAX);
  JMAX.setsize(1, NJMAX);
  JMAX2.setsize(1, NJMAX);
  for (i=1; i<=NJMAX; i++) in.read((char *)&JMAXE[i], sizeof(double));
  for (i=1; i<=NJMAX; i++) in.read((char *)&JMAX[i],  sizeof(double));
  for (i=1; i<=NJMAX; i++) in.read((char *)&JMAX2[i], sizeof(double));

}

