/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Test QLD program (for 2-dimensional distribution functions)
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

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>

#include "numerical.H"
#include "gaussQ.H"
#include "interp.H"

#include "massmodel.H"

#include "QPDistF.H"

bool QPDistF::MassEGrid    = true;
bool QPDistF::MassLinear   = true;
int QPDistF::ITERMAX       = 1000;
double QPDistF::ITERTOL    = 1.0e-6;
double QPDistF::FSIGE      = 1.2;
double QPDistF::FSIGK      = 2.0;

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

double QPDistF::kernel_x2(double x, double y, 
			  double x0, double y0, 
			  double sigmax, double sigmay)
{
  return ( (x-x0)*(x-x0)/(sigmax*sigmax)  - 1.0)/(sigmax*sigmax) *
    exp( 
	-0.5*(x-x0)*(x-x0)/(sigmax*sigmax) 
	-0.5*(y-y0)*(y-y0)/(sigmay*sigmay) 
	 )/(2.0*M_PI*sigmax*sigmay);
}

double QPDistF::kernel_xy(double x, double y, 
			  double x0, double y0, 
			  double sigmax, double sigmay)
{
  return (x-x0)*(y-y0)/(sigmax*sigmax*sigmay*sigmay) *
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

double QPDistF::kernel_y2(double x, double y, 
			  double x0, double y0, 
			  double sigmax, double sigmay)
{
  return ( (y-y0)*(y-y0)/(sigmay*sigmay)  - 1.0)/(sigmay*sigmay) *
    exp( 
	-0.5*(x-x0)*(x-x0)/(sigmax*sigmax) 
	-0.5*(y-y0)*(y-y0)/(sigmay*sigmay) 
	 )/(2.0*M_PI*sigmax*sigmay);
}



static AxiSymModPtr m;
static double find_rofm_mass;

static double find_rofm(double r)
{
  return find_rofm_mass - m->get_mass(r);
}

QPDistF::QPDistF(AxiSymModPtr T, AxiSymModPtr H, 
		 double rmmax, double remax, 
		 int egrid, int kgrid, int mgrid,
		 double sigma, double lambda, double alpha, double beta,
		 double gama, double roff, double eoff, double koff, 
		 double kmin, double kmax,
		 int nint, int numt)
{
				// Model
  t       = std::make_shared<DiskWithHalo>(T, H);
  nt      = true;
				// Parameters
  RMMAX   = rmmax;
  REMAX   = remax;
  EGRID   = egrid;
  KGRID   = kgrid;
  MGRID   = mgrid;
  SIGMA   = sigma;
  LAMBDA  = lambda;
  ALPHA   = alpha;
  BETA    = beta;
  GAMA    = gama;
  ROFF    = roff;
  EOFF    = eoff;
  KOFF    = koff;
  KMIN    = kmin;
  KMAX    = kmax;
  NINT    = nint;
  NUMT    = numt;

  NJMAX   = 100;
  verbose = 0;

  df_computed = false;
}

QPDistF::QPDistF(AxiSymModPtr T,
		 double rmmax, double remax, 
		 int egrid, int kgrid, int mgrid,
		 double sigma, double lambda, double alpha, double beta,
		 double gama, double roff, double eoff, double koff, 
		 double kmin, double kmax,
		 int nint, int numt)
{
				// Model
  t       = T;
  nt      = false;
				// Parameters
  RMMAX   = rmmax;
  REMAX   = remax;
  EGRID   = egrid;
  KGRID   = kgrid;
  MGRID   = mgrid;
  SIGMA   = sigma;
  LAMBDA  = lambda;
  ALPHA   = alpha;
  BETA    = beta;
  GAMA    = gama;
  ROFF    = roff;
  EOFF    = eoff;
  KOFF    = koff;
  KMIN    = kmin;
  KMAX    = kmax;
  NINT    = nint;
  NUMT    = numt;

  NJMAX   = 100;
  verbose = 0;

  df_computed = false;
}

QPDistF::QPDistF(AxiSymModPtr T, std::string file)
{
				// Model
  t  = T;
  nt = false;

				// Read in computed model
  read_state(file);
				// Initialize spherical orbit
  orb = std::make_shared<SphericalOrbit>(t);

  df_computed = true;
				// Check registered model with ID

				// Would like to add ID checking.Here
				// but haven't yet
}

QPDistF::~QPDistF(void)
{
  // Nothing
}

void QPDistF::set_verbose(void)
{
  if (verbose) 
    verbose=0;
  else
    verbose=1;
}

void QPDistF::set_verbose(int verbosity)
{
  const int vmin = 0;
  const int vmax = 5;

  if (verbosity<vmin) 
    verbose = vmin;
  else if (verbosity>vmax) 
    verbose = vmax;
  else 
    verbose = verbosity;

}

void QPDistF::compute_distribution(void)
{
  df_computed = true;

  //
  // Set-up model
  //

  double Rmin = t->get_min_radius();
  double Rmax = RMMAX;
  double Mmax = t->get_mass(Rmax);
  double Mmin = max<double>(t->get_mass(Rmin), 1.0e-6*Mmax);

  m = t;
				// Set up energy and kappa grid
  Emin = t->get_pot(Rmin);
  Emax = t->get_pot(Rmax);

  double dE = (Emax-Emin)/EGRID;
  double dK = (KMAX-KMIN)/KGRID;

  Egrid.resize(EGRID);
  Kgrid.resize(KGRID);

  sigma_E.resize(EGRID);
  sigma_K.resize(KGRID);

  if (MassEGrid) {
    double dM;

    if (MassLinear)
      dM = (Mmax - Mmin)/EGRID;
    else
      dM = (log(Mmax) - log(Mmin))/EGRID;

    if (verbose>3)
      cerr << "QPDistF: [Mmin, Mmax] = [" << Mmin 
	   << ", " << Mmax << "]" << endl;

    double rcur, rlo, rhi, mcur, mlo, mhi, target;
    double rmin = Rmin, rmax = Rmax;

    for (int i=0; i<EGRID; i++) {
      if (MassLinear)
	target = Mmin + dM*((double)i+0.5);
      else
	target = Mmin * exp(dM * i);
      
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

      if (verbose>3)
	cerr << "QPDistF: Egrid[ " << setw(9) << target
	     << "[" << setw(9) << rcur << "]] = " << setw(9) << Egrid[i]
	     << ",  eps = " << rhi - rlo << endl;

      if (i==0)
	dE = 2.0*(Egrid[i] - Emin);
      else
	dE = Egrid[i] - Egrid[i-1];

      sigma_E[i] = SIGMA * dE * FSIGE;
    }

  } else {

    for (int i=0; i<EGRID; i++) {
      double fac  = ((double)(i+1) - EOFF)/EGRID;
      double fac2 = pow(fac, GAMA);
      Egrid[i]    = Emin + (Emax-Emin)*fac2;
      sigma_E[i]  = SIGMA * dE * FSIGE * GAMA*fac2/fac;
    }

  }

  for (int i=0; i<KGRID; i++) {
    Kgrid[i] = KMIN + dK*((double)(i+1) - KOFF);
    sigma_K[i] = SIGMA * dK * FSIGK;
  }

				// Set up radial grid
  double Mtotal = Mmax - Mmin;
  double dM = Mtotal/MGRID;

  Eigen::VectorXd Rgrid(MGRID);
  Eigen::VectorXd Dgrid(MGRID);

  for (int i=0; i<MGRID; i++) {
    find_rofm_mass = Mtotal*pow(dM/Mtotal*((double)(i+1)-ROFF), BETA) + Mmin;

    Rgrid[i] = zbrent(find_rofm, 0.25*Rmin, 4.0*Rmax, 1.0e-8);
    Dgrid[i] = t->get_density(Rgrid[i]);
  }

  double vrmax, vv, th, pot, E, K, dt=0.5*M_PI/NUMT;
  LegeQuad wk(NINT);
  orb = std::make_shared<SphericalOrbit>(t);
  std::vector<Eigen::MatrixXd> basis(MGRID);
  for (int k=0; k<MGRID; k++) {

    basis[k].resize(EGRID, KGRID);
    basis[k].setZero();

				// Gravitational potential at R
    pot = t->get_pot(Rgrid[k]);
				// Maximum ("escape") velocity at R
    vrmax = sqrt(2.0*(Emax - pot));

    double ktmax = 0.0;

    if (t->dof()==2) {

      for (int iR=0; iR<NINT; iR++) {

	vv = vrmax * sqrt( wk.knot(iR) );
	E = 0.5*vv*vv + pot;

	for (int iT=1; iT<=NUMT; iT++) {
	  th = dt * ((double)iT - 0.5);
	  orb->new_orbit(E, 0.5);
	  K = Rgrid[k]*vv*sin(th)/orb->Jmax();

	  if (verbose>4) if (K>ktmax) ktmax = K;

	  for (int i=0; i<EGRID; i++) {
	    for (int j=0; j<KGRID; j++) {
	      basis[k](i, j) += 2.0*dt*wk.weight(iR)*vrmax*vrmax *
		kernel(E, K, Egrid[i], Kgrid[j], sigma_E[i], sigma_K[j]);
	    }
	  }
	  
	}
	if (verbose>4)
	  cout << "QPDistF::compute_distribution: Kmax = " << ktmax << endl;
      }
    } else if (t->dof()==3) {

      for (int ix=0; ix<NINT; ix++) {
	double x = wk.knot(ix);

	for (int iy=0; iy<NINT; iy++) {
	  double y = wk.knot(iy);

	  double E = pot + 0.5*vrmax*vrmax*(x*x + (1.0-x*x)*y*y);
	  double J = vrmax*sqrt(1.0 - x*x)*y*Rgrid[k];

	  orb->new_orbit(E, 0.5);
	  K = J/orb->Jmax();

	  if (verbose>4) if (K>ktmax) ktmax = K;

	  double fac = wk.weight(ix)*wk.weight(iy) * 4.0*M_PI *
	    vrmax*vrmax*vrmax * (1.0 - x*x)*y;

	  for (int i=0; i<EGRID; i++) {
	    for (int j=0; j<KGRID; j++) {
	      basis[k](i, j) += fac *
		kernel(E, K, Egrid[i], Kgrid[j], sigma_E[i], sigma_K[j]);
	    }
	  }
	}
      }
      
      if (verbose>4)
	cout << "QPDistF::compute_distribution: Kmax = " << ktmax << endl;
    } else {
      cerr << "QPDistF: dof=" << t->dof() << ", must be 2 or 3" << endl;
    }

    /*
      for (int iR=0; iR<NINT; iR++) {

	vv = vrmax * sqrt( wk.knot(iR) );
	E = 0.5*vv*vv + pot;

	for (int iT=1; iT<=NUMT; iT++) {
	  th = dt * ((double)iT - 0.5);
	  orb->new_orbit(E, 0.5);
	  K = Rgrid[k]*vv*sin(th)/orb->Jmax();

	  if (verbose>4) if (K>ktmax) ktmax = K;

	  for (i=0; i<EGRID; i++) {
	    for (j=0; j<KGRID; j++) {
	      basis[k][i][j] += 2.0*dt*wk.weight(iR)*vrmax*vrmax * 
		2.0*M_PI*vv*sin(th) *
		kernel(E, K, Egrid[i], Kgrid[j], sigma_E[i], sigma_K[j]);
	    }
	  }
	  
	}
	if (verbose>4)
	  cout << "QPDistF::compute_distribution: Kmax = " << ktmax << endl;
      }
    } else {
      cerr << "QPDistF: dof=" << t->dof() << ", must be 2 or 3" << endl;
    }
    */
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
  Eigen::MatrixXd C(N, N);
				// Constant vector
  Eigen::VectorXd D(N);

				// Data for linear constraints
  Eigen::MatrixXd A(MMAX, N);
				// Constant vector for linear constraints
  Eigen::VectorXd B(MMAX);

				
  Eigen::VectorXd XL(N);	// lower bounds
  Eigen::VectorXd XU(N);	// upper bounds
  X.resize(N);			// On return, X contains the optimal
				// solution vector

  Eigen::VectorXd U(MNN);	// Lagrange multipliers on return
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
  C.setZero();
  
  for (int i=0, ix=0; ix<EGRID; ix++) {
    for (int iy=0; iy<KGRID; iy++) {

      for (int j=0, jx=0; jx<EGRID; jx++) {
	for (int jy=0; jy<KGRID; jy++) {
	  
	  for (int k=0; k<MGRID; k++) {
	    C(i, j) += basis[k](ix, iy)*basis[k](jx, jy);
	    if (std::isnan(basis[k](ix, iy))) {
	      cout << "Basis NaN k=" << k << " ix=" << ix << " iy=" << iy
		   << endl;
	      exit(-1);
	    }
	    if (std::isnan(basis[k](jx, jy))) {
	      cout << "Basis NaN k=" << k << " jx=" << jx << " jy=" << jy
		   << endl;
	      exit(-1);
	    }
	  
	  }
	  j++;
	}
      }
      i++;
    }
  }

  Eigen::MatrixXd C0(C);

  if (LAMBDA>1.0e-20) {

    for (int i=0, ix=0; ix<EGRID; ix++) {
      for (int iy=0; iy<KGRID; iy++) {
	
	for (int j=0, jx=0; jx<EGRID; jx++) {
	  for (int jy=0; jy<KGRID; jy++) {
	    C(i, j) += LAMBDA * pow(Kgrid[iy]*Kgrid[jy], ALPHA);
	    j++;
	  }
	}
	i++;
      }
    }

  }

				// And constant vector
  D.setZero();
  for (int i=0, ix=0; ix<EGRID; ix++) {
    for (int iy=0; iy<KGRID; iy++) {
      for (int k=0; k<MGRID; k++)
	D[i] -= Dgrid[k] *  basis[k](ix, iy);
      i++;
    }
  }
  
				// Constant
  double constant=0.0;
  for (int k=0; k<MGRID; k++)
    constant += Dgrid[k] * Dgrid[k];
     
				// Limits
  for (int i=0, ix=0; ix<EGRID; ix++) {
    for (int iy=0; iy<KGRID; iy++) {
      XL[i] = 0.0;
      XU[i] = 1.0e8;
      i++;
    }
  }
  
  if (verbose>2) {

    cout.precision(4);

    cout << endl 
	 << "-------------------\n"
	 << "Objective function:\n"
	 << "-------------------\n";
      
    for (int i=0; i<C.rows(); i++) {
      for (int j=0; j<C.cols(); j++) 
	cout << setw(16) << C(i, j);
      cout << endl;
    }

    cout << endl 
	 << "----------------" << endl
	 << "Constant vector:" << endl
	 << "----------------" << endl;
      
    for (int i=0; i<D.size(); i++)
      cout << setw(4) << i << setw(16) << D[i];
    cout << endl;

  }


//======================================================================
// Convert matrices and vectors to arrays for passing to footron
//======================================================================

  double* c  = C.data();
  double* d  = D.data();
  double* a  = A.data();
  double* b  = B.data();
  double* xl = XL.data();
  double* xu = XU.data();
  double* x  = X.data();
  double* u  = U.data();

  ql0001_(&M, &ME, &MMAX, &N, &NMAX, &MNN, c, d, a, b, xl, xu, x, u, &IOUT,
	  &IFAIL, &IPRINT, WAR, &LWAR, IWAR, &LIWAR, &EPS);

				// Delete temp arrays
  delete [] WAR;
  delete [] IWAR;

//======================================================================
//======================================================================

				// Diagnostic output

  obj0 = 0.5*constant;
  obj = obj0;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      obj0 += 0.5 * X[i] * C0(i, j) * X[j];
      obj  += 0.5 * X[i] * C (i, j) *  X[j];
    }
    obj0 += D[i]*X[i];
    obj  += D[i]*X[i];
  }

  if (verbose>0) {
    cout << "----------------" << endl;
    cout << "Solution vector:" << endl;
    cout << "----------------" << endl;
    int i=0;
    for (int ix=0; ix<EGRID; ix++) {
      for (int iy=0; iy<KGRID; iy++) {
	if (X[i] > 1.0e-12)
	  cout << setw(4) << ix
	       << setw(4) << iy << "> "
	       << setw(16) << X[i] 
	       << setw(16) << Egrid[ix]
	       << setw(16) << Kgrid[iy]
	       << endl;
	else X[i] = 0.0;
	i++;
      }
    }

    cout << "Objective value0 : " << obj0 << endl;
    cout << "Objective value  : " << obj  << endl;
    cout << "Condition        : " << IFAIL << endl;
  }

  if (IFAIL) {
    char msg[] = "distribution function solution did not converge";
    bomb(msg);
  }

//======================================================================
// Set up JMax grid for distribution function
//======================================================================
  
  TOLE=0.005*(Emax-Emin);
  dE = (Emax-TOLE-Emin)/(NJMAX-1);
  JMAXE.resize(NJMAX);
  JMAX .resize(NJMAX);
  JMAX2.resize(NJMAX);
  for (int i=0; i<NJMAX; i++) {
    JMAXE[i] = Emin + TOLE + dE*(i-1);
    orb->new_orbit(JMAXE[i], 0.5);
    JMAX[i] = orb->Jmax();
  }
  Spline(JMAXE, JMAX, 1.0e30, 1.0e30, JMAX2);

}


void QPDistF::make_cdf(int ENUM, int KNUM, double KTOL)
{
  NUME = ENUM;
  NUMK = KNUM;
  TOLK = KTOL;
  
  double E, K;
  double dE = (Emax - Emin)/NUME;
  double dK = (1.0 - 2.0*TOLK)/NUMK;
  pdf = vector<double>(NUME*NUMK);
  for (int ix=0; ix<NUME; ix++) {
    E = Emin + dE*ix;
    for (int iy=0; iy<NUMK; iy++) {
      K = TOLK + dK*iy;
      orb->new_orbit(E, K);
      pdf[ix*NUMK+iy] = distf_EK(E, K)*orb->Jmax()/orb->get_freq(1) * dE*dK;
    }
  }

				// First sum columns 
  cdf = pdf;
				// [This can be done in one pass but
				// this construction is easier to
				// understand]
  for (int ix=0; ix<NUME; ix++)
    for (int iy=1; iy<NUMK; iy++)
      cdf[ix*NUMK+iy] += cdf[ix*NUMK+iy-1];
				// Next sum rows
  for (int iy=0; iy<NUMK; iy++)
    for (int ix=1; ix<NUME; ix++)
      cdf[ix*NUMK+iy] += cdf[(ix-1)*NUMK+iy];

				// Normalize
  double norm = 1.0/cdf[NUME*NUMK-1];
  for (int iy=0; iy<NUMK; iy++)
    for (int ix=1; ix<NUME; ix++)
      cdf[ix*NUMK+iy] *= norm;

				// Make multimap for realization
  for (int ix=0; ix<NUME; ix++)
    for (int iy=0; iy<NUMK; iy++)
      realz.insert(elem(cdf[ix*NUMK+iy], ip(ix, iy)));
}


pair<double, double> QPDistF::gen_EK(double r1, double r2)
{
  multimap<double, ip>::iterator jj = realz.upper_bound(r1);
  int ix = jj->second.first;
  int iy = jj->second.second;
  
  double dE = (Emax - Emin)/NUME;
  double dK = (1.0 - 2.0*TOLK)/NUMK;
  
  vector<double> x(2), y(2);

  double a = cdf[(ix-1)*NUMK+(iy-1)];
  double b = cdf[(ix-0)*NUMK+(iy-1)];
  double c = cdf[(ix-1)*NUMK+(iy-0)];
  double d = cdf[(ix-0)*NUMK+(iy-0)];

  // Locate intercepts: vert-horiz segment
  //
  if (r1 >= a && r1 < c) {
    x[0] = Emin + dE*(ix-1);
    y[0] = TOLK + dK*((r1-a)/(c-a)+iy-1);
  } else if (r1 >= c && r1 < d) {
    x[0] = Emin + dE*((r1-c)/(d-c)+ix-1);
    y[0] = TOLK + dK*(iy-0);
  } else {
    x[0] = Emin + dE*(ix-0);
    y[0] = TOLK + dK*(iy-0);
  }

  // Locate intercepts: horiz-vert segment
  //
  if (r1 >= a && r1 < b) {
    x[1] = Emin + dE*(ix-1);
    y[1] = TOLK + dK*((r1-a)/(b-a)+iy-1);
  } else if (r1 >= c && r1 < d) {
    x[1] = Emin + dE*((r1-b)/(d-b)+ix-1);
    y[1] = TOLK + dK*(iy-0);
  } else {
    x[1] = Emin + dE*(ix-0);
    y[1] = TOLK + dK*(iy-0);
  }

  return pair<double, double> (x[0] + (x[1] - x[0])*r2,
			       y[0] + (y[1] - y[0])*r2);
  
}


void QPDistF::dump_cdf(const string& file)
{
  ofstream out(file.c_str());
  if (!out) return;

  double E, K;
  double dE = (Emax -  Emin)/NUME;
  double dK = (1.0 - 2*TOLK)/NUMK;

  for (int ix=0; ix<NUME; ix++) {
    E = Emin + dE*ix;
    for (int iy=0; iy<NUMK; iy++) {
      K = TOLK + dK*iy;
      out << setw(18) << E << setw(18) << K 
	  << setw(18) << distf_EK(E, K)
	  << setw(18) << pdf[ix*NUMK+iy]
	  << setw(18) << cdf[ix*NUMK+iy] << endl;
    }
    out << endl;
  }

  // DEBUG
  ofstream tst("cdf.tmp");
  multimap<double, ip>::iterator jj = realz.begin();
  for (int j=0; jj!=realz.end(); j++, jj++) {
    tst << setw(4) << j 
	<< setw(4) << jj->second.first
	<< setw(4) << jj->second.second
	<< setw(18) << Emin + dE*jj->second.first
	<< setw(18) << TOLK + dK*jj->second.second
	<< setw(18) << jj->first
	<< endl;
  }
  // END
}


double QPDistF::distf(double E, double L)
{
  if (!df_computed) compute_distribution();

  double jmax;
  if (E < JMAXE[JMAXE.size()-1] && E > JMAXE[0])
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
  return dfde - dfdk*djmax/jmax*K + 2.0*distf_EK(E, K)*djmax/jmax;
}


double QPDistF::d2fdE2(double E, double L)
{
  if (!df_computed) compute_distribution();

  double jmax, djmax, djmax2;

  Splint3(JMAXE, JMAX, JMAX2, E, jmax, djmax, djmax2, 1);

  if (E > Emin+TOLE) {
    orb->new_orbit(E, 0.5);
    jmax = orb->Jmax();
  }
  double K = L/jmax;
  double dj  = djmax/jmax;
  double dj2 = dj*dj;

  return d2fdE2_EK(E, K) + 4.0*dfdE_EK(E, K)*dj - 
    2.0*d2fdEK_EK(E, K)*K*dj - dfdK_EK(E, K)*(K*djmax2 + 2.0*K*dj2) +
    d2fdK2_EK(E, K)*K*K*dj2 + distf_EK(E, K)*(2.0*dj2 + djmax2/jmax);
}


double QPDistF::dfdL(double E, double L)
{
  if (!df_computed) compute_distribution();

  double jmax;

  if (E < JMAXE[JMAXE.size()-1] && E > JMAXE[0])
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
  int i=0;

  for (int ix=0; ix<EGRID; ix++) {
    for (int iy=0; iy<KGRID; iy++) {
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
  int i=0;

  for (int ix=0; ix<EGRID; ix++) {
    for (int iy=0; iy<KGRID; iy++) {
      if (X[i] > 1.0e-10)
	ans += X[i] * 
	  kernel_x(E, K, Egrid[ix], Kgrid[iy], sigma_E[ix], sigma_K[iy]);
      i++;
    }
  }
  
  return ans;
}

double QPDistF::d2fdE2_EK(double E, double K)
{
  if (!df_computed) compute_distribution();

  double ans = 0.0;
  int i=0;

  for (int ix=0; ix<EGRID; ix++) {
    for (int iy=0; iy<KGRID; iy++) {
      if (X[i] > 1.0e-10)
	ans += X[i] * 
	  kernel_x2(E, K, Egrid[ix], Kgrid[iy], sigma_E[ix], sigma_K[iy]);
      i++;
    }
  }
  
  return ans;
}

double QPDistF::d2fdK2_EK(double E, double K)
{
  if (!df_computed) compute_distribution();

  double ans = 0.0;
  int i=0;

  for (int ix=0; ix<EGRID; ix++) {
    for (int iy=0; iy<KGRID; iy++) {
      if (X[i] > 1.0e-10)
	ans += X[i] * 
	  kernel_y2(E, K, Egrid[ix], Kgrid[iy], sigma_E[ix], sigma_K[iy]);
      i++;
    }
  }
  
  return ans;
}

double QPDistF::d2fdEK_EK(double E, double K)
{
  if (!df_computed) compute_distribution();

  double ans = 0.0;
  int i=0;

  for (int ix=0; ix<EGRID; ix++) {
    for (int iy=0; iy<KGRID; iy++) {
      if (X[i] > 1.0e-10)
	ans += X[i] * 
	  kernel_xy(E, K, Egrid[ix], Kgrid[iy], sigma_E[ix], sigma_K[iy]);
      i++;
    }
  }
  
  return ans;
}

double QPDistF::dfdK_EK(double E, double K)
{
  if (!df_computed) compute_distribution();

  double ans = 0.0;
  int i=0;

  for (int ix=0; ix<EGRID; ix++) {
    for (int iy=0; iy<KGRID; iy++) {
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
  std::ofstream out(name.c_str());
  if (!out) {
    std::cerr << "Couldn't open <" << name << "> to save state!" << std::endl;
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

  out.write((char *)Egrid.data(),   EGRID*sizeof(double));
  out.write((char *)sigma_E.data(), EGRID*sizeof(double));
  out.write((char *)Kgrid.data(),   EGRID*sizeof(double));
  out.write((char *)sigma_K.data(), KGRID*sizeof(double));
  out.write((char *)X.data(),       EGRID*KGRID*sizeof(double));
  out.write((char *)&obj0,          sizeof(double));
  out.write((char *)&obj ,          sizeof(double));

  out.write((char *)&NJMAX , sizeof(int));
  out.write((char *)&Emin, sizeof(double));
  out.write((char *)&Emax, sizeof(double));
  out.write((char *)&TOLE, sizeof(double));
  out.write((char *)JMAXE.data(), NJMAX*sizeof(double));
  out.write((char *)JMAX.data(),  NJMAX*sizeof(double));
  out.write((char *)JMAX2.data(), NJMAX*sizeof(double));
}

// Reinitialize the DF from a saved state

void QPDistF::read_state(string& name)
{
  std::ifstream in(name.c_str());
  if (!in) {
    std::cerr << "Couldn't open <" << name << "> to read state!" << std::endl;
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

  Egrid.resize(EGRID);
  in.read((char *)Egrid.data(), EGRID*sizeof(double));

  sigma_E.resize(EGRID);
  in.read((char *)sigma_E.data(), EGRID*sizeof(double));

  Kgrid.resize(KGRID);
  in.read((char *)Kgrid.data(), KGRID*sizeof(double));

  sigma_K.resize(KGRID);
  in.read((char *)sigma_K.data(), KGRID*sizeof(double));
  
  X.resize(EGRID*KGRID);
  in.read((char *)X.data(), EGRID*KGRID*sizeof(double));

  in.read((char *)&obj0, sizeof(double));
  in.read((char *)&obj , sizeof(double));

  in.read((char *)&NJMAX , sizeof(int));
  in.read((char *)&Emin, sizeof(double));
  in.read((char *)&Emax, sizeof(double));
  in.read((char *)&TOLE, sizeof(double));

  JMAXE.resize(NJMAX);
  JMAX .resize(NJMAX);
  JMAX2.resize(NJMAX);

  in.read((char *)JMAXE.data(), NJMAX*sizeof(double));
  in.read((char *)JMAX.data(),  NJMAX*sizeof(double));
  in.read((char *)JMAX2.data(), NJMAX*sizeof(double));
}

