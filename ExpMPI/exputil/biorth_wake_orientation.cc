
#include <string>
#include <iostream.h>
#include <math.h>

#include <Vector.h>
#include <CVector.h>
#include <biorth.h>
#include <biorth_wake.h>
#include <simann2.h>

double factrl(int n);
double plgndr(int l, int m, double x);
double rot_matrix(int l, int m, int n, double beta);
Matrix return_euler_slater(double, double, double, int);


void BiorthWake::orientation(int L, int M, 
			     Vector& phi, Vector& theta, Vector& psi,
			     Vector& cost)

{
  ll = L;
  mm = M;

  int n, m, l;

  if (L>lmax) {
    bomb ("wake_orientation: L out of bounds");
    exit(-1);
  }

  if (abs(M)>L) {
    bomb ("wake_orientation: M out of bounds");
    exit(-1);
  }

  phi.setsize(1, nmax);
  theta.setsize(1, nmax);
  psi.setsize(1, nmax);
  cost.setsize(1, nmax);

  ylm.setsize(-L, L);

				// For debugging only
#ifdef DEBUG
  test_transform();
#endif // DEBUG

  int loffset=0, moffset;
  for (l=0; l<L; l++) loffset += 2*l+1;

  double fac1, fac2, signflip=1.0, norm;

  for (n=1; n<=nmax; n++) {

    for (m=0, moffset=0; m<=l; m++) {

      fac1 = sqrt( (0.5*l+0.25)/M_PI );

      if (m==0) {
	ylm[0] = fac1 * expcoef[loffset+moffset][n];
	moffset++;
      }
      else {
	fac2 = fac1 * sqrt( factrl(l-m)/factrl(l+m) );
	ylm[-m] = fac2 * signflip *
	  Complex(expcoef[loffset+moffset][n], -expcoef[loffset+moffset+1][n]);
	ylm[ m] = fac2 * 
	  Complex(expcoef[loffset+moffset][n],  expcoef[loffset+moffset+1][n]);

	moffset +=2;
      }
      signflip *= -1.0;

    }

    norm = fabs(sqrt(ylm.Conjg()*ylm)) + 1.0e-10;

    ylm /= norm;

    get_transform(phi[n], theta[n], psi[n], cost[n]);

    cost[n] *= norm;

#ifdef DEBUG    
    check_orientation(phi[n], theta[n], psi[n]);
#endif      

  }

}


void BiorthWake::orientation_init(void)
{
  param = new double [ndim];
  psum = new double [ndim] - 1;
  ptry = new double [ndim] - 1;
  ambp = new double* [ndim+1] - 1;
  amby = new double  [ndim+1] - 1;

  for (int i=1; i<=ndim+1; i++) ambp[i] = new double [ndim] - 1;

  init_orientation = true;
}

BiorthWake* current;

				// World function accessible by SA
double sa_energy(double *params)
{
  return current->energy(params);
}

double BiorthWake::energy(double *params)
{
  Complex ansp=0.0, ansm=0.0;
  int n;

  for (n=-ll; n<=ll; n++) {
    ansp += exp(I*params[2]*n*(-1)) * exp(I*params[0]*mm*(-1)) * ylm[n]
      * rot_matrix(ll, mm, n, params[1]);  
  }
  
  if (mm != 0) {
    for (n=-ll; n<=ll; n++) {
      ansm += exp(I*params[2]*n*(-1)) * exp(I*params[0]*mm) * ylm[n]
	* rot_matrix(ll, -mm, n, params[1]);  
    }
  }
  
  return -(ansp.real()*ansp.real() + ansm.real()*ansm.real());
}

double BiorthWake::amoeba_energy(double *params)
{
  return energy(params+1);
}

void BiorthWake::modulo_param(double *params)
{
  int i, indx;

  params[1] += M_PI;

  for (i=0; i<3; i++) {

    indx = (int)(0.5*params[i]/M_PI);

    if (params[i]>=0.0)
      params[i] -= 2.0*M_PI*indx;
    else
      params[i] += 2.0*M_PI*(1-indx);
  }
  
  params[1] -= M_PI;

}


void BiorthWake::get_transform(double& phi, double& theta, double& psi, 
			       double& cost)
{
  if (!init_orientation) orientation_init();

  int i;

				// Initialize SimAnneal

				// Initial values
  for (i=0; i<=ndim; i++) param[i] = 0.5*M_PI;

  current = this;

  SimAnneal sa(sa_energy, ndim);

  if ( !sa ) {
    cerr << "problem initializing SimAnneal object\n";
    exit(1);
  }

  sa.melt();
  sa.anneal(iter);
  sa.optimum(param);
  modulo_param(param);


  for (i=1; i<=ndim; i++) 
    ambp[1][i] = ambp[2][i] = ambp[3][i] = ambp[4][i] = param[i-1];
  for (i=1; i<=ndim; i++) ambp[i][i] += 1.0e-2;
  for (i=1; i<=ndim+1; i++) amby[i] = energy( ambp[i] );

  amoeba();

  int which=1;
  double zmin = amby[1];
  for (i=2; i<=ndim+1; i++) 
    if (amby[i] < zmin) {
      which = i;
      zmin = amby[i];
    }

  phi = ambp[which][1];
  theta = ambp[which][2];
  psi = ambp[which][3];
  cost = zmin;

}


#define NMAX 5000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0
  
#define GET_PSUM for (j=1;j<=ndim;j++) { for (i=1,sum=0.0;i<=mpts;i++)\
					   sum += ambp[i][j]; psum[j]=sum;}


double BiorthWake::amotry(int ihi, double fac)
{
  int j;
  double fac1, fac2, ytry;

  fac1 = (1.0-fac)/ndim;
  fac2 = fac1-fac;
  for (j=1;j<=ndim;j++) ptry[j] = psum[j]*fac1-ambp[ihi][j]*fac2;
  ytry = amoeba_energy(ptry);
  ++nfunk;
  if (ytry < amby[ihi]) {
    amby[ihi] = ytry;
    for (j=1;j<=ndim;j++) {
      psum[j] += ptry[j]-ambp[ihi][j];
      ambp[ihi][j] = ptry[j];
    }
  }
  return ytry;
}

void BiorthWake::amoeba(void)
{
  int i, j, ilo, ihi, inhi, mpts=ndim+1;
  double ytry, ysave, sum, rtol;
  
  nfunk = 0;
  
  GET_PSUM
    for (;;) {
      ilo = 1;
      ihi = amby[1]>amby[2] ? (inhi=2,1) : (inhi=1,2);
      for (i=1;i<=mpts;i++) {
	if (amby[i] < amby[ilo]) ilo = i;
	if (amby[i] > amby[ihi]) {
	  inhi = ihi;
	  ihi = i;
	} else if (amby[i] > amby[inhi])
	  if (i != ihi) inhi = i;
      }
      rtol = 2.0*fabs(amby[ihi]-amby[ilo])/(fabs(amby[ihi])+fabs(amby[ilo]));
      if (rtol < tol) break;
      if (nfunk >= NMAX) break;

      ytry = amotry(ihi,-ALPHA);
      if (ytry <= amby[ilo])
	ytry = amotry(ihi,GAMMA);
      else if (ytry >= amby[inhi]) {
	ysave = amby[ihi];
	ytry = amotry(ihi,BETA);
	if (ytry >= ysave) {
	  for (i=1;i<=mpts;i++) {
	    if (i != ilo) {
	      for (j=1;j<=ndim;j++) {
		psum[j] = 0.5*(ambp[i][j]+ambp[ilo][j]);
		ambp[i][j] = psum[j];
	      }
	      amby[i] = amoeba_energy(psum);
	    }
	  }
	  nfunk += ndim;
	  GET_PSUM
	  }
      }
    }
}


#undef ALPHA
#undef BETA
#undef GAMMA
#undef NMAX

  
#ifdef DEBUG

#include <gaussQ.h>
#include <iomanip.h>

Complex BiorthWake::test_fct(double theta, double phi)
{
  return 
    sqrt( (0.5*ll + 0.25)/M_PI *
	  exp(lgamma(1.0+ll-mm) - lgamma(1.0+ll+mm)) ) * 
    plgndr(ll, mm, cos(theta)) * exp(I*phi*mm);
}

void BiorthWake::test_transform(void)
{
  const int NINT=20;

  double PHI=33;
  double THETA=66;
  double PSI=111;

  double onedeg = M_PI/180.0;

  PHI *= onedeg;
  THETA *= onedeg;
  PSI *= onedeg;

  Matrix trans = return_euler_slater(PHI, THETA, PSI, 1);

  Vector x0(1,3), x1(1,3);
  ylm.zero();
  
				// 2-d integral over theta and phi
  LegeQuad wk(NINT);

  int it, ip, m, n;
  double cost, sint, theta, phi, theta1, phi1, signflip, psi;
  Complex fac, fac2;

  for (it=1; it<=NINT; it++) {
    
    cost = 2.0*(wk.knot(it) - 0.5);
    sint = sqrt(1.0 - cost*cost);
    theta = acos(cost);

    for (ip=1; ip<=NINT; ip++) {
    
      phi = 2.0*M_PI*wk.knot(ip);

      x0[1] = sint*cos(phi);
      x0[2] = sint*sin(phi);
      x0[3] = cost;

      x1 = trans * x0;

      theta1 = acos(x1[3]/sqrt(x1*x1));
      phi1 = atan2(x1[2], x1[1]);

      fac = 4.0*M_PI * wk.weight(it) * wk.weight(ip) * test_fct(theta1, phi1);

      signflip = 1.0;

      for (m=0; m<=ll; m++) {
	
	fac2 = sqrt( (0.5*ll + 0.25)/M_PI *
		     exp(lgamma(1.0+ll-m) - lgamma(1.0+ll+m))
		     ) * plgndr(ll, m, cos(theta)) * exp(I*phi*m*(-1));

	ylm[m] += fac * fac2;
	if (m) ylm[-m] += fac * conjg(fac2) * signflip;

	signflip *= -1.0;
      }
    }
  }

  CMatrix rot(-ll, ll, -ll, ll);

  for (m=-ll; m<=ll; m++) {
    for (n=-ll; n<=ll; n++) {
      rot[m][n] = rot_matrix(ll, m, n, THETA) * 
	exp(I*PSI*n*(-1)) * exp(I*PHI*m*(-1));
    }
  }

  CVector ylm2 = rot * ylm;

  cout.precision(6);
  cout.setf(ios::scientific);

  for (m=-ll; m<=ll; m++)
    cout << setw(5) << m 
	 << setw(5) << ylm[m]
	 << setw(5) << ylm2[m] << endl;

  get_transform(phi, theta, psi, cost);

  cout << endl;
  cout << setw(15) << "phi" << setw(15) << PHI/onedeg 
       << setw(15) << phi/onedeg << endl;
  cout << setw(15) << "theta" << setw(15) << THETA/onedeg 
       << setw(15) << theta/onedeg << endl;
  cout << setw(15) << "psi" << setw(15) << PSI/onedeg 
       << setw(15) << psi/onedeg << endl;
  cout << setw(15) << "cost" << setw(15) << cost;
  cout << endl;

}

void BiorthWake::check_orientation(double phi, double theta, double psi)
{
  CMatrix rot(-ll, ll, -ll, ll);

  int m, n;

  for (m=-ll; m<=ll; m++) {
    for (n=-ll; n<=ll; n++) {
      rot[m][n] = rot_matrix(ll, m, n, theta) * 
	exp(I*psi*n*(-1)) * exp(I*phi*m*(-1));
    }
  }

  CVector ylm2 = rot * ylm;

  cout.precision(6);

  for (m=-ll; m<=ll; m++)
    cout << setw(5) << m 
	 << setw(5) << ylm[m]
	 << setw(5) << ylm2[m] << endl;

  cout << endl;

}

#endif // DEBUG
