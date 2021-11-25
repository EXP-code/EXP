
#include <string>
#include <iostream>
#include <cmath>

#include <biorth.H>
#include <biorth_wake.H>
#include <simann2.h>

double factrl(int n);
double plgndr(int l, int m, double x);
double rot_matrix(int l, int m, int n, double beta);
Eigen::Matrix3d return_euler_slater(double, double, double, int);

using namespace std;

void BiorthWake::orientation(int L, int M, 
			     Eigen::VectorXd& phi, Eigen::VectorXd& theta, Eigen::VectorXd& psi,
			     Eigen::VectorXd& cost)

{
  ll = L;
  mm = M;

  if (L>lmax) {
    bomb ("wake_orientation: L out of bounds");
    exit(-1);
  }

  if (abs(M)>L) {
    bomb ("wake_orientation: M out of bounds");
    exit(-1);
  }

  phi.resize(nmax);
  theta.resize(nmax);
  psi.resize(nmax);
  cost.resize(nmax);

  ylm.resize(2*L+1);

				// For debugging only
#ifdef DEBUG
  test_transform();
#endif // DEBUG

  int loffset=0, moffset;
  int l=0;
  for (; l<L; l++) loffset += 2*l+1;

  double fac1, fac2, signflip=1.0, norm;

  for (int n=0; n<nmax; n++) {

    for (int m=0, moffset=0; m<=l; m++) {

      fac1 = sqrt( (0.5*l+0.25)/M_PI );

      if (m==0) {
	ylm[L] = fac1 * expcoef(loffset+moffset, n);
	moffset++;
      }
      else {
	fac2 = fac1 * sqrt( factrl(l-m)/factrl(l+m) );
	ylm[L-m] = fac2 * signflip *
	  std::complex<double>(expcoef(loffset+moffset, n), -expcoef(loffset+moffset+1, n));
	ylm[L+m] = fac2 * 
	  std::complex<double>(expcoef(loffset+moffset, n),  expcoef(loffset+moffset+1, n));

	moffset +=2;
      }
      signflip *= -1.0;

    }

    auto conj = ylm.conjugate();

    norm = sqrt(fabs(conj.dot(ylm))) + 1.0e-10;

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
  param .resize(ndim);
  psum  .resize(ndim);
  ptry  .resize(ndim);
  ambp  .resize(ndim+1);
  amby  .resize(ndim+1);

  for (auto & v : ambp) v.resize(ndim);

  init_orientation = true;
}

double BiorthWake::energy(std::vector<double>& params)
{
  std::complex<double> ansp=0.0, ansm=0.0;

  for (int n=-ll; n<=ll; n++) {
    std::complex<double> cn = -n, cmm = -mm;
    ansp += exp(I*params[2]*cn) * exp(I*params[0]*cmm) * ylm[ll+n]
      * rot_matrix(ll, mm, n, params[1]);  
  }
  
  if (mm != 0) {
    for (int n=-ll; n<=ll; n++) {
      std::complex<double> cn = -n, cmm = mm;
      ansm += exp(I*params[2]*cn) * exp(I*params[0]*cmm) * ylm[ll+n]
	* rot_matrix(ll, -mm, n, params[1]);  
    }
  }
  
  return -(ansp.real()*ansp.real() + ansm.real()*ansm.real());
}

double BiorthWake::amoeba_energy(std::vector<double>& params)
{
  return energy(params);
}

void BiorthWake::modulo_param(std::vector<double>& params)
{
  params[1] += M_PI;

  for (int i=0; i<3; i++) {

    int indx = (int)(0.5*params[i]/M_PI);

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

				// Initialize SimAnneal

				// Initial values
  for (int i=0; i<ndim; i++) param[i] = 0.5*M_PI;

  auto F = [this](std::vector<double>& p) {return this->energy(p);};

  SimAnneal sa(F, ndim);

  if ( !sa ) {
    cerr << "problem initializing SimAnneal object\n";
    exit(1);
  }

  sa.melt();
  sa.anneal(iter);
  sa.optimum(param);
  modulo_param(param);


  for (int i=0; i<ndim; i++) 
    ambp[0][i] = ambp[1][i] = ambp[2][i] = ambp[3][i] = param[i];
  for (int i=0; i<ndim; i++)   ambp[i][i] += 1.0e-2;
  for (int i=0; i<ndim+1; i++) amby[i] = energy(ambp[i]);

  amoeba();

  int which=0;
  double zmin = amby[0];
  for (int i=1; i<ndim+1; i++) 
    if (amby[i] < zmin) {
      which = i;
      zmin = amby[i];
    }

  phi   = ambp[which][0];
  theta = ambp[which][1];
  psi   = ambp[which][2];
  cost  = zmin;

}

#define NMAX 5000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0
  
#define GET_PSUM for (int j=0;j<ndim;j++) { \
    for (int i=0, sum=0.0; i<mpts;i++)	    \
      sum += ambp[i][j]; psum[j]=sum;}


double BiorthWake::amotry(int ihi, double fac)
{
  double fac1 = (1.0-fac)/ndim;
  double fac2 = fac1-fac;
  for (int j=0; j<ndim; j++) ptry[j] = psum[j]*fac1-ambp[ihi][j]*fac2;
  double ytry = amoeba_energy(ptry);
  ++nfunk;
  if (ytry < amby[ihi]) {
    amby[ihi] = ytry;
    for (int j=0; j<ndim; j++) {
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
      ihi = amby[0]>amby[1] ? (inhi=1,0) : (inhi=0,1);
      for (int i=0; i<mpts; i++) {
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
	  for (int i=0 ; i<mpts; i++) {
	    if (i != ilo) {
	      for (int j=1; j<ndim; j++) {
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

#include <gaussQ.H>
#include <iomanip>

std::complex<double> BiorthWake::test_fct(double theta, double phi)
{
  return 
    sqrt( (0.5*ll + 0.25)/M_PI *
	  exp(lgamma(1.0+ll-mm) - lgamma(1.0+ll+mm)) ) * 
    plgndr(ll, mm, cos(theta)) * exp(I*phi*static_cast<double>(mm));
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

  Eigen::Matrix3d trans = return_euler_slater(PHI, THETA, PSI, 1);

  Eigen::Vector3d x0, x1;
  ylm.setZero();
  
				// 2-d integral over theta and phi
  LegeQuad wk(NINT);

  int m, n;
  double cost, sint, theta, phi, theta1, phi1, signflip, psi;
  std::complex<double> fac, fac2;

  for (int it=1; it<NINT; it++) {
    
    cost = 2.0*(wk.knot(it) - 0.5);
    sint = sqrt(1.0 - cost*cost);
    theta = acos(cost);

    for (int ip=0; ip<NINT; ip++) {
    
      phi = 2.0*M_PI*wk.knot(ip);

      x0[0] = sint*cos(phi);
      x0[1] = sint*sin(phi);
      x0[2] = cost;

      x1 = trans * x0;

      theta1 = acos(x1[2]/sqrt(x1.dot(x1)));
      phi1 = atan2(x1[1], x1[0]);

      fac = 4.0*M_PI * wk.weight(it) * wk.weight(ip) * test_fct(theta1, phi1);

      signflip = 1.0;

      for (m=0; m<=ll; m++) {
	
	fac2 = sqrt( (0.5*ll + 0.25)/M_PI *
		     exp(lgamma(1.0+ll-m) - lgamma(1.0+ll+m))
		     ) * plgndr(ll, m, cos(theta)) * exp(I*phi*static_cast<double>(-m));

	ylm[ll+m] += fac * fac2;
	if (m) ylm[ll-m] += fac * std::conj(fac2) * signflip;

	signflip *= -1.0;
      }
    }
  }

  Eigen::MatrixXcd rot(2*ll+1, 2*ll+1);

  for (m=-ll; m<=ll; m++) {
    for (n=-ll; n<=ll; n++) {
      rot(ll+m, ll+n) = rot_matrix(ll, m, n, THETA) * 
	exp(I*PSI*static_cast<double>(-n)) *
	exp(I*PHI*static_cast<double>(-m)) ;
    }
  }

  auto ylm2 = rot * ylm;

  cout.precision(6);
  cout.setf(ios::scientific);

  for (m=-ll; m<=ll; m++)
    cout << setw(5) << m 
	 << setw(5) << ylm[ll+m]
	 << setw(5) << ylm2[ll+m] << endl;

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
  Eigen::MatrixXcd rot(2*ll+1, 2*ll+1);

  for (int m=-ll; m<=ll; m++) {
    for (int n=-ll; n<=ll; n++) {
      rot(m+ll, n+ll) = rot_matrix(ll, m, n, theta) * 
	exp(I*psi*static_cast<double>(-n)) *
	exp(I*phi*static_cast<double>(-m)) ;
    }
  }

  auto ylm2 = rot * ylm;

  std::cout.precision(6);

  for (int m=-ll; m<=ll; m++)
    std::cout << setw(5) << m 
	 << setw(5) << ylm[ll+m]
	      << setw(5) << ylm2[ll+m] << std::endl;

  std::cout << std::endl;

}

#endif // DEBUG
