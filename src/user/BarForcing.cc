// This may look like C code, but it is really -*- C++ -*-

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include <Eigen/Eigen>

#include <biorth.H>
#include <model3d.H>
#include <isothermal.H>
#include <hernquist_model.H>
#include <gaussQ.H>

#include <localmpi.H>
#include <BarForcing.H>


int BarForcing::L0 = 2;
int BarForcing::M0 = 2;
double BarForcing::RMIN = 0.0;
double BarForcing::RMAX = 20.0;
double BarForcing::AMPLITUDE = 0.5;
double BarForcing::LENGTH = 0.5;

BarForcing::BarForcing(int NMAX, double Mass, double Length, 
		       double Corot, double Amp) : Perturbation(NMAX)
{
  ID = "BarForcing";
  user_omega = false;
  mass = Mass;
  length = Length;
  corot = Corot;
  amp = Amp;
  inertia = 1.0;

  LENGTH = length;

  l = 2;
  m = 2;
}


void BarForcing::compute_quad_parameters(double a21, double a32)
{
  const int N = 100;
  LegeQuad gq(N);

  double a1 = length;
  double a2 = a21*a1;
  double a3 = a32*a2;

  double geom = pow(a1*a2*a3, 1.0/3.0);

  double A12 = a1*a1/geom/geom;
  double A22 = a2*a2/geom/geom;
  double A32 = a3*a3/geom/geom;
  
  double u, d, t, denom, ans1=0.0, ans2=0.0;
  for (int i=0; i<N; i++) {
    t = 0.5*M_PI*gq.knot(i);
    u = tan(t);
    d = cos(t);
    d = 1.0/(d*d);
    
    denom = sqrt( (A12+u)*(A22+u)*(A32+u) );
    ans1 += d*gq.weight(i) /( (A12+u)*denom );
    ans2 += d*gq.weight(i) /( (A22+u)*denom );
  }
  ans1 *= 0.5*M_PI;
  ans2 *= 0.5*M_PI;
  
  inertia = 0.2*mass*(a1*a1 + a2*a2);

  
  if (myid==0) {
    cout << "====================================================\n";
    cout << "Computed quadrupole fit to homogenous ellipsoid\n";
    cout << "with Mass=" << mass << " A_1=" << a1 << " A_2=" << a2 
	 << " A_3=" << a3 << "\n"
	 << "with an exact fit to asymptotic quadrupole, e.g.\n"
	 << "     U_{22} = b1 r**2/( 1+(r/b5)**5 ) or\n"
	 << "            = b1 r**2/( 1+ r/b5 )**5\n";
    cout << "====================================================\n";

    cout << "V_1=" << ans1 << endl;
    cout << "V_2=" << ans2 << endl;
    cout << "I_3=" << 0.2*mass*(a1*a1 + a2*a2) << endl;
  }
    
  double rho = mass/(4.0*M_PI/3.0*a1*a2*a3);
  double b25 = 0.4*a1*a2*a3*(a2*a2 - a1*a1)/(ans1 - ans2);
  
  double b1 = M_PI*rho*sqrt(2.0*M_PI/15.0)*(ans1 - ans2);
  double b5 = pow(b25, 0.2);

  // double afac = 2.0 * b1;
  // Single M component
  double afac = b1;
  
  if (myid==0) {
    cout << "b1=" << b1 << endl;
    cout << "b5=" << b5 << endl;
    cout << "afac=" << afac << endl;
    cout << "====================================================\n" 
	 << flush;
  }

  LENGTH = b5;
  AMPLITUDE = afac*amp;
}

AxiSymBiorth *tst;

double testfunc(double r, int l, int m)
{
  return tst->potlR(3, 2, r);
}

double quadfunc(double r, int l, int m)
{
  double x = r/BarForcing::LENGTH;
  return BarForcing::AMPLITUDE * BarForcing::LENGTH*BarForcing::LENGTH
    * pow(x, l)/(1.0 + pow(x, 2.0*l+1.0));
}

double BarForcing::eval(double r)
{
  return quadfunc(r, L0, M0);
}

void BarForcing::compute_coefficients()
{
				// Get inner product (normalized)

  double rmin = max<double>(biorth->rb_to_r(biorth->rb_min()), RMIN);
  double rmax = min<double>(biorth->rb_to_r(biorth->rb_max()), RMAX);

  bcoef = scalar_prod(potential, rmin, rmax,
		      L0, M0, *biorth, quadfunc, nmax, NINT);
}

void BarForcing::compute_omega()
{
  if (!user_omega) 
    omega = sqrt(model->get_dpot(length*corot)/(length*corot));

  omega_computed = true;
}


