#include <filesystem>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <limits>
#include <string>
#include <cmath>

#include "EmpDeproj.H"
#include "gaussQ.H"

// EmpCylSL: Empirical cylindrical deprojection by numerical
// integration and finite difference

EmpDeproj::EmpDeproj(double H, double RMIN, double RMAX, int NUMR, int NINT,
		     std::function<double(double, double)> func, AbelType type)
{
  LegeQuad lq(NINT);

  std::vector<double> rr(NUMR), rl(NUMR), sigI(NUMR), rhoI(NUMR, 0.0);

  double Rmin = log(RMIN);
  double Rmax = log(RMAX);

  double dr = (Rmax - Rmin)/(NUMR-1);

  // Compute surface mass density, Sigma(R)
  //
  for (int i=0; i<NUMR; i++) {
    double r = Rmin + dr*i;

    // Save for finite difference
    //
    rl[i] = r;
    r = exp(r);
    rr[i] = r;

    // Interval by Legendre
    //
    sigI[i] = 0.0;
    for (int n=0; n<NINT; n++) {
      double y   = lq.knot(n);
      double y12 = 1.0 - y*y;
      double z   = y/sqrt(y12)*H;

      sigI[i] += lq.weight(n)*2.0*H*pow(y12, -1.5)*func(r, z);
    }
  }

  Linear1d surf(rl, sigI);
  surfRg = surf;
  
  // Now, compute Abel inverion integral
  //
  for (int i=0; i<NUMR; i++) {

    double surfS = type==AbelType::Subtraction ? surf.eval(rl[i]) : 0.0;
    double r = rr[i];

    // Interval by Legendre
    //
    rhoI[i] = 0.0;

    for (int n=0; n<NINT; n++) {

      double x   = lq.knot(n);
      double x12 = 1.0 - x*x;
      double R   = r/sqrt(x12);
      double lR  = log(R);

      switch (type) {
      case AbelType::Derivative:
	rhoI[i]   += lq.weight(n)/x12 * surf.deriv(lR)/R;
	break;
      case AbelType::Subtraction:
	rhoI[i]   += lq.weight(n)/(x*x*sqrt(x12)*r) * ( surf.eval(lR) - surfS );
	break;
      case AbelType::IBP:
	rhoI[i]   += lq.weight(n)/x12 * surf.eval(lR);
	break;
      }
    }
  }
  
  Linear1d intgr;
  if (type == AbelType::IBP) intgr = Linear1d(rl, rhoI);

  std::vector<double> rho(NUMR), mass(NUMR);
  for (int i=0; i<NUMR; i++) {
    if (type == AbelType::IBP)
      rho[i] = -intgr.deriv(rl[i])/rr[i]/M_PI;
    else
      rho[i] = -rhoI[i]/M_PI;
  }

  mass[0] = 0.0;
  for (int i=1; i<NUMR; i++) {
    double rlst = rr[i-1], rcur = rr[i];
    mass[i] = mass[i-1] + 2.0*M_PI*(rlst*rlst*rho[i-1] + rcur*rcur*rho[i])*(rcur - rlst);
  }

  // Debug
  //
  if (true) {
    std::string fname("deproject_sl.txt");
    std::ofstream out(fname);
    if (out) {
      for (int i=0; i<NUMR; i++)
	out << std::setw(18) << rl[i]
	    << std::setw(18) << rr[i]
	    << std::setw(18) << rho[i]
	    << std::setw(18) << mass[i]
	    << std::endl;
    }
  }

  // Finalize
  //
  densRg = Linear1d(rl, rho);
  massRg = Linear1d(rl, mass);
}
