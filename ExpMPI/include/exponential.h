// This may look like C code, but it is really -*- C++ -*-

//
// Exponential disk
//

const char rcsid_expon[] = "$Id$";

#ifndef _Expon_h

#define _Expon_h 1

#ifndef _Logic_
#include <logic.h>
#endif

#include <massmodel.h>

extern "C" double i0(double);
extern "C" double i1(double);
extern "C" double k0(double);
extern "C" double k1(double);

class ExponentialDisk : public AxiSymModel
{
private:
  
  double a;
  double rmin;
  double rmax;
  double den0;

public:

  ExponentialDisk(double RSCALE=1.0, double RMAX=20.0) { a=RSCALE; 
						   rmin = 1.0e-10;
						   rmax = RMAX;
						   den0 = 0.5/M_PI/a/a;
						   dist_defined = false;
						   dim = 2;
						   ModelID = "ExponentialDisk"; }
      

  // Required member functions

  double get_mass(const double r) { return 1.0 - exp(-r/a)*(1+r/a); }

  double get_density(const double r) { return den0*exp(-r/a); }

  double get_pot(const double r) {
    double y=0.5*r/a;
    return -M_PI*den0*r*(i0(y)*k1(y) - i1(y)*k0(y));
  }
							

  double get_dpot(const double r) {
    double y=0.5*r/a;
    return 2.0*M_PI*den0*y*(i0(y)*k0(y) - i1(y)*k1(y));
  }

  double get_dpot2(const double r) {
    double y=0.5*r/a;
    return M_PI*den0/a*(i0(y)*k0(y) + i1(y)*k1(y) - 
			2.0*y*(i1(y)*k0(y) - i0(y)*k1(y)));
  }

  void get_pot_dpot(const double r, double &ur, double &dur) {
    double y=0.5*r/a;
    double I0=i0(y), I1=i1(y), K0=k0(y), K1=k1(y);
    ur = -M_PI*den0*r*(I0*K1 - I1*K0);
    dur = 2.0*M_PI*den0*y*(I0*K0 - I1*K1);
  }
  
  // Addiional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }

  double distf(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }

  double dfde(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }

  double dfdl(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }
  
  double d2fde2(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }

};



#endif
