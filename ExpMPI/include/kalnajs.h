// This may look like C code, but it is really -*- C++ -*-

const char rcsid_kalnajs[] = "$Id$";

#ifndef _kalnajs_h

#define _kalnajs_h 1

#ifndef _Logic_
#include <logic.h>
#endif

// class MassModel;
// class DiskModel;

class KalnajsDisk : public AxiSymModel
{
private:
  
  double rmax;
  double mmax;
  double rhofac;
  double potfac1, potfac2;

  double omega0;
  double omega;
  double fac0, F;

public:

  KalnajsDisk(double RMOD=1.0, double MASS=1.0)
  {
    rmax = RMOD; mmax = MASS; 
    rhofac = mmax/(2.0*rmax*rmax*M_PI/3.0);
    potfac1 = 0.5*M_PI*M_PI*rhofac*rmax;
    potfac2 = 0.5*M_PI*rhofac*rmax;
    dim=2; ModelID = "KalnajsDisk"; 
  }
  


  // Required member functions

  double get_mass(const double r) { 
    if (r>rmax) return mmax;
    else return mmax*(1.0 - pow(1.0 - r*r/(rmax*rmax), 1.5));
  }
      
  double get_density(const double r) {
    if (r>rmax) return 0.0;
    else return rhofac*sqrt((1.0 - r*r/(rmax*rmax)));
  }

  double get_pot(const double r) { 
    double x = r/rmax;
    if (x>1.0) return potfac2*((x*x-2.0)*asin(1.0/x) - sqrt(x*x - 1.0));
    else return potfac1*(0.5*x*x - 1.0);
  }

  double get_dpot(const double r) {
    double x = r/rmax;
    if (x>1.0) return potfac2*2.0*(sqrt(x*x - 1.0)/r - x*asin(1.0/x)/rmax);
    else return potfac1*x/rmax;
  }

  double get_dpot2(const double r) {
    double x = r/rmax;
    if (x>1.0) potfac2*2.0*(2.0/(sqrt(x*x - 1.0)*rmax*rmax) - 
			    asin(1.0/x)/(rmax*rmax) -
			    sqrt(x*x - 1.0)/(r*r) );
    else return potfac1/(rmax*rmax);
  }
  
  void get_pot_dpot(const double r, double &ur, double &dur) {
    double x = r/rmax;
    if (x>1.0)  {ur = potfac2*((x*x-2.0)*asin(1.0/x) - sqrt(x*x - 1.0));
    dur = potfac2*2.0*(sqrt(x*x - 1.0)/r - x*asin(1.0/x)/rmax);}
    else {ur = potfac1*(0.5*x*x - 1.0); dur = potfac1*x/rmax;}
  }
  
  // Addiional member functions

  double get_min_radius(void) { return 0.0; }
  double get_max_radius(void) { return rmax; }

  void setup_df(double frac) { 
    if (fabs(frac)>=1.0) bomb("|frac| must less than 1");
    omega0 = sqrt(0.5*M_PI*M_PI*rhofac/rmax);
    omega  = omega0*frac;
    
    fac0 = (omega0*omega0-omega*omega)*rmax*rmax;
    F = rhofac/(2.0*M_PI*sqrt(fac0));
  }

  double distf(double E, double L) {
    double test = fac0 + 2.0*(E + omega*L);
    if (test>0.0)
      return F/sqrt(test);
    else
      return 0.0;
  }

  double dfde(double E, double L) {
    double test = fac0 + 2.0*(E + omega*L);
    if (test>0.0)
      return -F/(test*sqrt(test));
    else
      return 0.0;
  }

  double dfdl(double E, double L) {
    double test = fac0 + 2.0*(E + omega*L);
    if (test>0.0)
      return -F*omega/(test*sqrt(test));
    else
      return 0.0;
  }

  double d2fde2(double E, double L) {
    double test = fac0 + 2.0*(E + omega*L);
    if (test>0.0)
      return 3.0*F*omega/(test*test*sqrt(test));
    else
      return 0.0;
  }
  
};



#endif
