//
// Exponential disk
//

#ifndef _Expon_H
#define _Expon_H

#include "massmodel.H"

class ExponentialDisk : public AxiSymModel
{
private:
  
  double a;
  double m;
  double rmin;
  double rmax;
  double den0;

public:

  ExponentialDisk(double RSCALE=1.0, double RMAX=20.0, double DMASS=1.0) 
  { 
    a       = RSCALE; 
    rmin    = 1.0e-8;
    rmax    = RMAX;
    m       = DMASS;
    den0    = m*0.5/M_PI/a/a;
    dim     = 2;
    ModelID = "ExponentialDisk"; 

    dist_defined = false;
  }
      

  // Required member functions

  double get_mass(const double r) { return m*(1.0 - exp(-r/a)*(1.0+r/a)); }

  double get_density(const double r) { return den0*exp(-r/a); }

  double get_pot(const double r);

  double get_dpot(const double r);

  double get_dpot2(const double r);

  void get_pot_dpot(const double r, double &ur, double &dur);

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
