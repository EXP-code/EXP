//
// Toomre's Model 1
// with Kalnajs' m distribution function derived from inverting the
// defining integral equation with a separable form: e^m g(x)
//

#ifndef _toomre_H
#define _toomre_H

#include "massmodel.H"

class ToomreDisk : public AxiSymModel
{
private:
  
  double rmax;
  double m, e, x;
  double p0, p1, p2;

  void pdist(double E, double L);

public:

  inline ToomreDisk(int M=0, double RMAX=1.0e6) { 
    m = M+1.0e-10;
    rmax = RMAX;
    dist_defined = true;
    dim = 2;
    ModelID = "ToomreDisk"; }
      

  // Required member functions

  double get_mass(const double r) { return 1.0 - 1.0/sqrt(1.0 + r*r);}

  double get_density(const double r) { return pow(1.0 + r*r, -1.5)/(2.0*M_PI);}

  double get_pot(const double r) { return -1.0/sqrt(1.0 + r*r); }

  double get_dpot(const double r) { return r*pow(1.0 + r*r, -1.5); }

  double get_dpot2(const double r) { return (1.0 - 2.0*r*r)*
				       pow(1.0 + r*r, -2.5); }

  void get_pot_dpot(const double r, double &ur, double &dur) {
    ur = 1.0/sqrt(1.0 + r*r);
    dur = r*ur*ur*ur;
    ur *= -1.0;
  }
  
  // Addiional member functions

  double get_min_radius(void) { return 0.0; }
  double get_max_radius(void) { return rmax; }

  double distf(double E, double L);

  double dfde(double E, double L);

  double dfdl(double E, double L);
  
  double d2fde2(double E, double L);

};

#endif
