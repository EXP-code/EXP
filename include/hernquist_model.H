#ifndef _hernquistmodel_H
#define _hernquistmodel_H

#include <cmath>
#include "massmodel.H"

class HernquistSphere : public AxiSymModel
{
private:
  double rot;
  double rmin;
  double rmax;
  double rscl;
  double tolE;

public:

  HernquistSphere(double RSCL = 1.0, 
		   double RMIN = 1.0e-6, double RMAX = 1.0e6, 
		  double TOLE=1.0e-5)
    { rscl = RSCL; rmin = RMIN; rmax = RMAX; 
      dist_defined = true;
      ModelID = "HernquistSphere"; 
      dim = 3;
      tolE = TOLE; }
      

  // Required member functions

  double get_mass(const double r) { 
    return r*r*pow(r+rscl,-2.0);
  }

  double get_density(const double r) {
    if (r>0.0) return rscl/(2.0*M_PI*r)*pow(r+rscl,-3.0);
    else bomb("radius cannot be zero!");
    return 0.0;
  }

  double get_pot(const double r) { 
     return -1.0/(r+rscl);
  }

  double get_dpot(const double r) {
    return 1.0/(r+rscl)/(r+rscl);
  }

  double get_dpot2(const double r) {
    return -2.0/(r+rscl)/(r+rscl)/(r+rscl);
  }
  
  void get_pot_dpot(const double r, double &ur, double &dur) {
    ur = -1.0/(r+rscl); dur = ur*ur;
  }
  
  // Addiional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }

  double distf(double E, double L) {
    double q,cq;

    q = sqrt(-rscl*E);
    cq = sqrt(1.0-q*q);
    return 0.125/(M_SQRT2*M_PI*M_PI*M_PI) /(cq*cq*cq*cq*cq*pow(rscl, 1.5)) *
      ( 3.0*asin(q) + q*cq*(1.0 - 2.0*q*q)*(8.0*q*q*(q*q - 1.0) - 3.0) );
  }

  double dfde(double E, double L) {
    return 0.5*(distf(E+tolE,L)-distf(E-tolE,L))/tolE;
  }

  double dfdl(double E, double L) {
    return 0.0;
  }

  double d2fde2(double E, double L) {
    return (distf(E+tolE,L)+distf(E-tolE,L) - 2*distf(E,L))/(tolE*tolE);
  }


};

#endif
