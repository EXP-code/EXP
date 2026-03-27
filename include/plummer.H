#ifndef _plummer_H
#define _plummer_H

#include <cmath>
#include "massmodel.H"

class PlummerSphere : public AxiSymModel
{
  using AxiSymModel::get_mass;
  using AxiSymModel::get_density;
  using AxiSymModel::get_pot;

private:

  double rmin;
  double rmax;
  double rscl;
  double mass;

  double dfac;
  double ffac;

public:

  PlummerSphere(double RSCL = 1.0, double RMIN = 1.0e-6, double RMAX = 1.0e6,
		double MASS = 1.0)
  {
    rscl = RSCL; rmin = RMIN; rmax = RMAX; mass = MASS;
    dist_defined = true;
    ModelID = "PlummerSphere"; 
    dim = 3;

    dfac = 3.0*mass*rscl*rscl/(4.0*M_PI);
    ffac = 128.0/7.0 / (sqrt(8.0)*M_PI*M_PI) * 3.0*rscl*rscl/
      ( 4.0*M_PI*pow(mass, 4.0) );
  }

  // Required member functions

  double get_mass(const double r) { 
    double q = r/rscl;
    return mass*q*q*q/pow(1 + q*q, 1.5);
  }

  double get_density(const double r) {
    return dfac * pow( rscl*rscl + r*r, -2.5 );
  }

  double get_pot(const double r) { 
    return -mass / sqrt( rscl*rscl + r*r ); 
  }

  double get_dpot(const double r) {
    return mass*r * pow( rscl*rscl + r*r, -1.5 );
  }

  double get_dpot2(const double r) {
    return mass *  pow( rscl*rscl + r*r, -2.5 ) * ( rscl*rscl - 2.0*r*r );
  }
  
  void get_pot_dpot(const double r, double &ur, double &dur) {
    ur = -mass / sqrt( rscl*rscl + r*r ); 
    dur = -ur * r / ( rscl*rscl + r*r );
  }
  
  // Addiional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }

  double distf(double E, double L) {
    return ffac * pow(-E, 3.5);
  }

  double dfde(double E, double L) {
    return -3.5*ffac*pow(-E, 2.5);
  }

  double dfdl(double E, double L) {
    return 0.0;
  }

  double d2fde2(double E, double L) {
    return 3.5*2.5*ffac*pow(-E, 1.5);
  }

};

#endif


