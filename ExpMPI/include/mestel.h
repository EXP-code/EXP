// This may look like C code, but it is really -*- C++ -*-

const char rcsid_mestel[] = "$Id$";

#ifndef _mestel_h

#define _mestel_h 1

#ifndef _Logic_
#include <logic.h>
#endif

// class MassModel;
// class DiskModel;

class MestelDisk : public AxiSymModel
{
private:
  
  double rot;
  double rmin;
  double rmax;
  double sig2;
  double q;
  double F;

public:

  MestelDisk(double VROT = 1.0, 
		   double RMIN = 1.0e-6, double RMAX = 1.0e6)
    { rot = VROT*VROT; rmin = RMIN; rmax = RMAX; dist_defined = true;
      dim=2; ModelID = "MestelDisk"; }
      

  // Required member functions

  double get_mass(const double r) { 
    if (r>0.0) return rot*r; 
    else bomb("radius cannot be zero!");
    return 0;
  }

  double get_density(const double r) {
    if (r>0.0) return rot/(2.0*M_PI*r);
    else bomb("radius cannot be zero!");
    return 0;
  }

  double get_pot(const double r) { 
    if (r>0.0) return rot*log(r);
    else bomb("radius cannot be zero!");
    return 0;
  }

  double get_dpot(const double r) {
    if (r>0.0) return rot/r;
    else bomb("radius cannot be zero!");
    return 0;
  }

  double get_dpot2(const double r) {
    if (r>0.0) return -rot/(r*r);
    else bomb("radius cannot be zero!");
    return 0;
  }
  
  void get_pot_dpot(const double r, double &ur, double &dur) {
    if (r>0.0) {ur = rot*log(r); dur = rot/r;}
    else bomb("radius cannot be zero!");
  }
  
  // Addiional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }

  void setup_df(double sigma) { 
    sig2 = sigma*sigma;
    q = rot/sig2 - 1.0;
    F = rot/(4.0*M_PI) / ( sqrt(M_PI) * 
       exp(lgamma(0.5*(q+1.0)) + (2.0 + q)*log(sigma) + 0.5*q*log(2.0)) );
  }

  double distf(double E, double L) {
    L = fabs(L);
    if (L==0.0) return 0.0;
    return F*pow(L, q) * exp(-E/sig2);
  }

  double dfde(double E, double L) {
    L = fabs(L);
    if (L==0.0) return 0.0;
    return -F*pow(L, q) * exp(-E/sig2)/sig2;
  }

  double dfdl(double E, double L) {
    int sgn=1;
    if (L<0) {sgn=-1; L *= sgn;}
    if (L==0.0) return 0.0;
    return q*F*pow(L, q-1.0) * exp(-E/sig2) * sgn;
  }
  
  double d2fde2(double E, double L) {
    L = fabs(L);
    if (L<=0.0) return 0.0;
    return F*pow(L, q) * exp(-E/sig2)/sig2/sig2;
  }

};



#endif
