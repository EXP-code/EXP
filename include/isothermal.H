#ifndef _isothermal_H
#define _isothermal_H

#include <cmath>
#include "massmodel.H"
#include "interp.H"

class SingIsothermalSphere : public AxiSymModel
{
private:
  double rot;
  double rmin;
  double rmax;
  double F;

public:

  //! Constructor
  SingIsothermalSphere(double VROT = 1.0, 
		       double RMIN = 1.0e-6, double RMAX = 1.0e6)
  {
    rot = VROT*VROT; rmin = RMIN; rmax = RMAX; 
    dist_defined = true;
    ModelID = "SingIsothermalSphere"; 
    dim = 3;
    F = 1.79587122125166561689e-1/(VROT*VROT*VROT);
  }
  
  
  //@{
  //! Required member functions

  double get_mass(const double r) { 
    if (r>0.0) return rot*r; 
    else bomb("radius cannot be zero!");
    return 0.0;
  }

  double get_density(const double r) {
    if (r>0.0) return rot/(4.0*M_PI*r*r);
    else bomb("radius cannot be zero!");
    return 0.0;
  }

  double get_pot(const double r) { 
    if (r>0.0) return rot*log(r);
    else bomb("radius cannot be zero!");
    return 0.0;
  }

  double get_dpot(const double r) {
    if (r>0.0) return rot/r;
    else bomb("radius cannot be zero!");
    return 0.0;
  }

  double get_dpot2(const double r) {
    if (r>0.0) return -rot/(r*r);
    else bomb("radius cannot be zero!");
    return 0.0;
  }
  
  void get_pot_dpot(const double r, double &ur, double &dur) {
    if (r>0.0) {ur = rot*log(r); dur = rot/r;}
    else bomb("radius cannot be zero!");
  }
  //@}
  
  //@{
  //! Addiional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }

  double distf(double E, double L) {
    return F*exp(-2.0*E/rot);
  }

  double dfde(double E, double L) {
    return -2.0*F*exp(-2.0*E/rot)/rot;
  }

  double dfdl(double E, double L) {
    return 0.0;
  }

  double d2fde2(double E, double L) {
    return 4.0*F*exp(-2.0*E/rot)/(rot*rot);
  }
  //@}

};



class IsothermalSphere : public AxiSymModel
{
private:
  double vrot;
  double sigma;
  double rscale;
  double rmax;
  double F;

  int num;
  Eigen::VectorXd r, d, m, p;
  Eigen::VectorXd d2, m2, p2;

  double ret1, ret2;

public:

  IsothermalSphere(double RCORE = 1.0, double RMAX = 20.0, double VROT = 1.0, 
		   int NUM=2000, int dN=20);

  // Required member functions

  double get_mass(const double R) { 
    if (R<0.0) bomb("radius cannot be less than zero!");
    Splint1(r, m, m2, R, ret1, 1);
    return ret1;
  }

  double get_density(const double R) {
    if (R<0.0) bomb("radius cannot be less than zero!"); 
    Splint1(r, d, d2, R, ret1, 1);
   return ret1;
  }

  double get_pot(const double R) { 
    if (R<0.0) bomb("radius cannot be less than zero!");
    Splint1(r, p, p2, R, ret1, 1);
    return ret1;
  }

  double get_dpot(const double R) {
    if (R<0.0) bomb("radius cannot be less than zero!");
    if (R==0.0) return 0.0;
    Splint1(r, m, m2, R, ret1, 1);
    return ret1/(R*R);
  }

  double get_dpot2(const double R) {
    if (R<0.0) bomb("radius cannot be less than zero!");
    if (R==0.0) return 3.0*sigma*sigma/rscale;
    Splint1(r, m, m2, R, ret1, 1);
    Splint1(r, d, d2, R, ret2, 1);
    return -2.0*ret1/(R*R*R) + 4.0*M_PI*ret2;
  }
  
  void get_pot_dpot(const double R, double &ur, double &dur) {
    if (R<0.0) bomb("radius cannot be less than zero!");
    if (R==0)
      ur = dur = 0.0;
    else {
      Splint1(r, p, p2, R, ur, 1);
      Splint1(r, m, m2, R, dur, 1);
      dur *= 1.0/(R*R);
    }
  }
  
  // Addiional member functions

  double get_min_radius(void) { return 0.0; }
  double get_max_radius(void) { return rmax; }

  double distf(double E, double L) {
    return F*exp(-E/(sigma*sigma));
  }

  double dfde(double E, double L) {
    return -F*exp(-E/sigma*sigma)/(sigma*sigma);
  }

  double dfdl(double E, double L) {
    return 0.0;
  }

  double d2fde2(double E, double L) {
    return F*exp(-E/(sigma*sigma))/(sigma*sigma*sigma*sigma);
  }


};


class LowSingIsothermalSphere : public AxiSymModel
{
private:
  double rmin;
  double rmax;
  double F, Fconst, sigma;

  int num;
  Eigen::VectorXd u, d, m, p;
  Eigen::VectorXd d2, m2, p2;

  double ret1, ret2;

public:

  LowSingIsothermalSphere(double RMIN = 1.0e-3, double RMAX = 20.0, 
			  int NUM=2000, int dN=20);

  // Required member functions

  double get_mass(const double R) { 
    if (R<=0.0) bomb("radius cannot be zero!");
    if (R<=rmin) return m[1]*R/rmin;
    if (R>=rmax) return m[num];
    Splint1(u, m, m2, log(R), ret1, 1);
    return ret1;
  }

  double get_density(const double R) {
    if (R<=0.0) bomb("radius cannot be zero!"); 
    if (R<=rmin) return d[1]*rmin*rmin/(R*R);
    if (R>=rmax) return 0.0;
    Splint1(u, d, d2, log(R), ret1, 1);
   return ret1;
  }

  double get_pot(const double R) { 
    if (R<=0.0) bomb("radius cannot be zero!");
    if (R<=rmin) return u[1]*log(R)/log(rmin);
    if (R>=rmax) return -m[num]/R;
    Splint1(u, p, p2, log(R), ret1, 1);
    return ret1;
  }

  double get_dpot(const double R) {
    if (R<=0.0) bomb("radius cannot be zero!");
    if (R<=rmin) return u[1]/(R*log(rmin));
    if (R>=rmax) return m[num]/(R*R);
    Splint1(u, m, m2, log(R), ret1, 1);
    return ret1/(R*R);
  }

  double get_dpot2(const double R) {
    if (R<0.0) bomb("radius cannot be zero!");
    if (R<=rmin) return -u[1]/(R*R*log(rmin));
    if (R>=rmax) return -2.0*m[num]/(R*R*R);
    Splint1(u, m, m2, log(R), ret1, 1);
    Splint1(u, d, d2, log(R), ret2, 1);
    return -2.0*ret1/(R*R) + 4.0*M_PI*ret2;
  }
  
  void get_pot_dpot(const double R, double &ur, double &dur) {
    if (R<0.0) bomb("radius cannot be zero!");
    if (R<=rmin) {
      ur = u[1]*log(R)/log(rmin);
      dur = u[1]/(R*log(rmin));
    } else if (R>=rmax) {
      ur = -m[num]/R;
      dur = m[num]/(R*R);
    } else {
      Splint1(u, p, p2, log(R), ur, 1);
      Splint1(u, m, m2, log(R), dur, 1);
      dur *= 1.0/(R*R);
    }
  }
  
  // Addiional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }

  double distf(double E, double L) {
    if (E>p[num]) return 0.0;
    return F*exp(-E/(sigma*sigma)) - Fconst;
  }

  double dfde(double E, double L) {
    if (E>p[num]) return 0.0;
    return -F*exp(-E/sigma*sigma)/(sigma*sigma);
  }

  double dfdl(double E, double L) {
    return 0.0;
  }

  double d2fde2(double E, double L) {
    if (E>p[num]) return 0.0;
    return F*exp(-E/(sigma*sigma))/(sigma*sigma*sigma*sigma);
  }


};


#endif
