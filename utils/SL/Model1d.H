#ifndef _Model1d_H_
#define _Model1d_H_

#include <string>

#include <interp.H>

class OneDModel
{

public:
  bool dist_defined;

  virtual double get_mass(const double) = 0;
  virtual double get_density(const double) = 0;
  virtual double get_pot(const double) = 0;
  virtual double get_dpot(const double) = 0;
  virtual double get_dpot2(const double) = 0;
  virtual void   get_pot_dpot(const double, double&, double&) = 0;

  // Required members of mass model

  double get_mass(const double x, const double y, const double z) {
    return get_mass(z);}

  double get_density(const double x, const double y, const double z) {
    return get_density(z);}

  double get_pot(const double x, const double y, const double z) {
    return get_pot(z);}


  // Addiional member functions

  virtual double get_min_radius(void) = 0;
  virtual double get_max_radius(void) = 0;
  virtual double get_scale_height(void) = 0;
  virtual double distf(const double, const double V=0.0) = 0;
  virtual double dfde (const double, const double V=0.0) = 0;
  virtual double dfdv (const double, const double V=0.0) = 0;
};

class OneDModelTable : public OneDModel
{
protected:

  Spline1d mass, dens, pot;

  int even;
  int num;
  int numdf;
  double half_height;
  std::vector<double> params;

public:

  OneDModelTable() {};

  OneDModelTable(string filename, int PARM=0);

  OneDModelTable(int num, double *r, double *d, 
		 double *m, double *p, string ID = "" );

  // Required member functions

  double get_mass   (const double);
  double get_density(const double);
  double get_pot    (const double);
  double get_dpot   (const double);
  double get_dpot2  (const double);
  void   get_pot_dpot (const double, double&, double&);
  
  // Additional member functions

  const int get_num_param(void) { return params.size(); }
  const double get_param(int i) { return params[i-1]; }
  double get_scale_height(void) { return half_height; }

  double get_min_radius(void) { return mass.xlo(); }
  double get_max_radius(void) { return mass.xhi(); }
  int    grid_size(void) { return num; }

//  double distf(double E, double V);
//  double dfde(double E, double V);
//  double dfdv(double E, double V);
};

class LowIso : public OneDModelTable
{
private:
  double w0, Bfac, betak, gammak;
  double dispx, normx;

  void setup_model(void);

public:
  LowIso(string filename, double DISPX=0.159154943091895335768) :
    OneDModelTable(filename, 1) { 
      dispx = DISPX;
      setup_model();
    }

  double get_pot(const double);
  double get_dpot(const double);
  double get_dpot2(const double);
  void   get_pot_dpot(const double, double&, double&);

  double distf(const double E, const double V=0.0);
  double dfde(const double E, const double V=0.0);
  double dfdv(const double E, const double V=0.0);
};


class Sech2 : public OneDModel
{
private:
  double h;
  double dispz, dispx;
  double norm;

  static double HMAX;

public:

  Sech2(void)
  {
    dispz = 1.0;
    dispx = 1.0;
    //
    // Units: G = rho_o = 1
    //
    h = sqrt(dispz/(2.0*M_PI));
    norm = 1.0/( sqrt(2.0*M_PI*dispz) );
    dist_defined = true;
  }

  Sech2(const double DISPZ, const double DISPX=1.0)
  {
    dispz = DISPZ;
    dispx = DISPX;
    //
    // Units: G = rho_o = 1
    //
    h = sqrt(dispz/(2.0*M_PI));
    norm = 1.0/( sqrt(2.0*M_PI*dispz) );
    dist_defined = true;
  }
  
  double get_mass(const double z)
  {
    return 2.0*h/(1.0 + exp(-2.0*z/h));
  }

  double get_density(const double z)
  {
    double zz = fabs(z);
    double ret = 2.0*exp(-zz/h)/(1.0 + exp(-2.0*zz/h));
    return ret*ret;
  }

  double get_pot(const double z)
  {
    double zz = fabs(z);
    return 4.0*M_PI*h*(zz + h*log(1.0 + exp(-2.0*zz/h)) - h*M_LN2);
  }

  double get_dpot(const double z)
  {
    double zz = fabs(z);
    double ret = (1.0 - exp(-2.0*zz/h))/(1.0 + exp(-2.0*zz/h));
    return 4.0*M_PI*h* ret * z/(zz+1.0e-18);
  }

  double get_dpot2(const double z)
  {
    double zz = fabs(z);
    double ret = 2.0*exp(-zz/h)/(1.0 + exp(-2*zz/h));
    return 4.0*M_PI*ret*ret;
  }

  void get_pot_dpot(const double z, double& p, double& dp)
  {
    double zz = fabs(z);
    p = 4.0*M_PI*h*(zz + h*log(1.0 + exp(-2.0*zz/h)) - h*M_LN2);
    double ret = (1.0 - exp(-2.0*zz/h))/(1.0 + exp(-2.0*zz/h));
    dp = 4.0*M_PI*h* ret * z/(zz+1.0e-18);
  }

  double get_mass(const double x, const double y, const double z)
  {
    return get_mass(z);
  }

  double get_density(const double x, const double y, const double z)
  {
    return get_density(z);
  }

  double get_pot(const double x, const double y, const double z)
  {
    return get_pot(z);
  }

  double get_min_radius(void) { return 0.0; }
  double get_max_radius(void) { return HMAX*h; }
  double get_scale_height(void) { return h; }

  static void set_hmax(double hmax) { HMAX = hmax; }

  double distf(const double E, const double p)
  {
    return exp(-E/dispz - 0.5*p*p/dispx) * norm;
  }

  double dfde(const double E, const double p)
  {
    return -exp(-E/dispz - 0.5*p*p/dispx)/dispz * norm;
  }

  double dfdv(const double E, const double p)
  {
    return -exp(-E/dispz - 0.5*p*p/dispx)*p/dispx * norm;
  }

};



class Sech2mu : public OneDModel
{
private:
  double mu, h;
  double dispz, dispx;
  double dnorm, knorm, fnorm;

  static double HMAX;

public:

  Sech2mu(void)
  {
    mu    = 1.0;
    dispz = 1.0;
    dispx = 1.0;
    h     = 1.0;
    //
    // Units: G = mu = 1
    //
    dnorm = 0.25*mu/h;
    knorm = 2.0*M_PI*mu*h/dispz;
    fnorm = mu/(4.0*h*knorm*sqrt(2.0*M_PI*dispz));
    dist_defined = true;
  }

  Sech2mu(const double DISPZ, const double H, const double DISPX=1.0)
  {
    mu    = 1.0;
    dispz = DISPZ;
    dispx = DISPX;
    h     = H;
    //
    // Units: G = mu = 1
    //
    dnorm = 0.25*mu/h;
    knorm = 2.0*M_PI*mu*h/dispz;
    fnorm = mu/(4.0*h*knorm*sqrt(2.0*M_PI*dispz));
    dist_defined = true;
  }
  
  void setMu(double Sigma0)
  {
    mu = Sigma0;
    //
    // Units: G = 1, mu = Sigma0
    //
    dnorm = 0.25*mu/h;
    knorm = 2.0*M_PI*mu*h/dispz;
    fnorm = mu/(4.0*h*knorm*sqrt(2.0*M_PI*dispz));
    dist_defined = true;
  }

  double get_mass(const double z) {
    return mu*exp(z/h)/(1.0 + exp(z/h));
  }

  double get_density(const double z) {
    double zz  = fabs(z);
    double fac = exp(zz/h);
    return 4.0*dnorm/(fac + 1.0/fac + 2.0);
  }

  double get_pot(const double z) {
    double zz = fabs(z);
    return 2.0*dispz*(-M_LN2 + 0.5*zz/h + log(1.0 + exp(-zz/h)) ) - dispz*log(knorm);
  }

  double get_dpot(const double z) {
    double zz = fabs(z);
    double ret = (1.0 - exp(-zz/h))/(1.0 + exp(-zz/h));
    return 0.5*dispz*ret/h * z/(zz+1.0e-18);
  }

  double get_dpot2(const double z) {
    double zz = fabs(z);
    double ret = 2.0*exp(-0.5*zz/h)/(1.0 + exp(-zz/h));
    return 0.5*dispz*ret*ret/h/h;
  }

  void get_pot_dpot(const double z, double& p, double& dp) {
    double zz = fabs(z);
    p = 2.0*dispz*(-M_LN2 + 0.5*zz/h + log(1.0 + exp(-zz/h)) ) - dispz*log(knorm);
    double ret = (1.0 - exp(-zz/h))/(1.0 + exp(-zz/h));
    dp = 0.5*dispz*ret/h * z/(zz+1.0e-18);
  }

  double get_mass(const double x, const double y, const double z) {
    return get_mass(z);
  }

  double get_density(const double x, const double y, const double z) {
    return get_density(z);
  }

  double get_pot(const double x, const double y, const double z) {
    return get_pot(z);
  }

  double get_min_radius(void) { return 0.0; }
  double get_max_radius(void) { return HMAX*h; }
  double get_scale_height(void) { return h; }

  static void set_hmax(double hmax) { HMAX = hmax; }

  double distf(const double E, const double p) {
    return exp(-E/dispz - 0.5*p*p/dispx) * fnorm;
  }

  double dfde(const double E, const double p) {
    return -exp(-E/dispz - 0.5*p*p/dispx)/dispz * fnorm;
  }

  double dfdv(const double E, const double p) {
    return -exp(-E/dispz - 0.5*p*p/dispx)*p/dispx * fnorm;
  }

};


class Sech2Halo : public OneDModelTable
{
private:
  double h, rho0;
  double dispz, dispx;
  double dratio, hratio;
  double hh, rho0h;
  double norm, hmax;

  bool model_computed;

  static double HMAX;

  void reset();

public:

  static int NTABLE;
  static double OFFSET;
  static bool MU;

  Sech2Halo(void) {
    dispz = 1.0;
    dispx = 1.0;
    dratio = 0.0;
    hratio = 1.0;

    dist_defined = false;
    model_computed = false;
  }

  Sech2Halo(const double DISPZ, const double DRATIO, const double HRATIO,
	    const double DISPX=1.0);

  double get_pot(const double z) {
    if (!model_computed) reset();
    return OneDModelTable::get_pot(z) + 4.0*M_PI*hh*hh*rho0h*log(cosh(z/hh));
  }

  double get_dpot(const double z) {
    if (!model_computed) reset();
    return OneDModelTable::get_dpot(z) + 4.0*M_PI*hh*rho0h*tanh(z/hh);
  }

  double get_dpot2(const double z) {
    if (!model_computed) reset();
    double sech = 1.0/cosh(z/hh);
    return OneDModelTable::get_dpot2(z) + 4.0*M_PI*rho0h*sech*sech;
  }

  void get_pot_dpot(const double z, double& p, double& dp) {
    if (!model_computed) reset();
    OneDModelTable::get_pot_dpot(z, p, dp);
    p += 4.0*M_PI*hh*hh*rho0h*log(cosh(z/hh));
    dp += 4.0*M_PI*hh*rho0h*tanh(z/hh);
  }


  double get_min_radius(void) { return 0.0; }
  double get_max_radius(void) { return hmax*h; }
  double get_scale_height(void) { return h; }
  double get_scale_height_halo(void) { return hh; }
  double get_rho0(void) { return rho0; }
  double get_rho0_halo(void) { return rho0h; }

  static void set_hmax(double hmax) { HMAX = hmax; }
  
  double distf(const double E, const double p=0.0) {
    if (!model_computed) reset();
    return exp(-E/dispz - 0.5*p*p/dispx) * norm;
  }

  double dfde(const double E, const double p=0.0) {
    if (!model_computed) reset();
    return -exp(-E/dispz - 0.5*p*p/dispx)/dispz * norm;
  }

  double dfdv(const double E, const double p=0.0) {
    if (!model_computed) reset();
    return -exp(-E/dispz - 0.5*p*p/dispx)*p/dispx * norm;
  }

};

#endif

