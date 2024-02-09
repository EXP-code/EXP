#ifndef _mestel_H
#define _mestel_H

#include <massmodel.H>

//! Infiinte Mestel Disk
class MestelDisk : public AxiSymModel
{
protected:
  
  double vrot;
  double rot;
  double rmin;
  double rmax;
  double sig2;
  double q;
  double F;

public:

  //! Constructor
  MestelDisk(double VROT = 1.0, 
	     double RMIN = 1.0e-6, double RMAX = 1.0e6)
  {
    vrot = VROT;
    rot  = VROT*VROT;
    rmin = RMIN;
    rmax = RMAX;
    dim  = 2;

    dist_defined = true;
    ModelID = "MestelDisk";
  }
      

  //@{
  //! Required member functions
  virtual double get_mass(const double r) { 
    if (r>0.0) return rot*r; 
    else bomb("radius cannot be zero!");
    return 0;
  }

  virtual double get_density(const double r) {
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
  //@}
  
  //@{
  //! Addiional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }

  virtual void setup_df(double sigma) { 
    sig2 = sigma*sigma;
    q = rot/sig2 - 1.0;
    F = rot/(4.0*M_PI) / ( sqrt(M_PI) * 
       exp(lgamma(0.5*(q+1.0)) + (2.0 + q)*log(sigma) + 0.5*q*log(2.0)) );
  }

  virtual double distf(double E, double L) {
    L = fabs(L);
    if (L==0.0) return 0.0;
    return F*pow(L, q) * exp(-E/sig2);
  }

  virtual double dfde(double E, double L) {
    L = fabs(L);
    if (L==0.0) return 0.0;
    return -F*pow(L, q) * exp(-E/sig2)/sig2;
  }

  virtual double dfdl(double E, double L) {
    int sgn=1;
    if (L<0) {sgn=-1; L *= sgn;}
    if (L==0.0) return 0.0;
    return q*F*pow(L, q-1.0) * exp(-E/sig2) * sgn;
  }
  
  virtual double d2fde2(double E, double L) {
    L = fabs(L);
    if (L<=0.0) return 0.0;
    return F*pow(L, q) * exp(-E/sig2)/sig2/sig2;
  }
  //@}

};


//! Doubly tapered Mestel Disk
class TaperedMestelDisk : public MestelDisk
{
private:

  std::shared_ptr<Linear1d> interp;

protected:
  //! Taper parameters
  double nu, mu;

  //! Inner radius
  double Ri;

  //! Taper factor
  double Tifac;

  //! Inner taper function
  double Tinner(double Jp)
  {
    double fac = pow(Jp, nu);
    return fac/(Tifac + fac);
  }

  //! Outer taper function
  double Touter(double Jp)
  {
    return 1.0/(1.0 + pow(Jp/vrot, mu));
  }

  //! Deriv of inner taper function
  double dTinner(double Jp)
  {
    double fac  = pow(Jp, nu);
    double fac2 = Tifac + fac;
    return Tifac*nu/Jp/fac2;
  }

  //! Deriv of outer taper function
  double dTouter(double Jp)
  {
    double fac = pow(Jp/vrot, mu);
    double fac2 = 1.0 + fac;
    return -nu*fac/Jp/fac2;
  }


public:

  //! Constructor
  TaperedMestelDisk(double nu, double mu, double Ri,
		    double VROT = 1.0, 
		    double RMIN = 1.0e-6, double RMAX = 1.0e6) :
    nu(nu), mu(mu), Ri(Ri), MestelDisk(VROT, RMIN, RMAX)
    {
      Tifac = pow(Ri*vrot, nu);
      ModelID = "TaperedMestelDisk";
    }
      
  //! Overloaded density function
  double get_density(const double r) {
    if (r>0.0) return rot/(2.0*M_PI*r) * Tinner(r) * Touter(r);
    else bomb("radius cannot be zero!");
    return 0;
  }
  
  double get_mass(double R)
  {
    if (interp) {
      // Compute a grid spacing
      double dr = 0.01*Ri;
      int N = std::floor((rmax - rmin)/dr) + 1;
      dr = (rmax - rmin)/N;

      // Storage for table
      Eigen::VectorXd vecR(N+1), vecM(N+1);

      double lst = get_density(rmin);
      double cum = 0.0;

      // Make the table
      //
      vecR(0) = rmin;
      vecM(0) = 0.0;

      for (int i=1; i<=N; i++) {
	double rr  = rmin + dr*i;
	double cur = get_density(rr);
	cum += 0.5*dr*(lst + cur);
	lst = cur;

	vecR(i) = rr;
	vecM(i) = cum;
      }

      // Make the interpolation function
      //
      interp = std::make_shared<Linear1d>(vecR, vecM);
    }

    return interp->eval(R);
  }

  //@{
  //! Overloaded DF member functions

  double distf(double E, double L)
  {
    L = fabs(L);
    if (L==0.0) return 0.0;
    return F*pow(L, q)* Tinner(L)*Touter(L) * exp(-E/sig2);
  }

  double dfde(double E, double L)
  {
    L = fabs(L);
    if (L==0.0) return 0.0;
    return -F*pow(L, q)* Tinner(L)*Touter(L) * exp(-E/sig2)/sig2;
  }

  double dfdl(double E, double L)
  {
    int sgn=1;
    if (L<0) {sgn=-1; L *= sgn;}
    if (L==0.0) return 0.0;
    double Tfac = pow(L, q)*Tinner(L)*Touter(L);
    double dL = Tfac*(q/L + dTinner(L) + dTouter(L));
    return F* dL * exp(-E/sig2) * sgn;
  }
  
  double d2fde2(double E, double L)
  {
    L = fabs(L);
    if (L<=0.0) return 0.0;
    return F*pow(L, q)* Tinner(L)*Touter(L) * exp(-E/sig2)/sig2/sig2;
  }
  //@}

};



#endif
