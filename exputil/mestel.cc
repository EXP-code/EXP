#include <mestel.H>

// Implementation

MestelDisk::MestelDisk(double VROT, double RMIN, double RMAX)
{
  vrot = VROT;
  rot  = VROT*VROT;
  rmin = RMIN;
  rmax = RMAX;
  dim  = 2;
  
  setup_df(1.0);		// Default to sigma=1.0
  ModelID = "MestelDisk";
}
      

double MestelDisk::get_mass(const double r)
{ 
  if (r>0.0) return rot*r; 
  else bomb("radius cannot be zero!");
  return 0;
}

double MestelDisk::get_density(const double r)
{
  if (r>0.0) return rot/(2.0*M_PI*r);
  else bomb("radius cannot be zero!");
  return 0;
}

double MestelDisk::get_pot(const double r)
{ 
  if (r>0.0) return rot*log(r);
  else bomb("radius cannot be zero!");
  return 0;
}

double MestelDisk::get_dpot(const double r)
{
  if (r>0.0) return rot/r;
  else bomb("radius cannot be zero!");
  return 0;
}

double MestelDisk::get_dpot2(const double r)
{
  if (r>0.0) return -rot/(r*r);
  else bomb("radius cannot be zero!");
  return 0;
}

void MestelDisk::get_pot_dpot(const double r, double &ur, double &dur)
{
  if (r>0.0) {ur = rot*log(r); dur = rot/r;}
  else bomb("radius cannot be zero!");
}

void MestelDisk::setup_df(double sigma)
{ 
  sig2 = sigma*sigma;
  q = rot/sig2 - 1.0;
  F = rot/(4.0*M_PI) /
    ( sqrt(M_PI) *
      exp(lgamma(0.5*(q+1.0)) + (2.0 + q)*log(sigma) + 0.5*q*log(2.0)) );
  // 
  dist_defined = true;
}

double MestelDisk::distf(double E, double L)
{
  L = fabs(L);
  if (L==0.0) return 0.0;
  return F*pow(L, q) * exp(-E/sig2);
}

double MestelDisk::dfde(double E, double L)
{
  L = fabs(L);
  if (L==0.0) return 0.0;
  return -F*pow(L, q) * exp(-E/sig2)/sig2;
}

double MestelDisk::dfdl(double E, double L)
{
  int sgn = 1;
  if (L<0) {sgn=-1; L *= sgn;}
  if (L==0.0) return 0.0;
  return q*F*pow(L, q-1.0) * exp(-E/sig2) * sgn;
}
  
double MestelDisk::d2fde2(double E, double L)
{
  L = fabs(L);
  if (L<=0.0) return 0.0;
  return F*pow(L, q) * exp(-E/sig2)/sig2/sig2;
}

double TaperedMestelDisk::Tinner(double Jp)
{
  double fac = pow(Jp, nu);
  return fac/(Tifac + fac);
}

double TaperedMestelDisk::Touter(double Jp)
{
  return 1.0/(1.0 + pow(Jp/Tofac, mu));
}

double TaperedMestelDisk::dTinner(double Jp)
{
  double fac  = pow(Jp, nu);
  double fac2 = Tifac + fac;
  return nu*fac/Jp/(fac2*fac2);
}

double TaperedMestelDisk::dTouter(double Jp)
{
  double fac = pow(Jp/Tofac, mu);
  double fac2 = 1.0 + fac;
  return -mu*fac/Jp/(fac2*fac2);
}

TaperedMestelDisk::TaperedMestelDisk(double nu, double mu, double Ri, double Ro,
				     double VROT, double RMIN, double RMAX)
  : nu(nu), mu(mu), Ri(Ri), Ro(Ro), MestelDisk(VROT, RMIN, RMAX)
{
  Tifac = pow(Ri*vrot, nu);
  Tofac = Ro*vrot;
  ModelID = "TaperedMestelDisk";
}
      
double TaperedMestelDisk::get_density(const double r) {
  if (r>0.0) return rot/(2.0*M_PI*r) * Tinner(r) * Touter(r);
  else bomb("radius cannot be zero!");
  return 0;
}
  
double TaperedMestelDisk::get_mass(double R)
{
  if (not interp) {
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
      cum += 0.5*dr*(lst*(rr-dr) + cur*rr) * 2.0*M_PI;
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

double TaperedMestelDisk::distf(double E, double L)
{
  L = fabs(L);
  if (L==0.0) return 0.0;
  return F*pow(L, q)* Tinner(L)*Touter(L) * exp(-E/sig2);
}

double TaperedMestelDisk::dfde(double E, double L)
{
  L = fabs(L);
  if (L==0.0) return 0.0;
  return -F*pow(L, q)* Tinner(L)*Touter(L) * exp(-E/sig2)/sig2;
}

double TaperedMestelDisk::dfdl(double E, double L)
{
  int sgn=1;
  if (L<0) {sgn=-1; L *= sgn;}
  if (L==0.0) return 0.0;
  double Tfac = pow(L, q)*Tinner(L)*Touter(L);
  double dL = Tfac*(q/L + dTinner(L) + dTouter(L));
  return F* dL * exp(-E/sig2) * sgn;
}

double TaperedMestelDisk::d2fde2(double E, double L)
{
  L = fabs(L);
  if (L<=0.0) return 0.0;
  return F*pow(L, q)* Tinner(L)*Touter(L) * exp(-E/sig2)/sig2/sig2;
}
