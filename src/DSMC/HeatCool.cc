#include <fcntl.h>		// for open()

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

#include <climits>
#include <cfloat>
#include <cstdlib>
#include <cassert>

#include "global.H"
#include "HeatCool.H"

//
// Size of Romberg tableaux
//
#define MAXLEV 15


//
// xion tolerance 
//
#define TOLERANCE 1e-14
#define SMALLNUM  1.e-30

static double alphaj=-1.2;
static double fhydrogen=0.76;
static double jnu21=0.1;

bool   HeatCool::dbgCache = true;
bool   HeatCool::useCache = false;
double HeatCool::gp0_H;
double HeatCool::gp0_He;
double HeatCool::gp0_Hep;
double HeatCool::eps_H;
double HeatCool::eps_He;
double HeatCool::eps_Hep;

//
// rate coefficient of ionization of H by electron impact in cc/s 
// N.B. only accurate for 1e4 < T < 2e5 (Black 1981) 
// Modified as in Cen 1992 
//
double HeatCool::g_H(double t)
{
  if (t == 0.0) return 0.0;
  return 5.85e-11*sqrt(t)*exp(-157809.1 / t)/(1.0 + sqrt(t/1e5));
}

//
// rate coefficient of ionization of He by electron impact in cc/s 
// N.B. only accurate for 1e4 < T < 2e5 
// Modified as in Cen 1992 
//
double HeatCool::g_He(double t)
{
  if (t == 0.0) {
    return 0.0;
  }
  return 2.38e-11*sqrt(t)*exp(-285335.4 / t)/(1.0 + sqrt(t/1e5));
}

//
// rate coefficient of ionization of He+ by electron impact in cc/s 
// N.B. only accurate for 1e4 < T < 2e5 
// Modified as in Cen 1992 
//
double HeatCool::g_Hep(double t)
{
  if (t == 0.0) {
    return 0.0;
  }
  return 5.68e-12*sqrt(t)*exp(-631515.0 / t)/(1.0 + sqrt(t/1e5));
}

//
// radiative recombination rate coefficient for H (XXX units? most likely cc/s)
//
// Modified as in Cen 1992 
//
double HeatCool::a_Hp(double t)
{
  if (t > 0.0) {
    return (8.40e-11/sqrt(t))*pow(t/1e3, -0.2)/(1.0 + pow(t/1e6, 0.7));
  }
  return 0.0;
}

//
// radiative recombination rate coefficient for He
//
double HeatCool::a_Hep(double t)
{
  if (t > 0.0) {
    return 1.5e-10*pow(t, -0.6353);
  }
  return 0.0;
}

//
// radiative recombination rate coefficient for He+
// Modified as in Cen 1992 
//
double HeatCool::a_Hepp(double t)
{
  if (t > 0.0) {
    return (3.36e-10/sqrt(t))*pow(t/1e3, -0.2)/(1.0 + pow(t/1e6, 0.7));
  }
  return 0.0;
}

//
// dielectronic radiative recombination rate coefficient for He+
//
double HeatCool::a_p(double t)
{
  if (t == 0.0) {
    return 0.0;
  }
  return .0019*pow(t, -1.5)*exp(-4.7e5/t)*(exp(-9.4e4/t)*.3 + 1.0);
}

//
// cooling rate coefficient for recombination to H in erg*cc/s. 
// The cooling will scale as n_e*n_H+. 
// Modified as in Cen 1992 
//
double HeatCool::rate_Hp(double t)
{
  if (t == 0.0) {
    return 0.0;
  }
  return 8.70e-27*sqrt(t)*pow(t/1e3, -0.2)/(1.0 + pow(t/1e6, 0.7));
}

//
// cooling rate coefficient for recombination to He, collisional ionization
// of He+, collisional exitation of He+ and dielectronic recombination of 
// He+ in erg*cc/s. 
// The cooling will scale as n_e*n_He+. 
// Modified as in Cen 1992 
//
double HeatCool::rate_Hep(double t)
{
  if (t == 0.0) {
    return 0.0;
  }
  return 1.55e-26*pow(t, .3647)	// recombination 
    + 4.95e-22*sqrt(t)*exp(-631515.0/t)/(1.0 + sqrt(t/1e5)) // ionization 
    // exitation 
    + 5.54e-17*pow(t, -.397)*exp(-473638.0/t)/(1.0 + sqrt(t/1e5))
    // dielectric recombination 
    + 1.24e-13*pow(t, -1.5)*exp(-4.7e5/t)*(exp(-9.4e4/t)*.3 + 1.0);
}

//
// cooling rate coefficient for recombination to He+ in erg*cc/s. 
// The cooling will scale as n_e*n_He++. 
// Modified as in Cen 1992 
//
double HeatCool::rate_Hepp(double t)
{
  if (t == 0.0)
    return  0.0;
  return 3.48e-26*sqrt(t)*pow(t/1e3, -0.2)/(1.0 + pow(t/1e6, 0.7));
}

//
// cooling rate coefficient for collisional ionization and excitation of H in
// erg*cc/s. 
// The cooling will scale as n_e*n_H. 
// Modified as in Cen 1992 
//
double HeatCool::rate_H(double t)
{
  if (t == 0.0) return 0.0;
  // ionization 
  return 1.27e-21*sqrt(t)*exp(-157809.1/t)/(1.0 + sqrt(t/1e5))
    + 7.5e-19*exp(-118348.0/t)/(1.0 + sqrt(t/1e5));	// excitation 
}

//
// cooling rate coefficient for collisional excitation of H
// The cooling will scale as n_e*n_H. 
// Modified as in Cen 1992 
//
double HeatCool::rate_H_cex(double t)
{
  if (t == 0.0)
    return 0.0;
  return 7.5e-19*exp(-118348.0/t)/(1.0 + sqrt(t/1e5));	// excitation 
}

//
// cooling rate coefficient for collisional ionization of He in erg*cc/s. 
// The cooling will scale as n_e*n_He. 
// Modified as in Cen 1992 
//
double HeatCool::rate_He(double t)
{
  if (t == 0.0) {
    return 0.0;
  }
  return 9.38e-22*sqrt(t)*exp(-285335.4/t)/(1.0 + sqrt(t/1e5));
}

//
// cooling rate coefficient for bremsstrahlung for all ions in erg*cc/s. 
// The cooling will scale as n_e*(n_H+ + n_He+ + 4n_He++). 
//
double HeatCool::rate_br(double t)
{
  double gf = 1.1 + 0.34*exp(-pow(5.5-log10(t), 2.0)/3.0);
  return sqrt(t)*1.42e-27*gf;
}
//
// radiative emission coefficient for recombination to H in erg*cc/s. 
// (ie rate_Hp but includes potential energy as well) 
// The cooling will scale as n_e*n_H+. 
// Modified as in Cen 1992 
//
double HeatCool::radrate_Hp(double t)
{
  if (t == 0.0) {
    return 0.0;
  }
  return (8.70e-27 * t + 1.83e-21) / sqrt(t)
    * pow(t/1e3, -0.2)/(1.0 + pow(t/1e6, 0.7));
}


//
// radiative emission coefficient for recombination to He, 
// collisional ionization of He+, 
// collisional exitation of He+ and 
// dielectronic recombination of He+ in erg*cc/s. 
// The cooling will scale as n_e*n_He+. 
// Modified as in Cen 1992 
//
double HeatCool::radrate_Hep(double t)
{
  if (t == 0.0) {
    return 0.0;
  }
  return (1.55e-26 * t + 5.91e-21) * pow(t, -0.6353)	// recombination 
    + 4.95e-22*sqrt(t)*exp(-631515.0/t)/(1.0 + sqrt(t/1e5)) // ionization 
    // exitation 
    + 5.54e-17*pow(t, -.397)*exp(-473638.0/t)/(1.0 + sqrt(t/1e5))
    // dielectric recombination 
    + 1.99e-13*pow(t, -1.5)*exp(-4.7e5/t)*(exp(-9.4e4/t)*.3 + 1.0);
}

//
// radiative emission coefficient for recombination to He+ in erg*cc/s. 
// The cooling will scale as n_e*n_He++. 
// Modified as in Cen 1992 
//
double HeatCool::radrate_Hepp(double t)
{
  if (t == 0.0)
    return  0.0;
  return (3.48e-26 * t + 2.93e-20) / sqrt(t)
    * pow(t/1e3, -0.2)/(1.0 + pow(t/1e6, 0.7));
}

//
// calculate ionization equilibrium 
//
void HeatCool::xion(double t, double *x, double *x_1, double *x_2, double *x_3, double *p_ne)
{
  double a;
  int i;
  double zx;
  double zx1, zx2;
  double zx3;
  // recombination 
  double za_Hp, za_Hep, za_Hepp;
  // collisional ionization 
  double zg_H, zg_He, zg_Hep;	
  double n_e;
  double old_n_e;
  
  zx = .05;
  zx1 = 1.0;
  zx2 = 0.0;
  n_e = (1.0 - zx)*n_h;
  old_n_e = n_e;
  
  za_Hp = a_Hp(t);
  zg_H = g_H(t);
  za_Hep = a_Hep(t) + a_p(t);
  zg_He = g_He(t);
  za_Hepp = a_Hepp(t);
  zg_Hep = g_Hep(t);
  for (i = 1; i <= 50; ++i) {
    // solve for neutral H fraction given constant n_e 
    if (g0 + (zg_H + za_Hp)*n_e != 0.0) {
      zx = za_Hp*n_e/(g0 + (zg_H + za_Hp)*n_e);
    } else {
      zx = 1.0;
    }
    if (g1 + zg_He*n_e > 1e-50) {
      a = za_Hep*n_e/(g1 + zg_He*n_e);
    } else {
      zx1 = 1.0;
      zx2 = 0.0;
      zx3 = 0.0;
      n_e = (1.0 - zx)*n_h;
      if(old_n_e == 0.0 || fabs(n_e - old_n_e)/n_h < TOLERANCE)
	break;
      n_e = .5*(n_e + old_n_e);
      old_n_e = n_e;
      continue;
    }
    if (g2 + (zg_Hep + (a + 1.0)*za_Hepp)*n_e != 0.0) 
      {
	zx2 = za_Hepp*n_e/(g2 + (zg_Hep + (a + 1.0)*za_Hepp)*n_e);
	zx1 = zx2 * a;
	zx3 = zx2*(g2 + zg_Hep*n_e)/za_Hepp/n_e;
      } else {
	zx3 = 0.0;
	zx2 = 0.0;
	zx1 = 1.0;
      }
    n_e = (1.0 - zx + zx2*r + zx3*2.0*r)*n_h;
    if(fabs(n_e - old_n_e)/n_h < TOLERANCE)
      break;
    n_e = .5*(n_e + old_n_e);
    old_n_e = n_e;
  }
  if (i > 50)
    {
      ostringstream sout;
      sout << "<too many iterations in xion>" << endl;
      sout << "t: " << t << ", n_h: " << n_h << ", n_e: " << n_e << endl;
      throw (sout.str());
    }

  zx = min(1.,zx);
  *x = zx;
  *x_1 = zx1;
  *x_2 = zx2;
  *x_3 = zx3;
  *p_ne = n_e;
}

double HeatCool::fac1(double t)
{
  double tinv;
  double eps;
  double fac;
  
  tinv = 1./t;
  eps = sqrt(tinv-1.);
  fac = exp(4.-4.*atan(eps)/eps)/(1.-exp(-2.*M_PI/eps))*pow(t,(alphaj+3.));
  return fac;
}

double HeatCool::fac2(double t)
{
  double tinv;
  double eps;
  double fac;
  
  tinv = 1./t;
  eps = sqrt(tinv-1.);
  fac = exp(4.-4.*atan(eps)/eps)/(1.-exp(-2.*M_PI/eps))*pow(t,(alphaj+3.));
  fac = fac*(tinv-1.);
  return fac;
}


double HeatCool::romberg_o(std::function<double(double)> func,
			   double a, double b, double eps)
{
  double tllnew;
  double tll;
  double tlk[MAXLEV+1];
  int n = 1;
  int nsamples = 1;
  
  tlk[0] = tllnew = (b-a)*func((b-a)/2.0);
  tll = 0.5*HUGE_VAL;
  
  while((fabs(tllnew-tll) > fabs(tllnew)*eps) && (n < MAXLEV)) {
    // midpoint rule
    
    double deltax;
    double tlktmp;
    int i;
    
    nsamples *= 3;
    deltax = (b-a)/nsamples;
    tlktmp = tlk[0];
    tlk[0] = tlk[0]/3.0;
    
    for(i = 0; i < nsamples/3; i++) {
      tlk[0] += deltax*func(a + (3*i + 0.5)*deltax);
      tlk[0] += deltax*func(a + (3*i + 2.5)*deltax);
    }
    
    // Romberg extrapolation.
    
    for(i = 0; i < n; i++) {
      double tlknew = (pow(9.0, i+1)*tlk[i] - tlktmp)
	/(pow(9.0, i+1) - 1.0);
      
      tlktmp = tlk[i+1];
      tlk[i+1] = tlknew;
    }
    tll = tllnew;
    tllnew = tlk[n];
    n++;
  }
  
  assert(n < MAXLEV);
  return tllnew;
}

double HeatCool::spline_int(double* xa, double* ya, double* y2a, int n, double x)
  // x table 
  // y table 
  // y'' table 
  // number of entries 
  // evaluation point 
{
  int lo,hi,k;
  double h,b,a;
  double result;
  
  lo = 0;
  hi = n-1;
  while (hi-lo > 1) {
    k=(hi + lo) >> 1;	// bisect interval 
    if (xa[k] > x) hi = k;
    else lo = k;
  }
  h = xa[hi] - xa[lo];
  assert(h > 0.0);
  
  a = (xa[hi] - x)/h;
  b = (x - xa[lo])/h;
  result = a*ya[lo] + b*ya[hi] // linear interpolation 
    + ((a*a*a - a)*y2a[lo] + (b*b*b - b)*y2a[hi])*(h*h)/6.0;
  return result;
}

void
HeatCool::ionize()
{
  double a0;
  double planck;
  double ev;
  double e0_H;
  double e0_He;
  double e0_Hep;
  double gint;
  double eint;
  double at;
  double beta;
  double s;
  
  //  ionize -- compute photo-ionization rates and photo-ionization heating
  //
  //     input:
  //       * alphaj = index of UV background, J \propto (nu_L/nu)**alphaj
  //     output:
  //       * gp0_H, gp0_He, gp0_Hep = photo-ionization rates for H, He, Hep 
  
  a0     = 6.30e-18;
  planck = 6.63e-27;
  ev     = 1.60e-12;
  e0_H   = 13.60*ev;
  e0_He  = 24.59*ev;
  e0_Hep = 54.42*ev;
  
  //
  //  evaluate dimensionless integrals needed for H and He+ rates (see notes) 
  //
  gint = romberg_o(&fac1,0.,1.,1e-6);
  gp0_H = a0*gint/planck;
  gp0_Hep = gp0_H*pow((e0_H/e0_Hep),alphaj)/4.;
  eint = romberg_o(&fac2,0.,1.,1e-6);
  eps_H = a0*eint*(e0_H/planck);
  eps_Hep = eps_H*pow((e0_H/e0_Hep),(alphaj-1.))/4.;
  
  at = 7.83e-18;
  beta = 1.66;
  s = 2.05;
  
  gp0_He = (at/planck)*pow((e0_H/e0_He),alphaj)*(beta/(alphaj+s)+(1.-beta)/
						 (alphaj+s+1));
  eps_He = (e0_He/planck)*at*pow((e0_H/e0_He),alphaj)*(beta/(alphaj+s-1)+
						       (1-2*beta)/(alphaj+s)-(1-beta)/(alphaj+s+1));
  return;
}

void HeatCool::initialize()
{
  ionize();

  gp0_H *= jnu21 * 4.0 * M_PI * 1.e-21 ;
  gp0_He *= jnu21 * 4.0 * M_PI * 1.e-21 ;
  gp0_Hep *= jnu21 * 4.0 * M_PI * 1.e-21 ;
  eps_H *= jnu21 * 4.0 * M_PI * 1.e-21 ;
  eps_He *= jnu21 * 4.0 * M_PI * 1.e-21 ;
  eps_Hep *= jnu21 * 4.0 * M_PI * 1.e-21 ;
}


HeatCool::HeatCool(double nmin, double nmax, double tmin, double tmax,
		   unsigned nnum, unsigned tnum, string cache)
  : Nmin(nmin), Nmax(nmax), Tmin(tmin), Tmax(tmax), 
    Nnum(nnum), Tnum(tnum), Cache(cache), table(true)
{

  if (useCache && readCache()) return;

  dN = (log(Nmax) - log(Nmin))/(Nnum-1);
  dT = (log(Tmax) - log(Tmin))/(Tnum-1);

  double N, T;
  Rates r;

  for (unsigned n=0; n<Nnum; n++) {

    N = Nmin*exp(dN*n);

    vector<Rates> tmp;
    for (unsigned t=0; t<Tnum; t++) {

      T = Tmin*exp(dT*t);

      compute(N, T);
      
      if (crate > 0.0) 
	r.crate = log(crate);
      else
	r.crate = log(std::numeric_limits<double>::min());

      if (hrate > 0.0)
	r.hrate = log(hrate);
      else
	r.hrate = log(std::numeric_limits<double>::min());

      if (cmpcrate > 0.0)
	r.cmpcrate = cmpcrate;
      else
	r.cmpcrate = log(std::numeric_limits<double>::min());

      r.trate = trate;

      tmp.push_back(r);
    }
    data.push_back(tmp);
  }

  if (dbgCache) {
				// Try to make a lock file
    int fd = open(".HeatCoolDbg.lck", O_WRONLY | O_CREAT | O_EXCL, S_IRWXU);

    if (fd>=0) {	 // Owner writes the cache and clears the lock
      debugCache();
      remove(".HeatCoolDbg.lck");
    }

  }

  if (useCache) {
				// Try to make a lock file
    int fd = open(".HeatCoolExp.lck", O_WRONLY | O_CREAT | O_EXCL, S_IRWXU);

    if (fd<0) return;		// Someone else has it . . . 

				// Owner writes the cache and clears the lock
    writeCache();
    remove(".HeatCoolExp.lck");
  }

}


//
// Bilinear interpolation
//
void HeatCool::setPoint(double N, double T)
{
  if (!table) {
    compute(N, T);
    return;
  }

  unsigned n_indx, t_indx;
  double aN, bN, aT, bT;
  double lN=log(N/Nmin), lT=log(T/Tmin);

  if (N < Nmin) n_indx = 0;
  else if (N >= Nmax) n_indx = Nnum-2;
  else n_indx = min<unsigned>(Nnum-2, static_cast<unsigned>(floor(lN/dN)));

  aN = (dN*(n_indx+1) - lN)/dN;
  bN = (lN - dN*n_indx    )/dN;

  if (T <  Tmin) t_indx = 0;
  else if (T >= Tmax) t_indx = Tnum-2;
  else t_indx = min<unsigned>(Tnum-2, static_cast<unsigned>(floor(lT/dT)));

  aT = (dT*(t_indx+1) - lT)/dT;
  bT = (lT - dT*t_indx    )/dT;
  
  crate = 
    aN*(aT*data[n_indx  ][t_indx  ].crate + 
	bT*data[n_indx  ][t_indx+1].crate   ) +
    bN*(aT*data[n_indx+1][t_indx  ].crate + 
	bT*data[n_indx+1][t_indx+1].crate   ) ;

  crate = exp(crate);

  hrate = 
    aN*(aT*data[n_indx  ][t_indx  ].hrate + 
	bT*data[n_indx  ][t_indx+1].hrate   ) +
    bN*(aT*data[n_indx+1][t_indx  ].hrate + 
	bT*data[n_indx+1][t_indx+1].hrate   ) ;

  hrate = exp(hrate);

  cmpcrate = 
    aN*(aT*data[n_indx  ][t_indx  ].cmpcrate + 
	bT*data[n_indx  ][t_indx+1].cmpcrate   ) +
    bN*(aT*data[n_indx+1][t_indx  ].cmpcrate +
	bT*data[n_indx+1][t_indx+1].cmpcrate   ) ;

  cmpcrate = exp(cmpcrate);

  trate = 
    aN*(aT*data[n_indx  ][t_indx  ].trate + 
	bT*data[n_indx  ][t_indx+1].trate   ) +
    bN*(aT*data[n_indx+1][t_indx  ].trate +
	bT*data[n_indx+1][t_indx+1].trate   ) ;
}


double HeatCool::CoolRate(double N, double T)
{
  if (T<Tmin) return 0.0;

  if (!table) {
    compute(N, T);
    return crate;
  }

  unsigned n_indx, t_indx;
  double aN, bN, aT, bT;
  double lN=log(N/Nmin), lT=log(T/Tmin);

  if (N < Nmin) n_indx = 0;
  else if (N >= Nmax) n_indx = Nnum-2;
  else n_indx = min<unsigned>(Nnum-2, static_cast<unsigned>(floor(lN/dN)));

  aN = (dN*(n_indx+1) - lN)/dN;
  bN = (lN - dN*n_indx    )/dN;

  if (T <  Tmin) t_indx = 0;
  else if (T >= Tmax) t_indx = Tnum-2;
  else t_indx = min<unsigned>(Tnum-2, static_cast<unsigned>(floor(lT/dT)));

  aT = (dT*(t_indx+1) - lT)/dT;
  bT = (lT - dT*t_indx    )/dT;
  
  double ans = 
    aN*(aT*data[n_indx  ][t_indx  ].crate + 
	bT*data[n_indx  ][t_indx+1].crate   ) +
    bN*(aT*data[n_indx+1][t_indx  ].crate + 
	bT*data[n_indx+1][t_indx+1].crate   ) ;

  if (ans>-0.5) {
    cout << "Craziness in HeatCool" << endl;
  }

  return exp(ans);
}


void HeatCool::compute(double density, double temp)
{
  double h0, h1, h2;
  double x, x_1, x_2, x_3, f_e, n_e;
  double y;
  
  n_h = density*fhydrogen;
  
  //------------------------------------------------------- 
  //  crate is cooling rate in units of (erg cm^3/s) 
  //  hrate is heating rate 
  // 
  //  alphaj is spectral index J \propto \nu^-\alpha 
  // 
  //  t  is temperature in K 
  // 
  //  j_nu is photoionizing flux at Lyman limit 
  //  in units of 10^{-21} (per steradian) 
  // 
  //  n_H is density of hydrogen ions (Not helium) 
  // 
  //  y is helium abundance by mass 
  //  r is ratio of helium ions to hydrogen ions (all ionizations) 
  //  g0 is photoionization rate coefficient? 
  //  g1 is photoionization rate coefficient? 
  //  g2 is photoionization rate coefficient? 
  // 
  //  h0 is photoionization heating rate coefficient? 
  //  h1 is photoionization heating rate coefficient? 
  //  h2 is photoionization heating rate coefficient? 
  //  x is fraction of H that is NEUTRAL 
  //  x_1 is fraction of He that is NEUTRAL 
  //  x_2 is fraction of He that is singly ionized 
  
  y = 1.0 - fhydrogen;
  r = y / 4.0 / (1.0 - y);

  g0 = gp0_H;
  g1 = gp0_He;
  g2 = gp0_Hep;

  g0 = 0.0;
  g1 = 0.0;
  g2 = 0.0;
  
  h0 = eps_H;
  h1 = eps_He;
  h2 = eps_Hep;
  
  try {
    xion(temp, &x, &x_1, &x_2, &x_3, &n_e);
  } catch (string msg) {
    crate = hrate = cmpcrate = trate = 0.0;
    cerr << msg << ", no heating or cooling for this step" << endl;
    return;
  }

  f_e = 1.0 - x + x_2*r + (x_3)*2.0*r;
  crate = f_e*rate_Hp(temp)*(1.0 - x);
  crate += f_e*x_2*r*rate_Hep(temp);
  crate += f_e*(x_3)*r*rate_Hepp(temp);
  crate += f_e*x*rate_H(temp);
  crate += f_e*x_1*r*rate_He(temp);
  crate += f_e*(1.0 - x + x_2*r + (x_3)*4.0*r)*rate_br(temp);
  crate *= n_h*n_h;
  hrate = h0*x + h1*x_1*r + h2*x_2*r;
  hrate *= n_h;
  
  //
  // see Ikeuchi and Ostriker 1986, Ap. J. 301, 522. 
  //
  cmpcrate = 5.41e-36*f_e*temp*n_h;
  trate = crate + cmpcrate - hrate;
}


bool HeatCool::readCache()
{
  if (Cache.size()==0) return false;

  ifstream in((outdir + Cache).c_str());
  if (!in) {
    cerr << "HeatCool::readCache: could not open <" << outdir + Cache << ">"
	 << endl;
    return false;
  }

  // Read header info

  char cid[14];
  in.read((char*)cid, 14*sizeof(char));
  if (string(cid).compare("HeatCool data")) {
    cerr << "HeatCool::readCache: bad magic string <" << cid << ">"
	 << endl;
    return false;
  }
  
  double nmin, nmax, tmin, tmax;
  unsigned nnum, tnum;
  Rates r;

  in.read((char *)&nmin, sizeof(double));
  in.read((char *)&nmax, sizeof(double));
  in.read((char *)&tmin, sizeof(double));
  in.read((char *)&tmax, sizeof(double));
  in.read((char *)&nnum, sizeof(unsigned));
  in.read((char *)&tnum, sizeof(unsigned));

  if (fabs(Nmin-nmin)>1.0e-18 ||
      fabs(Nmax-nmax)>1.0e-18 ||
      fabs(Tmin-tmin)>1.0e-18 ||
      fabs(Tmax-tmax)>1.0e-18 ||
      Nnum != nnum            ||
      Tnum != tnum            ) {
    cerr << "HeatCool::readCache: header value mismatch"
	 << endl;
    return false;
  }

  in.read((char *)&dN, sizeof(double));
  in.read((char *)&dT, sizeof(double));

  for (unsigned n=0; n<Nnum; n++) {
    vector<Rates> tmp;
    for (unsigned t=0; t<Tnum; t++) {
      in.read((char *)&r.crate,    sizeof(double));
      in.read((char *)&r.hrate,    sizeof(double));
      in.read((char *)&r.cmpcrate, sizeof(double));
      in.read((char *)&r.trate,    sizeof(double));
      if (!in) {
	cerr << "HeatCool::readCache: error reading data values"
	     << endl;
	return false;
      }
    }
    data.push_back(tmp);
  }
  
  return true;
}


void HeatCool::writeCache()
{
  if (Cache.size()==0) return;

  ofstream out(string(outdir + Cache).c_str());
  if (!out) {
    cerr << "HeatCool::writeCache: could not open <" << outdir + Cache << ">"
	 << endl;
    return;
  }

  char cid[14] = "HeatCool data";

  out.write((const char*)cid,   14*sizeof(char));
  out.write((const char*)&Nmin, sizeof(double));
  out.write((const char*)&Nmax, sizeof(double));
  out.write((const char*)&Tmin, sizeof(double));
  out.write((const char*)&Tmax, sizeof(double));
  out.write((const char*)&Nnum, sizeof(unsigned));
  out.write((const char*)&Tnum, sizeof(unsigned));
  out.write((const char*)&dN,   sizeof(double));
  out.write((const char*)&dT,   sizeof(double));
  
  for (unsigned n=0; n<Nnum; n++) {
    for (unsigned t=0; t<Tnum; t++) {
      out.write((const char*)&data[n][t].crate,    sizeof(double));
      out.write((const char*)&data[n][t].hrate,    sizeof(double));
      out.write((const char*)&data[n][t].cmpcrate, sizeof(double));
      out.write((const char*)&data[n][t].trate,    sizeof(double));
    }
  }
  // Done!
}

void HeatCool::debugCache()
{
  ostringstream dfile;
  dfile << outdir << "HeatCool_" << runtag << ".debug";
  ofstream out(dfile.str().c_str());
  if (!out) {
    cerr << "HeatCool::writeCache: could not open <" << dfile.str() << ">"
	 << endl;
    return;
  }

  for (unsigned n=0; n<Nnum; n++) {

    double N = Nmin + exp(dN*n);

    for (unsigned t=0; t<Tnum; t++) {

      double T = Tmin * exp(dT*t);

      out << setw(18) << N
	  << setw(18) << T
	  << setw(18) << data[n][t].crate
	  << setw(18) << data[n][t].hrate
	  << setw(18) << data[n][t].cmpcrate
	  << setw(18) << data[n][t].trate
	  << endl;

    }

    out << endl;
  }
  // Done!
}
