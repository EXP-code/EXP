// This may look like C code, but it is really -*- C++ -*-

const char rcsid[] = "$Id$";

/*

*/

#define FRECS 16		// # of rectangles for frequency integrals
#define tol 1.0e-8		// tolerance for root finder
#define ZFRAC 0.001		// fraction below minimum grid for tp search
#define RMAXF 3.0		// Fraction beyond grid for r_apo search
#define tolnr 1.0e-10		// r_apo location using Newton-Raphson refinement
#define TOLEPI 1.0e-3		// Cut over for epicylic theory


#include <math.h>
#include <string>
#include <numerical.h>
#include <Vector.h>
#include <interp.h>
#include "massmodel.h"
#include "orbit.h"

using namespace std;

				// Global variables
static AxiSymModel *mm;
static double EE, JJ, KK;
static double ap, am, sp, sm;


void SphericalOrbit::compute_freq(void)
{
  double accum0, accum1, accum2, r, s, t, dt, ur, dudr, cost, tmp;
#ifdef DEBUG
  double tmp2;
#endif
  int i;
  double Ecirc(double), denom(double);

  mm = model;
  EE = energy;
  KK = kappa;

  /*  Find turning points  */

  double xmin, xmax;
  xmin = ZFRAC*model->get_min_radius();
  xmax = model->get_max_radius();
  if ( Ecirc(xmin)*Ecirc(xmax) < 0.0)
    r_circ = zbrent(Ecirc, xmin, xmax, tol);
  else				// Circular radius outside mass distribution
    r_circ = -0.5*model->get_mass(model->get_max_radius())/EE;

  dudr = model->get_dpot(r_circ);
  jmax = sqrt(r_circ*r_circ*r_circ*dudr);
  JJ = jmax*kappa;

//  r_apo = zbrent(denom,r_circ,model->get_max_radius(),tol);

  xmin = r_circ;
  xmax = model->get_max_radius();
  if ( denom(xmin)*denom(xmax) < 0.0)
    r_apo = zbrent(denom, xmin, xmax, tol);
  if (r_circ < model->get_max_radius()) {
    r_apo = zbrent(denom, xmin, RMAXF*xmax, tol);
  }
  else {	 		// Circular radius outside mass distribution
    double r0 = -0.5*model->get_mass(model->get_max_radius())/EE;
    r_apo = r0 + sqrt(r0*r0 - 0.5*JJ*JJ/EE);
				// Newton-Raphson refinement (slow . . .)
    double f, df;
    for (int i=0; i<200; i++) {
      f = 2.0*(EE - model->get_pot(r_apo)) - JJ*JJ/(r_apo*r_apo);
      if (fabs(f)<tolnr) break;
      df = -2.0*model->get_dpot(r_apo) + JJ*JJ/(r_apo*r_apo*r_apo);
      r_apo += -f/df;
    }

#ifdef DEBUG
    if (fabs(f)>tolnr) {
      cerr << "compute_freq: r_apo did not converge [f=" << f << "]\n";
    }
#endif
  }

  if (denom(ZFRAC*model->get_min_radius())*denom(r_circ) >= 0.0)
    r_peri = 0.99*ZFRAC*model->get_min_radius();
  else
    r_peri = zbrent(denom,ZFRAC*model->get_min_radius(),r_circ,tol);

  /*  Ok, prepare to do freq. integrals */

  ap = 0.5*(r_apo + r_peri);
  am = 0.5*(r_apo - r_peri);
  sp = ap/(r_apo*r_peri);
  sm = am/(r_apo*r_peri);

  /* Do using centered rectangle technique */
  accum0 = 0.0;
  accum1 = 0.0;
  accum2 = 0.0;
  dt = M_PI/FRECS;
  for (i=0, t=0.5*(dt-M_PI); i<FRECS; i++, t+=dt) {
    r = ap + am*sin(t);
    ur = model->get_pot(r);
    cost = cos(t);
#ifdef DEBUG
    tmp2 = 2.0*(energy-ur) - (jmax*jmax*kappa*kappa)/(r*r);
    if (tmp2 < 0.0) {
      cerr << "compute_freq: denominator [1] out of bounds\n";
      cerr << "compute_freq: E=" << EE << "  K=" << KK << 
	" r=" << r << " i=" << i << "/" << FRECS << " val=" << tmp2 << endl;

      tmp2 = 2.0*(energy-model->get_pot(r_peri)) - 
	(jmax*jmax*kappa*kappa)/(r_peri*r_peri);
      cerr << "compute_freq: denom(r_peri)=" << tmp2 << endl;

      tmp2 = 2.0*(energy-model->get_pot(r_apo)) - 
	(jmax*jmax*kappa*kappa)/(r_apo*r_apo);
      cerr << "compute_freq: denom(r_apo)=" << tmp2 << endl;

      tmp2 = (-model->get_dpot(r_peri) +
	(jmax*jmax*kappa*kappa)/(r_peri*r_peri*r_peri))*(r-r_peri);
      cerr << "compute_freq: denom_exp(r_peri)=" << tmp2 << endl;

      tmp2 = (-model->get_dpot(r_apo) +
	(jmax*jmax*kappa*kappa)/(r_apo*r_apo*r_apo))*(r-r_apo);
      cerr << "compute_freq: denom_exp(r_apo)=" << tmp2 << endl;

    }
#endif
    tmp = sqrt(2.0*(energy-ur) - (jmax*jmax*kappa*kappa)/(r*r));
    accum0 += cost * tmp;
    accum1 += cost / tmp;
    s = sp + sm*sin(t);
    ur = model->get_pot(1.0/s);
#ifdef DEBUG
    tmp2 = 2.0*(energy-ur) - (jmax*jmax*kappa*kappa*s*s);
    if (tmp2 < 0.0) {
      cerr << "compute_freq: denominator [2] out of bounds\n";
      cerr << "compute_freq: E=" << EE << "  K=" << KK << 
	" r=" << r << " i=" << i << "/" << FRECS << " val=" << tmp2 << endl;

      tmp2 = 2.0*(energy-model->get_pot(r_peri)) - 
	(jmax*jmax*kappa*kappa)/(r_peri*r_peri);
      cerr << "compute_freq: denom(r_peri)=" << tmp2 << endl;

      tmp2 = 2.0*(energy-model->get_pot(r_apo)) - 
	(jmax*jmax*kappa*kappa)/(r_apo*r_apo);
      cerr << "compute_freq: denom(r_apo)=" << tmp2 << endl;

      tmp2 = (-model->get_dpot(r_peri) +
	(jmax*jmax*kappa*kappa)/(r_peri*r_peri*r_peri))*(r-r_peri);
      cerr << "compute_freq: denom_exp(r_peri)=" << tmp2 << endl;

      tmp2 = (-model->get_dpot(r_apo) +
	(jmax*jmax*kappa*kappa)/(r_apo*r_apo*r_apo))*(r-r_apo);
      cerr << "compute_freq: denom_exp(r_apo)=" << tmp2 << endl;

    }
#endif
    accum2 += cost/sqrt(2.0*(energy-ur) - (jmax*jmax*kappa*kappa*s*s));
  }
  
  freq.setsize(1,3);
  action.setsize(1,3);

  freq[1] = M_PI/(am*accum1*dt);
  freq[2] = freq[1]*jmax*kappa * sm*accum2*dt/M_PI;
  freq[3] = 0.0;
  freq_defined = true;

  action[1] = am*accum0*dt/M_PI;
  action[2] = jmax*kappa;
  action[3] = 0.0;
  action_defined = true;
}

#undef FRECS
#undef tol

/* Function to iteratively locate radius of circular orbit with energy EE */
double Ecirc(double r)
{
  double ur, dudr, dif;
		
  mm->get_pot_dpot(r,ur,dudr);
  dif =  EE - 0.5*r*dudr - ur;
  return(dif);
}

/* Function to iteratively locate turning points for orbit (EE,JJ) */

double denom(double r)
{
  double ur = mm->get_pot(r);
  return 2.0*(EE-ur)*r*r - JJ*JJ;
}

/* for testing---frequencies in the epicyclic approx */

void SphericalOrbit::compute_freq_epi(void)
{
  double d2udr2;

  if (!freq_defined) compute_freq();
  
  d2udr2 = model->get_dpot2(r_circ);
  freq[1] = sqrt(3.0*jmax*jmax*kappa*kappa/(r_circ*r_circ*r_circ*r_circ) + 
		 d2udr2);
  freq[2] = jmax*kappa/(r_circ*r_circ);

}

double  rombe2(double a, double b, double (*f) (double), int n);

double dtp, dtm;

void SphericalOrbit::compute_angles(void)
{
  double accum1,accum2,r, s, t, sl, tl;
  double fw1(double t),ff(double t);
  int i;

  l1s = l2s = 0;

  angle_grid.t.setsize(1,2,0,recs-1);
  angle_grid.w1.setsize(1,2,0,recs-1);
  angle_grid.dw1dt.setsize(1,2,0,recs-1);
  angle_grid.f.setsize(1,2,0,recs-1);
  angle_grid.r.setsize(1,2,0,recs-1);
  
  if (!Gkn) {
    Gkn = new LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm =  M_PI;
  } 
  
  if (freq_defined) {
    mm = model;
    EE = energy;
    KK = kappa;
    JJ = jmax*kappa;
    ap = 0.5*(r_apo + r_peri);
    am = 0.5*(r_apo - r_peri);
    sp = ap/(r_apo*r_peri);
    sm = am/(r_apo*r_peri);
  }
  else
    compute_freq();

  accum1 = 0.0;
  accum2 = 0.0;
  tl = -0.5*M_PI;
  sl =  0.5*M_PI;
  for (i=0; i<recs; i++) {
    t = dtp + dtm*Gkn->knot(i+1);
    r = ap + am*sin(t);
    accum1 += rombe2(tl,t,fw1,nbsct);
    s = asin((1.0/r - sp)/sm);
    accum2 += rombe2(sl,s,ff, nbsct);

    angle_grid.t[1][i] = t;
    angle_grid.w1[1][i] = freq[1]*accum1;
    angle_grid.dw1dt[1][i] = freq[1]*fw1(t);
    angle_grid.f[1][i]  = freq[2]*accum1 + JJ*accum2;
    angle_grid.r[1][i]  = r;

    tl = t;
    sl = s;
  }

  Spline(angle_grid.w1[1], angle_grid.t[1], 1.0e30, 1.0e30,
	 angle_grid.t[2]);

  Spline(angle_grid.w1[1], angle_grid.dw1dt[1], 1.0e30, 1.0e30,
	 angle_grid.dw1dt[2]);

  Spline(angle_grid.w1[1], angle_grid.f[1], 1.0e30, 1.0e30,
	 angle_grid.f[2]);

  Spline(angle_grid.w1[1], angle_grid.r[1], 1.0e30, 1.0e30,
	 angle_grid.r[2]);

  angle_grid.num = recs;

  angle_defined = true;
}

#define tol 1.0e-8       /* tolerance for radicand */

double fw1(double t)
{
  double r, ur, dudr, radcand, sgn;

  r = ap + am*sin(t);
  ur = mm->get_pot(r);
  radcand = 2.0*(EE-ur)-JJ*JJ/(r*r);

  if (radcand < tol) {
  /* values near turning points */
    if (t<0)
      sgn = -1.0;
    else
      sgn = 1.0;
    r = ap+sgn*am;
    dudr = mm->get_dpot(r);
    return sqrt( am*sgn/(dudr-JJ*JJ/(r*r*r)) );
  }
  return am*cos(t)/sqrt(radcand);
}

double ff(double t)
{
  double s, ur, dudr, radcand, sgn;

  s = sp + sm*sin(t);
  ur = mm->get_pot(1.0/s);
  radcand = 2.0*(EE-ur) - JJ*JJ*s*s;

  if (radcand < tol) {
  /* values near turning points */
    if (t<0)
      sgn = -1.0;
    else
      sgn = 1.0;
    s = sp+sgn*sm;
    dudr = mm->get_dpot(1.0/s);
    return sqrt( -sm*s*s*sgn/(dudr-JJ*JJ*s*s*s) );
  }
  return sm*cos(t)/sqrt(radcand);
}

/* epicyclic test */

void SphericalOrbit::compute_angles_epi(void)
{
  double r,t,a,a2,fac,ur;
  double Jcirc(double r);
  int i;

  angle_grid.t.setsize(1,2,0,recs-1);
  angle_grid.w1.setsize(1,2,0,recs-1);
  angle_grid.dw1dt.setsize(1,2,0,recs-1);
  angle_grid.f.setsize(1,2,0,recs-1);
  angle_grid.r.setsize(1,2,0,recs-1);
  
  if (!Gkn) {
    Gkn = new LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm =  M_PI;
  } 
  else if (Gkn->get_n() != recs) {
    delete Gkn;
    Gkn = new LegeQuad(recs);
  } 
  
  if (freq_defined) {
    mm = model;
    EE = energy;
    KK = kappa;
  }
  else
    compute_freq_epi();

  JJ = jmax*kappa;
  ur = mm->get_pot(r_circ);
  a2 = 2.0*(EE-0.5*JJ*JJ/(r_circ*r_circ)-ur)/(freq[1]*freq[1]);
  a = sqrt(a2);

  for (i=0; i<recs; i++) {
    t = dtp + dtm*Gkn->knot(i+1);
    r = r_circ - a*cos(t);
    fac = sqrt(fabs(a2 - (r-r_circ)*(r-r_circ)));
    angle_grid.t[1][i] = t;
    angle_grid.w1[1][i] = t;
    angle_grid.dw1dt[1][i] = 1.0;
    angle_grid.f[1][i]  = -2.0*freq[2]/(freq[1]*r_circ) * fac;
    angle_grid.r[1][i]  = r;

  }

  Spline(angle_grid.w1[1], angle_grid.t[1], 1.0e30, 1.0e30,
	 angle_grid.t[2]);

  Spline(angle_grid.w1[1], angle_grid.dw1dt[1], 1.0e30, 1.0e30,
	 angle_grid.dw1dt[2]);

  Spline(angle_grid.w1[1], angle_grid.f[1], 1.0e30, 1.0e30,
	 angle_grid.f[2]);

  Spline(angle_grid.w1[1], angle_grid.r[1], 1.0e30, 1.0e30,
	 angle_grid.r[2]);

  angle_grid.num = recs;

  angle_defined = true;

}


double Jcirc(double r)
{
  double dudr = mm->get_dpot(r);
  return (JJ*JJ - r*r*r*dudr);
}
#undef tol



void SphericalOrbit::compute_biorth(void)
{
  
  if (!angle_defined) compute_angles();

  angle_grid.fr.setsize(0, recs-1, 1, nmax);

  for (int j=0; j<recs; j++) {
    if (RECUR) {
      biorth->potl(nmax, l, biorth->r_to_rb(angle_grid.r[1][j]), 
		   angle_grid.fr[j]);
    } else
      for (int i=1; i<=nmax; i++) {
	angle_grid.fr[j][i] = 
	  biorth->potl(i, l, biorth->r_to_rb(angle_grid.r[1][j]));
      }
  }

  biorth_defined = true;
}


double SphericalOrbit::pot_trans(int l1, int l2, double (*func)(double))
{

  if (!Gkn) {
    Gkn = new LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm =  M_PI;
  }

  if (Gkn->get_n() != recs) {
    delete Gkn;
    Gkn = new LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm =  M_PI;
  } 
  
  if (!angle_defined) compute_angles();

  double accum = 0.0;

  if (kappa < 1.0-TOLEPI) {
    for (int i=0; i<angle_grid.num; i++)
      accum += Gkn->weight(i+1)*angle_grid.dw1dt[1][i]*
	cos(angle_grid.w1[1][i]*l1 + angle_grid.f[1][i]*l2)*
	  func(angle_grid.r[1][i]);

    accum *= dtm/M_PI;
  }
  else {
    if (l1 == 0)
      accum = func(angle_grid.r[1][(angle_grid.num-1)/2]);
    else
      accum = 0.0;
  }

  return accum;
}

double SphericalOrbit::pot_trans(int l1, int l2, int n)
{

  if (!Gkn) {
    Gkn = new LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm =  M_PI;
  }
  else if (Gkn->get_n() != recs) {
    delete Gkn;
    Gkn = new LegeQuad(recs);
  }

  if (!angle_defined) compute_angles();
  if (!biorth_defined) compute_biorth();

  if (l1s==0 && l2s==0) cosvec.setsize(0, angle_grid.num-1);
  if (l1 != l1s || l2 != l2s) {
    for (int i=0; i<angle_grid.num; i++)
      cosvec[i] = cos(angle_grid.w1[1][i]*l1 + angle_grid.f[1][i]*l2);
    l1s = l1;
    l2s = l2;
  }

  double accum = 0.0;

  if (kappa < 1.0-TOLEPI) {
    for (int i=0; i<angle_grid.num; i++)
      accum += Gkn->weight(i+1)*angle_grid.dw1dt[1][i] * cosvec[i] * 
	angle_grid.fr[i][n];

    accum *= dtm/M_PI;
  }
  else {
    if (l1 == 0)
      accum = angle_grid.fr[(angle_grid.num-1)/2][n];
    else
      accum = 0.0;
  }

  return accum;
}

void SphericalOrbit::pot_trans(int l1, int l2, Vector& t)
{

  if (!Gkn) {
    Gkn = new LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm =  M_PI;
  }
  else if (Gkn->get_n() != recs) {
    delete Gkn;
    Gkn = new LegeQuad(recs);
  }

  if (!angle_defined) compute_angles();
  if (!biorth_defined) compute_biorth();

  if (l1s==0 && l2s==0) cosvec.setsize(0, angle_grid.num-1);
  if (l1 != l1s || l2 != l2s) {
    for (int i=0; i<angle_grid.num; i++)
      cosvec[i] = cos(angle_grid.w1[1][i]*l1 + angle_grid.f[1][i]*l2);
    l1s = l1;
    l2s = l2;
  }

  t.zero();
  int nm = t.gethigh();
  double tmpi;

  if (kappa < 1.0-TOLEPI) {
    for (int i=0; i<angle_grid.num; i++) {
      tmpi = Gkn->weight(i+1)*angle_grid.dw1dt[1][i] * cosvec[i];
      for (int n=1; n<=nm; n++)
	t[n] += tmpi * angle_grid.fr[i][n];
    }
    t *= dtm/M_PI;
  }
  else {
    if (l1 == 0)
      for (int n=1; n<=nm; n++)
	t[n] = angle_grid.fr[(angle_grid.num-1)/2][n];
  }

}



