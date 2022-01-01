#define FRECS 16		// # of rectangles for frequency integrals
#define TOLEPI 1.0e-3		// Cut over for epicylic theory

#include <cstdlib>
#include <sstream>
#include <cmath>

#include <numerical.H>
#include <interp.H>
#include <massmodel.H>
#include <orbit.H>
#include <biorth.H>

				// Global variables
static AxiSymModPtr mm;
static double EE, JJ, KK;
static double ap, am, sp, sm;

double SphericalOrbit::tol = 1.0e-8;  // root finder tolerance
double SphericalOrbit::ZFRAC = 0.001; // fraction below minimum grid for tp search
double SphericalOrbit::RMAXF = 3.0;
double SphericalOrbit::tolnr = 1.0e-10;	// r_apo location using Newton-Raphson refinement

void SphericalOrbit::compute_freq(void)
{
  const string rname("compute_freq");
  double accum0, accum1, accum2, r, s, t, dt, ur, dudr, cost, tmp;
#ifdef DEBUG
  double tmp2;
#endif
  int i;
  double Ecirc(double), denom(double);

  mm = model;
  EE = energy;
  KK = kappa;

  //  Find turning points

  double xmin, xmax;
  xmin = ZFRAC*model->get_min_radius();
  xmax = model->get_max_radius();

  if ( Ecirc(xmin)*Ecirc(xmax) < 0.0) {
    try {
      r_circ = zbrent(Ecirc, xmin, xmax, tol);
    }
    catch (const char *error) {
      cerr << "SphericalOrbit::compute_freq @ r_circ: model=" 
	   << model->ModelID         << endl
	   << " energy="   << energy << endl
	   << " kappa="    << kappa  << endl
	   << " xmin="     << xmin   << endl
	   << " xmax="     << xmax   << endl;
      throw error;
    }
  }
  else {		  // Circular radius outside mass distribution
    r_circ = -0.5*model->get_mass(model->get_max_radius())/EE;
    if (r_circ < model->get_min_radius()) {
      r_circ = model->get_min_radius();
    }
    std::cerr << "SphericalOrbit::compute_freq warning:"
	      << " Circular radius outside mass distribution" << std::endl;
  }

  dudr = model->get_dpot(r_circ);
  if (dudr>0.0) jmax = sqrt(r_circ*r_circ*r_circ*dudr);
  else          jmax = 0.0;

  JJ = jmax*kappa;

  xmin = r_circ;
  xmax = model->get_max_radius();
#ifdef DEBUG
  bool used_asymp = false;
#endif
  if ( denom(xmin)*denom(xmax) < 0.0) {
    try {
      r_apo = zbrent(denom, xmin, xmax, tol);
    }
    catch (const char *error) {
      cerr << "SphericalOrbit::compute_freq @ r_apo[1]: model=" 
	   << model->ModelID        << endl
	   << " energy=" << energy  << endl
	   << " kappa="  << kappa   << endl
	   << " xmin="  << xmin     << endl
	   << " xmax="  << xmax     << endl;
      throw error;
    }
  }
  if (r_circ < model->get_max_radius()) {
    try {
      r_apo = zbrent(denom, xmin, RMAXF*xmax, tol);
    }
    catch (const char *error) {
      cerr << "SphericalOrbit::compute_freq @ r_apo[2]: model=" 
	   << model->ModelID        << endl
	   << " energy=" << energy  << endl
	   << " pot(max)=" << model->get_pot(xmax) << endl
	   << " kappa="  << kappa   << endl
	   << " xmin="   << xmin    << endl
	   << " xmax="   << xmax*RMAXF << endl
	   << " rc check=" << energy - model->get_pot(r_circ) - jmax*jmax/(2.0*r_circ*r_circ) << endl;
      throw error;
    }
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
      ostringstream msg;
      msg << "r_apo did not converge [f=" << f << "]";
      warn(rname, msg.str());
    }
    used_asymp = true;
#endif
  }

  if (denom(ZFRAC*model->get_min_radius())*denom(r_circ) >= 0.0) {
    r_peri = 0.99*ZFRAC*model->get_min_radius();
  }
  else {
    try {
      r_peri = zbrent(denom, ZFRAC*model->get_min_radius(), r_circ, tol);
    }
    catch (const char *error) {
      cerr << "SphericalOrbit::compute_freq @ r_peri: model=" << model->ModelID
	   << " energy=" << energy
	   << " kappa="  << kappa
	   << " xmin="  << ZFRAC*model->get_min_radius()
	   << " xmax="  << r_circ << endl;
      throw error;
    }
  }

  //  Ok, prepare to do freq. integrals 

  ap = 0.5*(r_apo + r_peri);
  am = 0.5*(r_apo - r_peri);
  sp = ap/(r_apo*r_peri);
  sm = am/(r_apo*r_peri);

  // Do using centered rectangle technique 
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
      ostringstream msg;

      msg << "\t\tdenominator [1] out of bounds" << endl;
      msg << "E=" << EE << "  K=" << KK << 
	" r=" << r << " i=" << i << "/" << FRECS << " val=" << tmp2 << endl;

      tmp2 = 2.0*(energy-model->get_pot(r_peri)) - 
	(jmax*jmax*kappa*kappa)/(r_peri*r_peri);
      msg << "\t\tdenom(r_peri)=" << tmp2 << endl;

      tmp2 = 2.0*(energy-model->get_pot(r_apo)) - 
	(jmax*jmax*kappa*kappa)/(r_apo*r_apo);
      msg << "\t\tdenom(r_apo)=" << tmp2 << endl;

      tmp2 = (-model->get_dpot(r_peri) +
	(jmax*jmax*kappa*kappa)/(r_peri*r_peri*r_peri))*(r-r_peri);
      msg << "\t\tdenom_exp(r_peri)=" << tmp2 << endl;

      tmp2 = (-model->get_dpot(r_apo) +
	(jmax*jmax*kappa*kappa)/(r_apo*r_apo*r_apo))*(r-r_apo);
      msg << "\t\tdenom_exp(r_apo)=" << tmp2 << endl;

      warn(rname, msg.str());
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
      ostringstream msg;

      msg << "\t\t denominator [2] out of bounds" << endl;
      msg << "\t\tE=" << EE << "  K=" << KK << 
	" r=" << r << " i=" << i << "/" << FRECS << " val=" << tmp2 << endl;

      tmp2 = 2.0*(energy-model->get_pot(r_peri)) - 
	(jmax*jmax*kappa*kappa)/(r_peri*r_peri);
      msg << "\t\tdenom(r_peri)=" << tmp2 << endl;

      tmp2 = 2.0*(energy-model->get_pot(r_apo)) - 
	(jmax*jmax*kappa*kappa)/(r_apo*r_apo);
      msg << "\t\tdenom(r_apo)=" << tmp2 << endl;

      tmp2 = (-model->get_dpot(r_peri) +
	(jmax*jmax*kappa*kappa)/(r_peri*r_peri*r_peri))*(r-r_peri);
      msg << "\t\tdenom_exp(r_peri)=" << tmp2 << endl;

      tmp2 = (-model->get_dpot(r_apo) +
	(jmax*jmax*kappa*kappa)/(r_apo*r_apo*r_apo))*(r-r_apo);
      msg << "\t\tcompute_freq: denom_exp(r_apo)=" << tmp2 << endl;

      warn(rname, msg.str());
    }
#endif
    accum2 += cost/sqrt(2.0*(energy-ur) - (jmax*jmax*kappa*kappa*s*s));
  }
  
  freq[0] = M_PI/(am*accum1*dt);
  freq[1] = freq[0]*jmax*kappa * sm*accum2*dt/M_PI;
  freq[2] = 0.0;
  freq_defined = true;

  action[0] = am*accum0*dt/M_PI;
  action[1] = jmax*kappa;
  action[2] = 0.0;
  action_defined = true;
}

#undef FRECS
#undef tol

// Function to iteratively locate radius of circular orbit with energy EE 
double Ecirc(double r)
{
  double ur, dudr, dif;
		
  mm->get_pot_dpot(r, ur, dudr);
  dif =  EE - 0.5*r*dudr - ur;
  return(dif);
}

// Function to iteratively locate turning points for orbit (EE,JJ) 

double denom(double r)
{
  double ur = mm->get_pot(r);
  return 2.0*(EE-ur)*r*r - JJ*JJ;
}

// for testing---frequencies in the epicyclic approx 

void SphericalOrbit::compute_freq_epi(void)
{
  double d2udr2;

  if (!freq_defined) compute_freq();
  
  d2udr2 = model->get_dpot2(r_circ);
  freq[0] = sqrt(3.0*jmax*jmax*kappa*kappa/(r_circ*r_circ*r_circ*r_circ) + 
		 d2udr2);
  freq[1] = jmax*kappa/(r_circ*r_circ);

}

extern 
double  rombe2(double a, double b, double (*f) (double), int n);

double dtp, dtm;

void SphericalOrbit::compute_angles(void)
{
  double accum1,accum2,r, s, t, sl, tl;
  double fw1(double t), ff(double t);
  int i;

  l1s = l2s = 0;

  angle_grid.t.    resize(2, recs);
  angle_grid.w1.   resize(2, recs);
  angle_grid.dw1dt.resize(2, recs);
  angle_grid.f.    resize(2, recs);
  angle_grid.r.    resize(2, recs);
  
  if (Gkn.n == 0) {
    Gkn = LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm = M_PI;
  } 
  else if (Gkn.n != recs) {
    Gkn = LegeQuad(recs);
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
  for (int i=0; i<recs; i++) {
    t = dtp + dtm*Gkn.knot(i);
    r = ap + am*sin(t);
    accum1 += rombe2(tl,t,fw1,nbsct);
    s = asin((1.0/r - sp)/sm);
    accum2 += rombe2(sl,s,ff, nbsct);

    angle_grid.t(0, i) = t;
    angle_grid.w1(0, i) = freq[0]*accum1;
    angle_grid.dw1dt(0, i) = freq[0]*fw1(t);
    angle_grid.f(0, i)  = freq[1]*accum1 + JJ*accum2;
    angle_grid.r(0, i)  = r;

    tl = t;
    sl = s;
  }


  Eigen::VectorXd work(angle_grid.t.row(0).size());
  
  Spline(angle_grid.w1.row(0), angle_grid.t.row(0), 1.0e30, 1.0e30, work);
  angle_grid.t.row(1) = work;

  Spline(angle_grid.w1.row(0), angle_grid.dw1dt.row(0), 1.0e30, 1.0e30, work);
  angle_grid.dw1dt.row(1) = work;

  Spline(angle_grid.w1.row(0), angle_grid.f.row(0), 1.0e30, 1.0e30, work);
  angle_grid.f.row(1) = work;

  Spline(angle_grid.w1.row(0), angle_grid.r.row(0), 1.0e30, 1.0e30, work);
  angle_grid.r.row(1) = work;

  angle_grid.num = recs;

  angle_defined = true;
}

#define tol 1.0e-8       // tolerance for radicand 

double fw1(double t)
{
  double r, ur, dudr, radcand, sgn;

  r = ap + am*sin(t);
  ur = mm->get_pot(r);
  radcand = 2.0*(EE-ur)-JJ*JJ/(r*r);

  if (radcand < tol) {
  // values near turning points 
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
  // values near turning points 
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

//
// epicyclic test 
//
void SphericalOrbit::compute_angles_epi(void)
{
  double r,t,a,a2,fac,ur;
  double Jcirc(double r);
  int i;

  angle_grid.t.    resize(2, recs);
  angle_grid.w1.   resize(2, recs);
  angle_grid.dw1dt.resize(2, recs);
  angle_grid.f.    resize(2, recs);
  angle_grid.r.    resize(2, recs);
  
  if (Gkn.n == 0) {
    Gkn = LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm = M_PI;
  }
  else if (Gkn.n != recs) {
    Gkn = LegeQuad(recs);
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
  a2 = 2.0*(EE-0.5*JJ*JJ/(r_circ*r_circ)-ur)/(freq[0]*freq[0]);
  a = sqrt(a2);

  for (int i=0; i<recs; i++) {
    t = dtp + dtm*Gkn.knot(i);
    r = r_circ - a*cos(t);
    fac = sqrt(fabs(a2 - (r-r_circ)*(r-r_circ)));
    angle_grid.t(0, i) = t;
    angle_grid.w1(0, i) = t;
    angle_grid.dw1dt(0, i) = 1.0;
    angle_grid.f(0, i)  = -2.0*freq[1]/(freq[0]*r_circ) * fac;
    angle_grid.r(0, i)  = r;

  }

  Eigen::VectorXd work(angle_grid.t.row(0).size());

  Spline(angle_grid.w1.row(0), angle_grid.t.row(0), 1.0e30, 1.0e30, work);
  angle_grid.t.row(1) = work;

  Spline(angle_grid.w1.row(0), angle_grid.dw1dt.row(0), 1.0e30, 1.0e30, work);
  angle_grid.dw1dt.row(1) = work;

  Spline(angle_grid.w1.row(0), angle_grid.f.row(0), 1.0e30, 1.0e30, work);
  angle_grid.f.row(1) = work;

  Spline(angle_grid.w1.row(0), angle_grid.r.row(0), 1.0e30, 1.0e30, work);
  angle_grid.r.row(1) = work;

  angle_grid.num = recs;
  angle_defined  = true;
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

  angle_grid.fr.resize(recs, nmax);

  Eigen::VectorXd work(nmax);

  for (int j=0; j<recs; j++) {
    if (RECUR) {
      biorth->potl(nmax, l, biorth->r_to_rb(angle_grid.r(0, j)), work);
      angle_grid.fr.row(j) = work;
    } else
      for (int i=0; i<nmax; i++) {
	angle_grid.fr(j, i) = 
	  biorth->potl(i, l, biorth->r_to_rb(angle_grid.r(0, j)));
      }
  }

  biorth_defined = true;
}


double SphericalOrbit::pot_trans(int l1, int l2, double (*func)(double))
{

  if (Gkn.n == 0) {
    Gkn = LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm = M_PI;
  }
  else if (Gkn.n != recs) {
    Gkn = LegeQuad(recs);
  }
  
  if (!angle_defined) compute_angles();

  double accum = 0.0;

  if (kappa < 1.0-TOLEPI) {
    for (int i=0; i<angle_grid.num; i++)
      accum += Gkn.weight(i)*angle_grid.dw1dt(0, i)*
	cos(angle_grid.w1(0, i)*l1 + angle_grid.f(0, i)*l2)*
	  func(angle_grid.r(0, i));

    accum *= dtm/M_PI;
  }
  else {
    if (l1 == 0)
      accum = func(angle_grid.r(0, (angle_grid.num-1)/2));
    else
      accum = 0.0;
  }

  return accum;
}

double SphericalOrbit::pot_trans(int l1, int l2, int n)
{

  if (Gkn.n == 0) {
    Gkn = LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm = M_PI;
  }
  else if (Gkn.n != recs) {
    Gkn = LegeQuad(recs);
  }

  if (!angle_defined)  compute_angles();
  if (!biorth_defined) compute_biorth();

  if (l1s==0 && l2s==0) cosvec.resize(angle_grid.num);
  if (l1 != l1s || l2 != l2s) {
    for (int i=0; i<angle_grid.num; i++)
      cosvec[i] = cos(angle_grid.w1(0, i)*l1 + angle_grid.f(0, i)*l2);
    l1s = l1;
    l2s = l2;
  }

  double accum = 0.0;

  if (kappa < 1.0-TOLEPI) {
    for (int i=0; i<angle_grid.num; i++)
      accum += Gkn.weight(i)*angle_grid.dw1dt(0, i) * cosvec[i] * 
	angle_grid.fr(i, n);

    accum *= dtm/M_PI;
  }
  else {
    if (l1 == 0)
      accum = angle_grid.fr((angle_grid.num-1)/2, n);
    else
      accum = 0.0;
  }

  return accum;
}

void SphericalOrbit::pot_trans(int l1, int l2, Eigen::VectorXd& t)
{

  if (Gkn.n == 0) {
    Gkn = LegeQuad(recs);
    dtp = -0.5*M_PI;
    dtm = M_PI;
  }
  else if (Gkn.n != recs) {
    Gkn = LegeQuad(recs);
  }

  if (!angle_defined) compute_angles();
  if (!biorth_defined) compute_biorth();

  if (l1s==0 && l2s==0) cosvec.resize(angle_grid.num);
  if (l1 != l1s || l2 != l2s) {
    for (int i=0; i<angle_grid.num; i++)
      cosvec[i] = cos(angle_grid.w1(0, i)*l1 + angle_grid.f(0, i)*l2);
    l1s = l1;
    l2s = l2;
  }

  t.setZero();
  int nm = t.size();
  double tmpi;

  if (kappa < 1.0-TOLEPI) {
    for (int i=0; i<angle_grid.num; i++) {
      tmpi = Gkn.weight(i)*angle_grid.dw1dt(0, i) * cosvec[i];
      for (int n=0; n<nm; n++)
	t[n] += tmpi * angle_grid.fr(i, n);
    }
    t *= dtm/M_PI;
  }
  else {
    if (l1 == 0)
      for (int n=0; n<nm; n++)
	t[n] = angle_grid.fr((angle_grid.num-1)/2, n);
  }
}

