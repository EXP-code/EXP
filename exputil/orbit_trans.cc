#define FRECS 16		// # of rectangles for frequency integrals
#define TOLEPI 1.0e-3		// Cut over for epicylic theory

#include <functional>
#include <cstdlib>
#include <sstream>
#include <cmath>

#include "numerical.H"
#include "interp.H"
#include "massmodel.H"
#include "orbit.H"
#include "biorth.H"

int    SphericalOrbit::Nseg = 40;     // Number of trial segments for root finding 
double SphericalOrbit::tol = 1.0e-8;  // root finder tolerance
double SphericalOrbit::ZFRAC = 0.001; // fraction below minimum grid for tp search
double SphericalOrbit::RMAXF = 3.0;
double SphericalOrbit::tolnr = 1.0e-10;	// r_apo location using Newton-Raphson refinement

std::tuple<double, double, bool> SphericalOrbit::search
(std::function<double(double)> func, double rmin, double rmax)
{
  bool use_log = false;
  if (rmin > 0.0) {
    use_log = true;
    rmin = log(rmin);
    rmax = log(rmax);
  }
  
  double dx = (rmax - rmin)/Nseg;
  bool cond = false;
  
  // Loop variables
  double xend = rmin, xbeg, ylst;
  double ycur = func(use_log ? exp(xend) : xend);

  // Now, increment the end point and test for sign change in the
  // interval
  for (int n=1; n<=Nseg; n++) {
    xbeg = xend;
    ylst = ycur;
    xend += dx;
    ycur = func(use_log ? exp(xend) : xend);

    // Found an interval, return it
    if (ycur*ylst <= 0.0) {
      cond = true;
      break;
    }
  }
  
  if (use_log) {		// Scale back to linear
    xbeg = exp(xbeg);
    xend = exp(xend);
  }
				// Done!
  return std::tuple<double, double, bool>(xbeg, xend, cond);
}

void SphericalOrbit::compute_freq(void)
{
  const string rname("compute_freq");
  double accum0, accum1, accum2, r, s, t, dt, ur, dudr, cost, tmp;
#ifdef DEBUG
  double tmp2;
#endif
  int i;

  //  Find turning points
  //
  double xmin, xmax;
  xmin = ZFRAC*model->get_min_radius();
  xmax = model->get_max_radius();

  // Functor whose zero locates radius of circular orbit
  //
  auto Ecirc = [&](double r)
  {
    double ur, dudr;
    model->get_pot_dpot(r, ur, dudr);
    return  energy - 0.5*r*dudr - ur;
  };

  // New strategy: break into Nseg segments to find good starting
  // points
  //
  {
    auto [xbeg, xend, cond] = search(Ecirc, xmin, xmax);

    if ( Ecirc(xbeg)*Ecirc(xend) < 0.0) {
      try {
	r_circ = zbrent(Ecirc, xbeg, xend, tol);
      }
      catch (const std::exception& error) {
	std::ostringstream mesg;
	mesg << "SphericalOrbit::compute_freq: " << error.what() << std::endl
	     << "SphericalOrbit::compute_freq @ r_circ: model=" 
	     << model->ModelID           << std::endl
	     << "** energy="   << energy << std::endl
	     << "** kappa="    << kappa  << std::endl
	     << "** xmin="     << xbeg   << std::endl
	     << "** xmax="     << xend   << std::endl
	     << "** f(xmin)="  << Ecirc(xbeg) << std::endl
	     << "** f(xmax)="  << Ecirc(xend) << std::endl
	     << "** cond="     << std::boolalpha << cond << std::endl;
	throw std::runtime_error(mesg.str());
      }
    }
    else {		  // Circular radius outside mass distribution
      r_circ = -0.5*model->get_mass(model->get_max_radius())/energy;
      if (r_circ < model->get_min_radius()) {
	r_circ = model->get_min_radius();
      }
      /*
	std::string message("SphericalOrbit::compute_freq warning: Circular radius outside mass distribution");
	throw std::runtime_error(message);
      */
    }
  }

  dudr = model->get_dpot(r_circ);
  if (dudr>0.0) jmax = sqrt(r_circ*r_circ*r_circ*dudr);
  else          jmax = 0.0;

  // Functor whose zero locates turning points for orbit
  //
  auto denom = [&](double r)
  {
    double ur = model->get_pot(r);
    double JJ = jmax*kappa;
    return 2.0*(energy-ur)*r*r - JJ*JJ;
  };

  xmin = r_circ;
  xmax = model->get_max_radius();
#ifdef DEBUG
  bool used_asymp = false;
#endif

  {
    auto [xbeg, xend, cond] = search(denom, xmin, xmax);

    if ( denom(xbeg)*denom(xend) < 0.0) {
      try {
	r_apo = zbrent(denom, xbeg, xend, tol);
      }
      catch (const std::exception& error) {
	std::ostringstream mesg;
	mesg << "SphericalOrbit::compute_freq: " << error.what() << std::endl
	     << "SphericalOrbit::compute_freq @ r_apo[1]: model=" 
	     << model->ModelID          << std::endl
	     << "** energy=" << energy  << std::endl
	     << "** kappa="  << kappa   << std::endl
	     << "** energy=" << energy  << std::endl
	     << "** kappa="  << kappa   << std::endl
	     << "** J="      << kappa*jmax << std::endl
	     << "** rmin="   << xbeg    << std::endl
	     << "** rmax="   << xend    << std::endl
	     << "** cond="   << std::boolalpha << cond << std::endl;

	throw std::runtime_error(mesg.str());
      }
    }
  }
  
  if (r_circ < model->get_max_radius()) {

    auto [xbeg, xend, cond] = search(denom, r_circ, RMAXF*xmax);

    if (denom(xbeg)*denom(xend) < 0.0) {
      try {
	r_apo = zbrent(denom, xbeg, xend, tol);
      }
      catch (const std::exception& error) {
	std::ostringstream mesg;
	mesg << "SphericalOrbit::compute_freq: " << error.what() << std::endl
	     << "SphericalOrbit::compute_freq @ r_apo[2]: model=" 
	     << model->ModelID            << std::endl
	     << "** energy="   << energy  << std::endl
	     << "** pot(max)=" << model->get_pot(xend) << std::endl
	     << "** kappa="    << kappa   << std::endl
	     << "** J="        << kappa*jmax << std::endl
	     << "** xmin="     << xbeg    << std::endl
	     << "** xmax="     << xend    << std::endl
	     << "** cond="     << std::boolalpha << cond << std::endl
	     << "** rc check=" << energy - model->get_pot(r_circ) - jmax*jmax/(2.0*r_circ*r_circ) << std::endl;
	throw std::runtime_error(mesg.str());
      }
    } else {
      r_apo = r_circ;
    }
  }
  else {	 		// Circular radius outside mass distribution
    double JJ = kappa*jmax;
    double r0 = -0.5*model->get_mass(model->get_max_radius())/energy;
    r_apo = r0 + sqrt(r0*r0 - 0.5*JJ*JJ/energy);
				// Newton-Raphson refinement (slow . . .)
    double f, df;
    for (int i=0; i<200; i++) {
      f = 2.0*(energy - model->get_pot(r_apo)) - JJ*JJ/(r_apo*r_apo);
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

  {
    auto [xbeg, xend, cond] = search(denom, ZFRAC*model->get_min_radius(), r_circ);

    if (not cond) {
      // Assume peri is inside the range
      r_peri = 0.99*ZFRAC*model->get_min_radius();

    } else {
      // Find the peri from the bounded interval
      try {
	r_peri = zbrent(denom, 0.999*xbeg, 1.001*xend, tol);
      }
      catch (const std::exception& error) {
	std::ostringstream mesg;
	mesg << "SphericalOrbit::compute_freq: " << error.what() << std::endl
	     << "SphericalOrbit::compute_freq @ r_peri: model=" << model->ModelID
	     << std::endl
	     << "** energy=" << energy     << std::endl
	     << "** kappa="  << kappa      << std::endl
	     << "** J="      << kappa*jmax << std::endl
	     << "** rmin="   << xbeg       << std::endl
	     << "** rmax="   << xend       << std::endl
	     << "** cond="   << std::boolalpha << cond << std::endl;
	
	throw std::runtime_error(mesg.str());
      }
    }
  }

  //  Ok, prepare to do freq. integrals 
  //
  ap = 0.5*(r_apo + r_peri);
  am = 0.5*(r_apo - r_peri);
  sp = ap/(r_apo*r_peri);
  sm = am/(r_apo*r_peri);

  // Do using centered rectangle technique 
  //
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
      msg << "E=" << energy << "  K=" << kappa << 
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
      msg << "\t\tE=" << energy << "  K=" << kappa << 
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
    double denom = 2.0*(energy-ur) - jmax*jmax*kappa*kappa*s*s;
    if (denom>0.0) accum2 += cost/sqrt(denom);
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

// for testing---frequencies in the epicyclic approx 
//
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
double  rombe2(double a, double b, std::function<double(double)>, int n);

double dtp, dtm;

void SphericalOrbit::compute_angles(void)
{
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

  double JJ = kappa * jmax;

  if (freq_defined) {
    ap = 0.5*(r_apo + r_peri);
    am = 0.5*(r_apo - r_peri);
    sp = ap/(r_apo*r_peri);
    sm = am/(r_apo*r_peri);
  }
  else
    compute_freq();


  const double tol = 1.0e-8;       // tolerance for radicand 

  auto fw1 = [&](double t)
  {
    double r  = ap + am*sin(t);
    double ur = model->get_pot(r);
    double radcand = 2.0*(energy-ur)-JJ*JJ/(r*r);

    double eps1 = fabs( 0.5*M_PI + t);
    double eps2 = fabs(-0.5*M_PI + t);

    if (radcand < tol or eps1<1.0e-3 or eps2<1.0e-3) {
      // value near pericenter
      if (t<0) {
	radcand = fabs(JJ*JJ/(r_peri*r_peri*r_peri) - model->get_dpot(r_peri));
	// radcand = JJ*JJ/(r_peri*r_peri*r_peri) - model->get_dpot(r_peri);
	if (radcand < 0.0) {
	  std::ostringstream sout;
	  sout << "fw1: impossible radicand with t=" << t;
	  throw std::runtime_error(sout.str());
	}
	return sqrt(am/radcand);
      }
      // value near apocenter
      else {
	radcand = fabs(model->get_dpot(r_apo) - JJ*JJ/(r_apo*r_apo*r_apo));
	// radcand = model->get_dpot(r_apo) - JJ*JJ/(r_apo*r_apo*r_apo);
	if (radcand < 0.0) {
	  std::ostringstream sout;
	  sout << "fw1: impossible radicand with t=" << t;
	  throw std::runtime_error(sout.str());
	}
	return sqrt(am/radcand);
      }
    }

    return am*cos(t)/sqrt(radcand);
  };

  auto ff = [&](double t)
  {
    double s = sp + sm*sin(t);
    double ur = model->get_pot(1.0/s);
    double radcand = 2.0*(energy-ur) - JJ*JJ*s*s;

    double eps1 = fabs( 0.5*M_PI + t);
    double eps2 = fabs(-0.5*M_PI + t);

    if (radcand < tol or eps1<1.0e-3 or eps2<1.0e-3) {
      // value near apocenter
      if (t<0) {
	radcand = fabs(model->get_dpot(r_apo) - JJ*JJ/(r_apo*r_apo*r_apo));
	// radcand = model->get_dpot(r_apo) - JJ*JJ/(r_apo*r_apo*r_apo);
	if (radcand < 0.0) {
	  std::ostringstream sout;
	  sout << "ff: impossible radicand with t=" << t;
	  throw std::runtime_error(sout.str());
	}
	return sqrt(sm/radcand)/r_apo;
      }
      // value near pericenter
      else {
	radcand = fabs(JJ*JJ/(r_peri*r_peri*r_peri) - model->get_dpot(r_peri));
	// radcand = JJ*JJ/(r_peri*r_peri*r_peri) - model->get_dpot(r_peri);
	if (radcand < 0.0) {
	  std::ostringstream sout;
	  sout << "ff: impossible radicand with t=" << t;
	  throw std::runtime_error(sout.str());
	}
	return sqrt(sm/radcand)/r_peri;
      }
    }

    return sm*cos(t)/sqrt(radcand);
  };

  double accum1 = 0.0;
  double accum2 = 0.0;
  double tl = -0.5*M_PI;
  double sl =  0.5*M_PI;

  for (int i=0; i<recs; i++) {
    double t = dtp + dtm*i/(recs-1);
    double r = ap + am*sin(t);

    if (i>0) accum1 += rombe2(tl, t, fw1, nbsct);

    double arg = (1.0/r - sp)/sm, s;

    if (arg>1.0)       s =  0.5*M_PI;
    else if (arg<-1.0) s = -0.5*M_PI;
    else               s = asin((1.0/r - sp)/sm);

    if (i>0) accum2 += rombe2(sl, s, ff, nbsct);

    angle_grid.t(0, i)     = t;
    angle_grid.w1(0, i)    = freq[0]*accum1;
    angle_grid.dw1dt(0, i) = freq[0]*fw1(t);
    angle_grid.f(0, i)     = freq[1]*accum1 + JJ*accum2;
    angle_grid.r(0, i)     = r;

    tl = t;
    sl = s;
  }


  Eigen::VectorXd work(angle_grid.t.row(0).size());
  
  const double bc = 1.0e30;

  Spline(angle_grid.w1.row(0), angle_grid.t.row(0),     bc, bc, work);
  angle_grid.t.row(1) = work;

  Spline(angle_grid.w1.row(0), angle_grid.dw1dt.row(0), bc, bc, work);
  angle_grid.dw1dt.row(1) = work;

  Spline(angle_grid.w1.row(0), angle_grid.f.row(0),     bc, bc, work);
  angle_grid.f.row(1) = work;

  Spline(angle_grid.w1.row(0), angle_grid.r.row(0),     bc, bc, work);
  angle_grid.r.row(1) = work;

  angle_grid.num = recs;

  angle_defined = true;
}

void SphericalOrbit::compute_angles_old(void)
{
  double accum1,accum2,r, s, t, sl, tl;
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

  
  double JJ = kappa * jmax;

  if (freq_defined) {
    ap = 0.5*(r_apo + r_peri);
    am = 0.5*(r_apo - r_peri);
    sp = ap/(r_apo*r_peri);
    sm = am/(r_apo*r_peri);
  }
  else
    compute_freq();


  const double tol = 1.0e-8;       // tolerance for radicand 

  auto fw1 = [&](double t)
  {
    double r  = ap + am*sin(t);
    double ur = model->get_pot(r);
    double radcand = 2.0*(energy-ur)-JJ*JJ/(r*r);

    double eps1 = fabs( 0.5*M_PI + t);
    double eps2 = fabs(-0.5*M_PI + t);

    if (radcand < tol or eps1<1.0e-3 or eps2<1.0e-3) {
      // value near pericenter
      if (t<0) {
	radcand = fabs(JJ*JJ/(r_peri*r_peri*r_peri) - model->get_dpot(r_peri));
	// radcand = JJ*JJ/(r_peri*r_peri*r_peri) - model->get_dpot(r_peri);
	if (radcand < 0.0) {
	  std::ostringstream sout;
	  sout << "fw1: impossible radicand with t=" << t;
	  throw std::runtime_error(sout.str());
	}
	return sqrt(am/radcand);
      }
      // value near apocenter
      else {
	radcand = fabs(model->get_dpot(r_apo) - JJ*JJ/(r_apo*r_apo*r_apo));
	// radcand = model->get_dpot(r_apo) - JJ*JJ/(r_apo*r_apo*r_apo);
	if (radcand < 0.0) {
	  std::ostringstream sout;
	  sout << "fw1: impossible radicand with t=" << t;
	  throw std::runtime_error(sout.str());
	}
	return sqrt(am/radcand);
      }
    }

    return am*cos(t)/sqrt(radcand);
  };

  auto ff = [&](double t)
  {
    double s = sp + sm*sin(t);
    double ur = model->get_pot(1.0/s);
    double radcand = 2.0*(energy-ur) - JJ*JJ*s*s;

    double eps1 = fabs( 0.5*M_PI + t);
    double eps2 = fabs(-0.5*M_PI + t);

    if (radcand < tol or eps1<1.0e-3 or eps2<1.0e-3) {
      // value near apocenter
      if (t<0) {
	radcand = fabs(model->get_dpot(r_apo) - JJ*JJ/(r_apo*r_apo*r_apo));
	// radcand = model->get_dpot(r_apo) - JJ*JJ/(r_apo*r_apo*r_apo);
	if (radcand < 0.0) {
	  std::ostringstream sout;
	  sout << "ff: impossible radicand with t=" << t;
	  throw std::runtime_error(sout.str());
	}
	return sqrt(sm/radcand)/r_apo;
      }
      // value near pericenter
      else {
	radcand = fabs(JJ*JJ/(r_peri*r_peri*r_peri) - model->get_dpot(r_peri));
	// radcand = JJ*JJ/(r_peri*r_peri*r_peri) - model->get_dpot(r_peri);
	if (radcand < 0.0) {
	  std::ostringstream sout;
	  sout << "ff: impossible radicand with t=" << t;
	  throw std::runtime_error(sout.str());
	}
	return sqrt(sm/radcand)/r_peri;
      }
    }

    return sm*cos(t)/sqrt(radcand);
  };

  accum1 = 0.0;
  accum2 = 0.0;
  tl = -0.5*M_PI;
  sl =  0.5*M_PI;
  for (int i=0; i<recs; i++) {
    t = dtp + dtm*Gkn.knot(i);
    r = ap + am*sin(t);
    accum1 += rombe2(tl, t, fw1, nbsct);
    s = asin((1.0/r - sp)/sm);
    accum2 += rombe2(sl, s, ff, nbsct);

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

//
// epicyclic test 
//
void SphericalOrbit::compute_angles_epi(void)
{
  double r,t,a,a2,fac,ur;
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

  if (not freq_defined) compute_freq_epi();

  double JJ = jmax*kappa;
  ur = model->get_pot(r_circ);
  a2 = 2.0*(energy-0.5*JJ*JJ/(r_circ*r_circ)-ur)/(freq[0]*freq[0]);
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


double SphericalOrbit::pot_trans(int l1, int l2,
				 std::function<double(double)> func)
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

