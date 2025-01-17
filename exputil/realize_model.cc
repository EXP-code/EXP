/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Generate a point for a 3D spherical model
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/20/91
 *
 ***************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <memory>
#include <cmath>

#include <localmpi.H>
#include <massmodel.H>
#include <interp.H>

#ifdef DEBUG
#include <orbit.H>
static SphericalOrbit orb;
#endif

#include <libvars.H>
using namespace __EXP__;

// Multimass Monte-Carlo rejection type.  Fixed-radius is the
// probabilistically correct choice to recover the number density
// profile but will be slow for very steep number densities.  If in
// doubt, use fixed-radius.
//
#define FIXED_RADIUS 1

bool     AxiSymModel::gen_EJ    = true;
int      AxiSymModel::numr      = 50;
int      AxiSymModel::numj      = 50;
int      AxiSymModel::gen_N     = 400;
int      AxiSymModel::gen_E     = 400;
int      AxiSymModel::gen_K     = 200;
double   AxiSymModel::gen_tolE  = 0.01;
double   AxiSymModel::gen_tolK  = 0.02;
double   AxiSymModel::gen_rmin  = 0.0;
int      AxiSymModel::gen_logr  = 1;
double   AxiSymModel::gen_kmin  = 0.0;
int      AxiSymModel::gen_itmax = 20000;

const bool verbose = true;
const double ftol = 0.01;

AxiSymModel::PSret AxiSymModel::gen_point_2d()
{
  if (!dist_defined) {
    throw std::runtime_error("AxiSymModel: must define distribution before realizing!");
  }

  int it;			// Iteration counter
  double r=0.0, pot, vmax, vr=0.0, vt=0.0, eee, fmax, vv;
  double phi=0.0, sinp, cosp;
  double T, w1;
  double tol = 1.0e-5;
  double rmin = max<double>(get_min_radius(), gen_rmin);

  double Emin = get_pot(rmin);
  double Emax = get_pot(get_max_radius());

  if (gen_EJ) {

    if (gen_firstime) {

      gen_orb = SphericalOrbit(shared_from_this());

      double dx = (Emax - Emin - 2.0*tol)/(numr-1);
      double dy = (1.0 - gen_kmin - 2.0*tol)/(numj-1);

#ifdef DEBUG
      std::cout << "gen_point_2d[" << ModelID << "]: " << get_max_radius() << std::endl;
#endif
      gen_fomax = 0.0;

      for (int j=0; j<numr; j++) {
	double xxx = Emin + tol + dx*j;

	for (int k=0; k<numj; k++) {
	  double yyy = gen_kmin + tol + dy*k;

	  gen_orb.new_orbit(xxx, yyy);
	    
	  double zzz = distf(xxx, gen_orb.AngMom())*gen_orb.Jmax()
	    /gen_orb.get_freq(1);
	  gen_fomax = zzz>gen_fomax ? zzz : gen_fomax;
	}
      }

      gen_firstime = false;
    }

    
    // Trial velocity point
    //
    for (it=0; it<gen_itmax; it++) {

      double xxx = Emin + tol + (Emax - Emin - 2.0*tol)*Unit(random_gen);
      double yyy = gen_kmin + tol + (1.0 - gen_kmin - 2.0*tol)*Unit(random_gen);

      gen_orb.new_orbit(xxx, yyy);

      double zzz = distf(xxx, gen_orb.AngMom()) * gen_orb.Jmax()/gen_orb.get_freq(1);

      if (Unit(random_gen) > zzz/gen_fomax ) continue;

      w1 = 2.0*M_PI*Unit(random_gen);
      T = w1/gen_orb.get_freq(1);
      
      r = gen_orb.get_angle(6, T);
      phi = 2.0*M_PI*Unit(random_gen) - gen_orb.get_angle(5, T);

      pot = get_pot(r);

      vt = gen_orb.AngMom()/r;
      vv = 2.0*(xxx  - pot) - vt*vt;
      if (vv<=0.0)
	vr = 0.0;
      else
	vr = sqrt(vv);

      if (w1 > M_PI) vr *= -1.0;
      
      break;
    }

  }
  else {

    if (gen_firstime) {

      double tol = 1.0e-2;
      double dx = (1.0 - 2.0*tol)/(numr-1);
      double dy = (1.0 - 2.0*tol)/(numj-1);
      double dr;

      gen_mass.resize(gen_N);
      gen_rloc.resize(gen_N);
      gen_fmax.resize(gen_N);

#ifdef DEBUG
      std::cout << "gen_point_2d[" << ModelID << "]: " << get_max_radius() << std::endl;
#endif

      if (rmin <= 1.0e-16) gen_logr = 0;

      if (gen_logr)
	dr = (log(get_max_radius()) - log(rmin))/(gen_N-1);
      else
	dr = (get_max_radius() - rmin)/(gen_N-1);

      for (int i=0; i<gen_N; i++) {

	if (gen_logr) {
	  gen_rloc[i] = log(rmin + dr*i);
	  r = exp(gen_rloc[i]);
	}
	else {
	  gen_rloc[i] = rmin + dr*i;
	  r = gen_rloc[i];
	}

	gen_mass[i] = get_mass(r);
	
	pot = get_pot(r);
	vmax = sqrt(2.0*fabs(Emax - pot));

	fmax = 0.0;
	for (int j=0; j<numr; j++) {
	  double xxx = sqrt(tol + dx*j);

	  for (int k=0; k<numj; k++) {
	    double yyy = 0.5*M_PI*(tol + dy*k);

	    vr = vmax*xxx*cos(yyy);
	    vt = vmax*xxx*sin(yyy);
	    eee = pot + 0.5*(vr*vr + vt*vt);

	    double zzz = distf(eee, gen_rloc[i]*vt);
	    fmax = zzz>fmax ? zzz : fmax;
	  }
	}
	gen_fmax[i] = fmax*(1.0 + ftol);

      }
      gen_firstime = false;
    }

    r = odd2(Unit(random_gen)*gen_mass[gen_N-1], gen_mass, gen_rloc, 0);
    fmax = odd2(r, gen_rloc, gen_fmax, 1);
    if (gen_logr) r = exp(r);

    pot = get_pot(r);
    vmax = sqrt(2.0*fabs(Emax - pot));

    // Trial velocity point
    //
    for (it=0; it<gen_itmax; it++) {

      double xxx = sqrt(Unit(random_gen));
      double yyy = 0.5*M_PI*Unit(random_gen);

      vr = vmax*xxx*cos(yyy);
      vt = vmax*xxx*sin(yyy);
      eee = pot + 0.5*(vr*vr + vt*vt);

      if (Unit(random_gen) > distf(eee, r*vt)/fmax ) continue;

      if (Unit(random_gen)<0.5) vr *= -1.0;
    
      phi = 2.0*M_PI*Unit(random_gen);

      break;
    }

  }

                
  Eigen::VectorXd out(7);

  static unsigned totcnt = 0, toomany = 0;
  totcnt++;
  
  if (it==gen_itmax) {
    if (verbose)
      std::cerr << "Velocity selection failed [" << myid << "]: r="
		<< std::setw(12) << r
		<< std::setw(12) << fmax
		<< " %=" << std::setw(12)
		<< static_cast<double>(++toomany)/totcnt << std::endl;
    out.setZero();
    return {out, 1};
  }

  cosp = cos(phi);
  sinp = sin(phi);

  out[0] = 1.0;
  out[1] = r * cosp;
  out[2] = r * sinp;
  out[3] = 0.0;
  out[4] = vr * cosp - vt * sinp;
  out[5] = vr * sinp + vt * cosp;
  out[6] = 0.0;
    
  return {out, 0};
}


AxiSymModel::PSret AxiSymModel::gen_point_2d(double r)
{
  if (!dist_defined) {
    throw std::runtime_error("AxiSymModel: must define distribution before realizing!");
  }

  int it;			// Iteration counter
  double pot, vmax, vr=0.0, vt=0.0, eee, fmax;
  double phi=0.0, sinp, cosp;
  double rmin = max<double>(get_min_radius(), gen_rmin);
  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {

    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);

    gen_mass.resize(gen_N);
    gen_rloc.resize(gen_N);
    gen_fmax.resize(gen_N);

#ifdef DEBUG
    std::cout << "gen_point_2d[" << ModelID << "]: " << rmin
	 << ", " << get_max_radius() << std::endl;
#endif

    double dr = (get_max_radius() - rmin)/gen_N;

    for (int i=0; i<gen_N; i++) {
      gen_rloc[i] = rmin + dr*i;
      gen_mass[i] = get_mass(gen_rloc[i]);

      pot = get_pot(gen_rloc[i]);
      vmax = sqrt(2.0*fabs(Emax - pot));

      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	double xxx = sqrt(tol + dx*j);

	for (int k=0; k<numj; k++) {
	  double yyy = 0.5*M_PI*(tol + dy*k);

	  vr = vmax*xxx*cos(yyy);
	  vt = vmax*xxx*sin(yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  double zzz = distf(eee, gen_rloc[i]*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*(1.0 + ftol);

    }
    gen_firstime = false;
  }

  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  pot = get_pot(r);
  vmax = sqrt(2.0*fabs(Emax - pot));
  
  // Trial velocity point
  //
  for (it=0; it<gen_itmax; it++) {

    double xxx = sqrt(Unit(random_gen));
    double yyy = 0.5*M_PI*Unit(random_gen);

    vr = vmax*xxx*cos(yyy);
    vt = vmax*xxx*sin(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (Unit(random_gen) > distf(eee, r*vt)/fmax ) continue;
    
    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    phi = 2.0*M_PI*Unit(random_gen);

    break;
  }

  
  Eigen::VectorXd out(7);

  static unsigned totcnt = 0, toomany = 0;
  totcnt++;

  if (it==gen_itmax) {
    if (verbose)
      std::cerr << "Velocity selection failed [" << myid << "]: r="
		<< std::setw(12) << r
		<< std::setw(12) << fmax
		<< " %=" << std::setw(12)
		<< static_cast<double>(++toomany)/totcnt << std::endl;
    out.setZero();
    return {out, 1};
  }

  cosp = cos(phi);
  sinp = sin(phi);

  out[0] = 1.0;
  out[1] = r * cosp;
  out[2] = r * sinp;
  out[3] = 0.0;
  out[4] = vr * cosp - vt * sinp;
  out[5] = vr * sinp + vt * cosp;
  out[6] = 0.0;
    
  return {out, 0};
}


AxiSymModel::PSret AxiSymModel::gen_point_3d()
{
  if (!dist_defined) {
    throw std::runtime_error("AxiSymModel: must define distribution before realizing!");
  }

#ifdef DEBUG
  static ofstream tout("gen3d.ktest");
#endif

  double r, pot, vmax, vr=0.0, vt, eee, vt1=0.0, vt2=0.0, fmax;
  double phi, sint, cost, sinp, cosp;

  double rmin = max<double>(get_min_radius(), gen_rmin);
  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {

#ifdef DEBUG
    orb = SphericalOrbit(this);
#endif

    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);
    double dr;

    gen_mass.resize(gen_N);
    gen_rloc.resize(gen_N);
    gen_fmax.resize(gen_N);

    if (rmin <= 1.0e-16) gen_logr = 0;
    
    if (gen_logr)
      dr = (log(get_max_radius()) - log(rmin))/(gen_N-1);
    else
      dr = (get_max_radius() - rmin)/(gen_N-1);


    for (int i=0; i<gen_N; i++) {

      if (gen_logr) {
	gen_rloc[i] = log(rmin) + dr*i;
	r = exp(gen_rloc[i]);
      }
      else {
	gen_rloc[i] = rmin + dr*i;
	r = gen_rloc[i];
      }

      gen_mass[i] = get_mass(r);

      pot = get_pot(r);
      vmax = sqrt(2.0*fabs(Emax - pot));

      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	double xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  double yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  double zzz = distf(eee, r*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*(1.0 + ftol);

    }

    // Debug
    //
    if (myid==0) {
      std::ofstream test("test.grid");
      if (test) {

	test << "# Rmin=" << rmin
	     << "  Rmax=" << get_max_radius()
	     << std::endl;
	
	for (int i=0; i<gen_N; i++) {
	  test << std::setw(15) << gen_rloc[i]
	       << std::setw(15) << gen_mass[i]
	       << std::setw(15) << gen_fmax[i]
	       << std::endl;
	}
      }
    }

    gen_firstime = false;
  }

  r = odd2(Unit(random_gen)*gen_mass[gen_N-1], gen_mass, gen_rloc, 0);
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  if (gen_logr) r = exp(r);
  
  pot = get_pot(r);
  vmax = sqrt(2.0*fabs(Emax - pot));

  // Trial velocity point
  //
  int it;			// Iteration counter
  double angmom[3];

  for (it=0; it<gen_itmax; it++) {

    double xxx = -2.0*cos(acos(Unit(random_gen))/3.0 - 2.0*M_PI/3.0);
    double yyy = (1.0 - xxx*xxx)*Unit(random_gen);

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    // Debug
    /*
    if (sqrt(vr*vr + vt*vt)>0.99*vmax) {
      std::cout << "Check: df val = " << distf(eee, r*vt)/fmax 
	   << "  v/vmax = " << sqrt(vr*vr + vt*vt)/vmax 
	   << "  eee = " << eee
	   << std::endl;
    }
    */

    if (Unit(random_gen) > distf(eee, r*vt)/fmax ) continue;

    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    // Orientation of tangential velocity vector
    //
    double azi = 2.0*M_PI*Unit(random_gen);
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  Eigen::VectorXd out(7);

  static unsigned totcnt = 0, toomany = 0;
  totcnt++;

  if (it==gen_itmax) {
    if (verbose)
      std::cerr << "Velocity selection failed [" << myid << "]: r="
		<< std::setw(12) << r
		<< std::setw(12) << fmax
		<< " %=" << std::setw(12)
		<< static_cast<double>(++toomany)/totcnt << std::endl;
    out.setZero();
    return {out, 1};
  }

  if (Unit(random_gen)>=0.5) vr *= -1.0;

  phi = 2.0*M_PI*Unit(random_gen);
  cost = 2.0*(Unit(random_gen) - 0.5);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  out[0] = 1.0;
  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    
#ifdef DEBUG
  eee = pot + 0.5*(out[4]*out[4]+out[5]*out[5]+out[6]*out[6]);
  orb.new_orbit(eee, 0.5);
  angmom[0] = out[2]*out[6] - out[3]*out[5];
  angmom[1] = out[3]*out[4] - out[1]*out[6];
  angmom[2] = out[1]*out[5] - out[2]*out[4];
  tout << std::setw(15) << eee
       << std::setw(15) << sqrt(angmom[0]*angmom[0]+
				angmom[1]*angmom[1]+
			   angmom[2]*angmom[2])/orb.Jmax()
       << std::setw(15) << r
       << std::setw(15) << sqrt(out[4]*out[4]+out[5]*out[5]+out[6]*out[6])
       << std::endl;
#endif


  return {out, 0};
}


AxiSymModel::PSret AxiSymModel::gen_point_3d(double r, double theta, double phi)
{
  if (!dist_defined) {
    throw std::runtime_error("AxiSymModel: must define distribution before realizing!");
  }

  int it;			// Iteration counter
  double pot, vmax, vr=0.0, vt=0.0, eee, fmax;
  double rmin = max<double>(get_min_radius(), gen_rmin);
  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {

    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);

    gen_mass.resize(gen_N);
    gen_rloc.resize(gen_N);
    gen_fmax.resize(gen_N);

#ifdef DEBUG
    std::cout << "gen_point_3d[" << ModelID << "]: " << rmin
	 << ", " << get_max_radius() << std::endl;
#endif

    double dr = (get_max_radius() - rmin)/gen_N;

    for (int i=0; i<gen_N; i++) {
      gen_rloc[i] = rmin + dr*i;
      gen_mass[i] = get_mass(gen_rloc[i]);

      pot = get_pot(gen_rloc[i]);
      vmax = sqrt(2.0*fabs(Emax - pot));

      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	double xxx = sqrt(tol + dx*j);

	for (int k=0; k<numj; k++) {
	  double yyy = 0.5*M_PI*(tol + dy*k);

	  vr = vmax*xxx*cos(yyy);
	  vt = vmax*xxx*sin(yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  double zzz = distf(eee, gen_rloc[i]*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*(1.0 + ftol);

    }
    gen_firstime = false;
  }

  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  pot  = get_pot(r);
  vmax = sqrt(2.0*fabs(Emax - pot));
  
  // Trial velocity point
  //
  for (it=0; it<gen_itmax; it++) {

    double xxx = pow(Unit(random_gen), 0.3333333333333333);
    double yyy = 0.5*M_PI*Unit(random_gen);

    vr = vmax*xxx*cos(yyy);
    vt = vmax*xxx*sin(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (Unit(random_gen) > distf(eee, r*vt)/fmax ) continue;
    
    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    break;
  }

  Eigen::VectorXd out(7);

  if (std::isnan(vr) or std::isnan(vt)) {
    if (verbose) std::cout << "NaN found in AxiSymModel::gen_point_3d with r="
			   << r << " theta=" << theta << " phi=" << phi
			   << std::endl;
    out.setZero();
    return {out, 1};
  }
  
  static unsigned totcnt = 0, toomany = 0;
  totcnt++;

  if (it==gen_itmax) {
    if (verbose)
      std::cerr << "Velocity selection failed [" << myid << "]: r="
		<< std::setw(12) << r
		<< std::setw(12) << fmax
		<< " %=" << std::setw(12)
		<< static_cast<double>(++toomany)/totcnt << std::endl;
    out.setZero();
    return {out, 1};
  }

  double tv  = 2.0*M_PI*Unit(random_gen);
  double vth = vt*cos(tv);
  double vph = vt*sin(tv);

  double cost = cos(theta);
  double sint = sin(theta);

  double cosp = cos(phi);
  double sinp = sin(phi);

  out[0] = 1.0;
  out[1] = r*sint*cosp;
  out[2] = r*sint*sinp;
  out[3] = r*cost;
  out[4] = vr*sint*cosp - vph*sinp + vth*cost*cosp;
  out[5] = vr*sint*sinp + vph*cosp + vth*cost*sinp;
  out[6] = vr*cost - vth*sint;
    
  return {out, 0};
}


AxiSymModel::PSret AxiSymModel::gen_point_3d(double Emin, double Emax, 
				double Kmin, double Kmax)
{
  if (!dist_defined) {
    throw std::runtime_error("AxiSymModel: must define distribution before realizing!");
  }

#ifdef DEBUG
  double angmom[3];
  static ofstream tout("gen3d.Etest");
#endif

  double r, vr, vt, vt1, vt2, E, K, J, jmax, w1t, eee, pot;
  double phi, sint, cost, sinp, cosp;
  double rmin = max<double>(get_min_radius(), gen_rmin);
  double rmax = get_max_radius();

  if (gen_firstime_E) {

#ifdef DEBUG
    orb = SphericalOrbit(shared_from_this());
#endif
    gen_orb = SphericalOrbit(shared_from_this());

    Emin_grid = get_pot(rmin)*(1.0 - gen_tolE);
    Emax_grid = get_pot(rmax)*(1.0 + gen_tolE);

    dEgrid = (Emax_grid - Emin_grid)/(gen_E-1);
    dKgrid = (1.0 - 2.0*gen_tolK)/(gen_K-1);

    for (int i=0; i<gen_E; i++) Egrid.push_back(Emin_grid + dEgrid*i);
    for (int i=0; i<gen_K; i++) Kgrid.push_back(gen_tolK  + dKgrid*i);

    double ans, K, angles = 2.0*pow(2.0*M_PI, 3.0), tfac;
    vector<double> tmass;
    ANGLE_GRID *agrid;

    for (int i=0; i<gen_E; i++) {
      ans = 0.0;

      vector<WRgrid> wrvec;
      for (int j=0; j<gen_K; j++) {
				// Trapezoidal rule factor
	if (j==0 || j==gen_K) tfac = 0.5;
	else tfac = 1.0;

	K = gen_tolK + dKgrid*j;
	gen_orb.new_orbit(Egrid[i], K);
	ans += K*gen_orb.Jmax()*gen_orb.Jmax() /
	  gen_orb.get_freq(1) * distf(Egrid[i], gen_orb.AngMom()) * tfac;

	WRgrid wr;
	agrid = gen_orb.get_angle_grid();

	for (int n=0; n<agrid->num; n++) {
	  wr.w1.push_back(agrid->w1(0, n));
	  wr.r.push_back(agrid->r(0, n));
	}
	wrvec.push_back(wr);

      } // End Kappa loop

      Jmax.push_back(gen_orb.Jmax());
      tmass.push_back(ans*angles*dKgrid*dEgrid);
      Rgrid.push_back(wrvec);

    } // Energy E loop

    EgridMass.push_back(0.0);
    for (int i=1; i<gen_E; i++) 
      EgridMass.push_back(EgridMass[i-1]+tmass[i]);

    // Debug
    
    ofstream test("test.grid");
    if (test) {

      test << "# Emin=" << Emin_grid
	   << "  Emax=" << Emax_grid
	   << std::endl;

      for (int i=0; i<gen_E; i++) {
	test << std::setw(15) << Egrid[i]
	     << std::setw(15) << EgridMass[i]
	     << std::setw(15) << Jmax[i]
	     << std::endl;
      }
    }

    gen_firstime_E = false;

#ifdef DEBUG
    tout << "# Emin=" << Emin_grid << "  Emax=" << Emax_grid
	 << "  Mass frac=" << 
      ( odd2(Emax, Egrid, EgridMass, 1) - 
	odd2(Emin, Egrid, EgridMass, 1) )/EgridMass[gen_E-1] << std::endl;
#endif    
  }

  // Enforce limits
  //
  Emin = max<double>(Emin, Emin_grid);
  Emin = min<double>(Emin, Emax_grid);
  Emax = max<double>(Emax, Emin_grid);
  Emax = min<double>(Emax, Emax_grid);

  double Mmin = odd2(Emin, Egrid, EgridMass, 1);
  double Mmax = odd2(Emax, Egrid, EgridMass, 1);
  double kmin = max<double>(Kmin, gen_tolK);
  double kmax = min<double>(Kmax, 1.0 - gen_tolK);

  E = odd2(Mmin + (Mmax-Mmin)*Unit(random_gen), EgridMass, Egrid, 0);
  K = sqrt(kmin*kmin + (kmax*kmax - kmin*kmin)*Unit(random_gen));

  int indxE = int( (E - Emin_grid)/dEgrid );
  int indxK = int( (K - gen_tolK)/dKgrid );
  
  indxE = max<int>(0, min<int>(gen_E-2, indxE));
  indxK = max<int>(0, min<int>(gen_K-2, indxK));

  double cE[2], cK[2];

  cE[1] = (E - (Emin_grid + dEgrid*indxE))/dEgrid;
  cE[0] = 1.0 - cE[1];

  cK[1] = (K - (gen_tolK + dKgrid*indxK))/dKgrid;
  cK[0] = 1.0 - cK[1];


  r = 0.0;
  J = 0.0;
  jmax = 0.0;
  w1t = M_PI*Unit(random_gen);

  for (int ie=0; ie<2; ie++) {
    J += cE[ie]*Jmax[indxE+ie] * K;
    jmax += cE[ie]*Jmax[indxE+ie];
    for (int ik=0; ik<2; ik++) {
      r += cE[ie]*cK[ik]*odd2(w1t, Rgrid[indxE+ie][indxK+ik].w1, 
			     Rgrid[indxE+ie][indxK+ik].r, 0);
    }
  }
      
  pot = get_pot(r);
  vt = J/r;
  int ierr = 0;

  // Interpolation check (should be rare).  Error condition is set.
  //
  if (2.0*(E - pot) - vt*vt < 0.0) {
    ierr = 1;
    if (E < pot) E = pot;
    vt = sqrt(E - pot);
  }
  vr = sqrt( 2.0*(E - pot) - vt*vt );

  if (Unit(random_gen)<0.5) vr *= -1.0;
    
  // Orientation of tangential velocity vector
  //
  double azi = 2.0*M_PI*Unit(random_gen);
  vt1 = vt*cos(azi);
  vt2 = vt*sin(azi);

  Eigen::VectorXd out(7);

  phi = 2.0*M_PI*Unit(random_gen);
  cost = 2.0*(Unit(random_gen) - 0.5);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  out[0] = (Mmax - Mmin)/EgridMass[gen_E-1];
  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    

#ifdef DEBUG
  if (ierr) return {out, ierr};

  eee = pot + 0.5*(out[4]*out[4]+out[5]*out[5]+out[6]*out[6]);
  orb.new_orbit(eee, 0.5);
  angmom[0] = out[2]*out[6] - out[3]*out[5];
  angmom[1] = out[3]*out[4] - out[1]*out[6];
  angmom[2] = out[1]*out[5] - out[2]*out[4];
  tout << std::setw(15) << eee
       << std::setw(15) << E
       << std::setw(15) << sqrt(angmom[0]*angmom[0]+
			   angmom[1]*angmom[1]+
			   angmom[2]*angmom[2])/orb.Jmax()
       << std::setw(15) << K
       << std::setw(15) << r
       << std::setw(15) << sqrt(out[4]*out[4]+out[5]*out[5]+out[6]*out[6])
       << std::setw(15) << orb.Jmax()
       << std::setw(15) << jmax
       << std::endl;
#endif

  return {out, ierr};
}



AxiSymModel::PSret AxiSymModel::gen_point_jeans_3d()
{
  double r, d, vr, vt, vt1, vt2, vv, vtot;
  double phi, sint, cost, sinp, cosp;
  double rmin = max<double>(get_min_radius(), gen_rmin);

  if (gen_firstime_jeans) {
    double dr;

    gen_mass.resize(gen_N);
    gen_rloc.resize(gen_N);
    gen_fmax.resize(gen_N);
    Eigen::VectorXd work(gen_N);
    Eigen::VectorXd work2(gen_N);

    if (rmin <= 1.0e-16) gen_logr = 0;
    
    if (gen_logr)
      dr = (log(get_max_radius()) - log(rmin))/(gen_N-1);
    else
      dr = (get_max_radius() - rmin)/(gen_N-1);


    for (int i=0; i<gen_N; i++) {

      if (gen_logr) {
	gen_rloc[i] = log(rmin) + dr*i;
	r = exp(gen_rloc[i]);
      }
      else {
	gen_rloc[i] = rmin + dr*i;
	r = gen_rloc[i];
      }

      gen_mass[i] = get_mass(r);
      work[i] = get_dpot(r) * get_density(r);
      if (gen_logr) work[i] *= r;
    }

    Trapsum(gen_rloc, work, work2);

    for (int i=0; i<gen_N; i++)
      gen_fmax[i] = 3.0*(work2[gen_N-1] - work2[i]);


    // Debug
    //
    std::ofstream test("test.grid");
    if (test) {

      test << "# [Jeans] Rmin=" << rmin
	   << "  Rmax=" << get_max_radius()
	   << std::endl;

      for (int i=0; i<gen_N; i++) {
	if (gen_logr) r = exp(gen_rloc[i]);
	else r = gen_rloc[i];

	d = get_density(r);
	if (d>0.0 && gen_fmax[i]>=0.0)
	  vtot = sqrt(gen_fmax[i]/d);
	else
	  vtot = 0.0;

	test << std::setw(15) << gen_rloc[i]
	     << std::setw(15) << gen_mass[i]
	     << std::setw(15) << gen_fmax[i]
	     << std::setw(15) << vtot
	     << std::endl;
      }
    }

    gen_firstime_jeans = false;
  }

  r = odd2(Unit(random_gen)*gen_mass[gen_N-1], gen_mass, gen_rloc, 0);
  vv = odd2(r, gen_rloc, gen_fmax);

  if (gen_logr) r = exp(r);

  d = get_density(r);
  if (d>0.0 && vv >=0.0)
    vtot = sqrt(vv/d);
  else
    vtot = 0.0;
    
  
  double xxx = -2.0*cos(acos(Unit(random_gen))/3.0 - 2.0*M_PI/3.0);
  double yyy = (1.0 - xxx*xxx)*Unit(random_gen);

  vr = vtot*xxx;
  vt = vtot*sqrt(yyy);

  // Orientation of tangential velocity vector
  //
  double azi = 2.0*M_PI*Unit(random_gen);
  vt1 = vt*cos(azi);
  vt2 = vt*sin(azi);

  Eigen::VectorXd out(7);

  if (Unit(random_gen)>=0.5) vr *= -1.0;

  phi = 2.0*M_PI*Unit(random_gen);
  cost = 2.0*(Unit(random_gen) - 0.5);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  out[0] = 1.0;
  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    
  return {out, 0};
}

int AxiSymModel::gen_velocity(double* pos, double* vel)
{
  if (dof()!=3)
      bomb( "AxiSymModel: gen_velocity only implemented for 3d model!" );
  
  if (!dist_defined) {
    throw std::runtime_error("AxiSymModel: must define distribution before realizing!");
  }

  double r, pot, vmax, vr=0.0, vt, eee, vt1=0.0, vt2=0.0, fmax;
  double phi, sint, cost, sinp, cosp;
  double rmin = max<double>(get_min_radius(), gen_rmin);

  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {

    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);
    double dr;

    gen_mass.resize(gen_N);
    gen_rloc.resize(gen_N);
    gen_fmax.resize(gen_N);

    if (rmin <= 1.0e-16) gen_logr = 0;
    
    if (gen_logr)
      dr = (log(get_max_radius()) - log(rmin))/(gen_N-1);
    else
      dr = (get_max_radius() - rmin)/(gen_N-1);


    for (int i=0; i<gen_N; i++) {

      if (gen_logr) {
	gen_rloc[i] = log(rmin) + dr*i;
	r = exp(gen_rloc[i]);
      }
      else {
	gen_rloc[i] = rmin + dr*i;
	r = gen_rloc[i];
      }

      gen_mass[i] = get_mass(r);

      pot = get_pot(r);
      vmax = sqrt(2.0*fabs(Emax - pot));

      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	double xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  double yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  double zzz = distf(eee, r*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*(1.0 + ftol);

    }
    gen_firstime = false;

#ifdef DEBUG
    ofstream tout("test.dfgrid");
    tout << "# Rmin=" << rmin 
	 << "  Rmax=" << get_max_radius() << std::endl;
    for (int i=0; i<gen_N; i++)
      tout << std::setw(18) << gen_rloc[i]
	   << std::setw(18) << gen_mass[i]
	   << std::setw(18) << gen_fmax[i]
	   << std::endl;
#endif
  }

  r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
  if (gen_logr) r = log(r);
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  if (gen_logr) r = exp(r);
  
  pot = get_pot(r);
  vmax = sqrt(2.0*fabs(Emax - pot));

  // Trial velocity point
  //
  int it;			// Iteration counter

  for (it=0; it<gen_itmax; it++) {

    double xxx = -2.0*cos(acos(Unit(random_gen))/3.0 - 2.0*M_PI/3.0);
    double yyy = (1.0 - xxx*xxx)*Unit(random_gen);

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (Unit(random_gen) > distf(eee, r*vt)/fmax ) continue;

    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    // Orientation of tangential velocity vector
    //
    double azi = 2.0*M_PI*Unit(random_gen);
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  static unsigned totcnt = 0, toomany = 0;
  totcnt++;

  if (it==gen_itmax) {
    if (verbose)
      std::cerr << "Velocity selection failed [" << myid << "]: r="
		<< std::setw(12) << r
		<< std::setw(12) << fmax
		<< " %=" << std::setw(12)
		<< static_cast<double>(++toomany)/totcnt << std::endl;
    return 1;
  }

  if (Unit(random_gen)>=0.5) vr *= -1.0;

  phi = atan2(pos[1], pos[0]);
  cost = pos[2]/(r+1.0e-18);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  vel[0] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  vel[1] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  vel[2] = vr * cost      - vt1 * sint;

  return 0;
}

AxiSymModel::PSret SphericalModelMulti::gen_point()
{
  if (!real->dist_defined || !fake->dist_defined) {
    throw std::runtime_error("SphericalModelMulti: input distribution functions must be defined before realizing!");
  }

  double r, pot, vmax;
  double vr=0.0, vt=0.0, eee=0.0, vt1=0.0, vt2=0.0, fmax, emax;
  double mass, phi, sint, cost, sinp, cosp;

  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {

    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);
    double dr;

    gen_mass.resize(gen_N);
    gen_rloc.resize(gen_N);
    gen_fmax.resize(gen_N);

    gen_mass.setZero();
    gen_rloc.setZero();
    gen_fmax.setZero();

    std::vector<int> ibeg(numprocs);
    std::vector<int> iend(numprocs);

    int dN = gen_N/numprocs;
    for (int n=0; n<numprocs; n++) {
      ibeg[n] = dN*n;
      iend[n] = dN*(n+1);
    }
    iend[numprocs-1] = gen_N;

    std::vector<double> gen_emax(gen_N, 0.0), gen_vmax(gen_N, 0.0);

    if (rmin_gen <= 1.0e-16) gen_logr = 0;
    
    if (gen_logr)
      dr = (log(rmax_gen) - log(rmin_gen))/(gen_N-1);
    else
      dr = (rmax_gen - rmin_gen)/(gen_N-1);

    for (int i=ibeg[myid]; i<iend[myid]; i++) {

      if (gen_logr) {
	gen_rloc[i] = log(rmin_gen) + dr*i;
	r = exp(gen_rloc[i]);
      }
      else {
	gen_rloc[i] = rmin_gen + dr*i;
	r = gen_rloc[i];
      }

      gen_mass[i] = fake->get_mass(r);

      pot  = get_pot(r);
      vmax = sqrt(2.0*fabs(Emax - pot));

      emax = pot;
      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	double xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  double yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);
	  
	  double zzz = fake->distf(eee, r*vt);
	  if (zzz>fmax) {
	    emax = eee;
	    fmax = zzz;
	  }
	}
      }
      gen_emax[i] = emax;
      gen_vmax[i] = vmax;
      gen_fmax[i] = fmax*(1.0 + ftol);
    }

    // These are for diagnostic output only
    //
    if (myid==0) {		// Node 0 receives
      MPI_Reduce(MPI_IN_PLACE, gen_emax.data(), gen_N, MPI_DOUBLE, MPI_SUM,
		 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, gen_vmax.data(), gen_N, MPI_DOUBLE, MPI_SUM,
		 0, MPI_COMM_WORLD);
    } else {			// Nodes >0 send
      MPI_Reduce(gen_emax.data(), 0, gen_N, MPI_DOUBLE, MPI_SUM,
		 0, MPI_COMM_WORLD);
      MPI_Reduce(gen_vmax.data(), 0, gen_N, MPI_DOUBLE, MPI_SUM,
		 0, MPI_COMM_WORLD);
    }

    // All processes need these . . .
    //
    MPI_Allreduce(MPI_IN_PLACE, gen_rloc.data(), gen_N, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, gen_mass.data(), gen_N, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, gen_fmax.data(), gen_N, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

    // Debug
    //
    if (myid==0) {

      std::ofstream test("test_multi.grid");
      if (test) {

	test << "# Rmin=" << rmin_gen
	     << "  Rmax=" << rmax_gen
	     << std::endl;

	test << std::left 
	     << std::endl << std::setfill('-') // Separator
	     << std::setw(15) << "#"	
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::endl << std::setfill(' ') // Labels
	     << std::setw(15) << "# radius"
	     << std::setw(15) << "+ mass"
	     << std::setw(15) << "+ Emax"
	     << std::setw(15) << "+ Vmax"
	     << std::setw(15) << "+ Fmax"
	     << std::setw(15) << "+ Phi(r)"
	     << std::setw(15) << "+ F_real(Phi)"
	     << std::setw(15) << "+ F_fake(Phi)"
	     << std::setw(15) << "+ Ratio"
	     << std::endl
	     << std::setw(15) << "# [1]" // Column number
	     << std::setw(15) << "+ [2]"
	     << std::setw(15) << "+ [3]"
	     << std::setw(15) << "+ [4]"
	     << std::setw(15) << "+ [5]"
	     << std::setw(15) << "+ [6]"
	     << std::setw(15) << "+ [7]"
	     << std::setw(15) << "+ [8]"
	     << std::setw(15) << "+ [9]"
	     << std::endl << std::setfill('-') // Separator
	     << std::setw(15) << "#"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::setw(15) << "+"
	     << std::endl << std::setfill(' ');

	for (int i=0; i<gen_N; i++) {
	  double r = exp(gen_rloc[i]);
	  double p = get_pot(r);
	  test << std::setw(15) << gen_rloc[i]
	       << std::setw(15) << gen_mass[i]
	       << std::setw(15) << gen_emax[i]
	       << std::setw(15) << gen_vmax[i]
	       << std::setw(15) << gen_fmax[i]
	       << std::setw(15) << p
	       << std::setw(15) << real->distf(p, 0.5)
	       << std::setw(15) << fake->distf(p, 0.5)
	       << std::setw(15) << real->get_density(r)/fake->get_density(r)
	       << std::endl;
	}
      } else {
	std::cerr << "Error opening <test_multi.grid>" << std::endl;
      }
    }

    gen_firstime = false;
  }

				// Diagnostics
  int reject=0, negmass=0;
				// Save DF evaluations
  double uuu, vvv, maxv3 = 0.0;

  // Iteration counter
  int it;			// Needs to be defined outside the loop
				// to easily check the iteration count

  // +---This keeps radius fixed.  As long as only a few states are
  // |   rejected for numerical issues, reselecting radius will be 
  // |   faster.
  // v
#if  FIXED_RADIUS>0
  // Generate a radius (outside the loop)
  mass = gen_mass[0] + Unit(random_gen)*(gen_mass[gen_N-1]-gen_mass[0]);
  r    = odd2(mass, gen_mass, gen_rloc, 0);
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  if (gen_logr) r = exp(r);
  
  pot = get_pot(r);
  vmax = sqrt(2.0*max<double>(Emax - pot, 0.0));

  // Trial velocity point
  //
  for (it=0; it<gen_itmax; it++) {
#else
  // Trial phase-space point
  //
  for (it=0; it<gen_itmax; it++) {

    // Generate a radius (inside the loop)
    mass = gen_mass[0] + Unit(random_gen)*(gen_mass[gen_N-1]-gen_mass[0]);
    r = odd2(mass, gen_mass, gen_rloc, 0);
    fmax = odd2(r, gen_rloc, gen_fmax, 1);
    if (gen_logr) r = exp(r);
  
    pot = get_pot(r);
    vmax = sqrt(2.0*max<double>(Emax - pot, 0.0));
#endif

    double xxx = 2.0*sin(asin(Unit(random_gen))/3.0);
    double yyy = (1.0 - xxx*xxx)*Unit(random_gen);

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (fmax<=0.0) continue;

    uuu = real->distf(eee, r*vt);
    vvv = fake->distf(eee, r*vt);

    if (noneg and (uuu<=0.0 or vvv<=0.0)) {
      negmass++;
      continue;
    }

    if (Unit(random_gen) > vvv/fmax ) {
      reject++;
      maxv3 = std::max<double>(maxv3, vvv);
      continue;
    }

    if (Unit(random_gen) < 0.5) vr *= -1.0;
    
    // Orientation of tangential velocity vector
    //
    double azi = 2.0*M_PI*Unit(random_gen);
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  Eigen::VectorXd out(7);

  static unsigned totcnt = 0, toomany = 0;
  totcnt++;

  if (it==gen_itmax) {
    if (verbose) {
      std::cerr << "Velocity selection failed [" << myid << "]: r="
		<< std::setw(12) << r << "["
		<< std::setw(12) << maxv3/fmax << "]"
		<< " reject="  << std::setw( 6) << reject
		<< " negmass=" << std::setw( 6) << negmass
		<< " %=" << std::setw(12)
		<< static_cast<double>(++toomany)/totcnt << std::endl;
    }
    out.setZero();
    return {out, 1};
  }

  if (Unit(random_gen)>=0.5) vr *= -1.0;

  phi  = 2.0*M_PI*Unit(random_gen);
  cost = 2.0*(Unit(random_gen) - 0.5);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  eee = std::min<double>(eee, gen_ecut);

  double dfr = real->distf(eee, r*vt);
  double dfn = fake->distf(eee, r*vt);
  double rat = dfr/dfn;

  // Deep debug
  //
#ifdef DEBUG
  if (rat <= 0.0) {
    std::cout << "[" << std::setw(3) << myid << "] Bad mass: rat=" << std::setw(16) << rat
	      << " df(M)=" << std::setw(16) << dfr
	      << " df(N)=" << std::setw(16) << dfn
	      << " dF(M)=" << std::setw(16) << uuu
	      << " dF(N)=" << std::setw(16) << vvv
	      << " r="     << std::setw(16) << r
	      << " E="     << std::setw(16) << eee
	      << std::endl;
  }
#endif
  
  out[0] = rat;
  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    
  int ierr = 0;
  for (int i=0; i<7; i++) {
    if (std::isnan(out[i]) || std::isinf(out[i])) ierr = 1;
  }

  return {out, ierr};
}


AxiSymModel::PSret SphericalModelMulti::gen_point(double radius)
{
  if (!real->dist_defined || !fake->dist_defined) {
    throw std::runtime_error("SphericalModelMulti: input distribution functions must be defined before realizing!");
  }

  double r, pot, vmax;
  double vr=0.0, vt=0.0, eee=0.0, vt1=0.0, vt2=0.0, fmax;

  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {
    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);
    double dr;

    gen_mass.resize(gen_N);
    gen_rloc.resize(gen_N);
    gen_fmax.resize(gen_N);

    if (rmin_gen <= 1.0e-16) gen_logr = 0;
    
    if (gen_logr)
      dr = (log(rmax_gen) - log(rmin_gen))/(gen_N-1);
    else
      dr = (rmax_gen - rmin_gen)/(gen_N-1);

    for (int i=0; i<gen_N; i++) {

      if (gen_logr) {
	gen_rloc[i] = log(rmin_gen) + dr*i;
	r = exp(gen_rloc[i]);
      }
      else {
	gen_rloc[i] = rmin_gen + dr*i;
	r = gen_rloc[i];
      }

      gen_mass[i] = fake->get_mass(r);

      pot = get_pot(r);
      vmax = sqrt(2.0*fabs(Emax - pot));

      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	double xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  double yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  double zzz = fake->distf(eee, r*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*(1.0 + ftol);

    }

    // Debug
    //
    ofstream test("test_multi.grid");
    if (test) {

      test << "# Rmin=" << rmin_gen
	   << "  Rmax=" << rmax_gen
	   << std::endl;

      for (int i=0; i<gen_N; i++) {
	test << std::setw(15) << gen_rloc[i]
	     << std::setw(15) << gen_mass[i]
	     << std::setw(15) << gen_fmax[i]
	     << std::setw(15) << get_pot(exp(gen_rloc[i]))
	     << std::setw(15) << fake->distf(get_pot(exp(gen_rloc[i])), 0.5)
	     << std::endl;
      }
    }

    gen_firstime = false;
  }

  r    = max<double>(rmin_gen, min<double>(rmax_gen, radius));
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  if (gen_logr) r = exp(r);
  
  pot  = get_pot(r);
  vmax = sqrt(2.0*max<double>(Emax - pot, 0.0));

  // Trial velocity point
  //
				// Diagnostic counters
  int reject=0;
  int it;			// Iteration counter
  for (it=0; it<gen_itmax; it++) {

    double xxx = 2.0*sin(asin(Unit(random_gen))/3.0);
    double yyy = (1.0 - xxx*xxx)*Unit(random_gen);

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (fmax<=0.0) continue;

    if (Unit(random_gen) > fake->distf(eee, r*vt)/fmax ) {
      reject++;
      continue;
    }

    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    // Orientation of the tangential velocity vector
    //
    double azi = 2.0*M_PI*Unit(random_gen);
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  Eigen::VectorXd out(7);

  if (it==gen_itmax) {
    if (verbose) {
      static unsigned totcnt = 0;
      std::cerr << "Velocity selection failed [" << std::setw(7) << ++totcnt
		<< "," << std::setw(4) << myid << "]: r="
		<< std::setw(12) << r
		<< std::setw(12) << fmax
		<< " reject="  << std::setw( 6) << reject
		<< std::endl;
    }
    out.setZero();
    return {out, 1};
  }

  if (Unit(random_gen)>=0.5) vr *= -1.0;

  double phi  = 2.0*M_PI*Unit(random_gen);
  double cost = 2.0*(Unit(random_gen) - 0.5);
  double sint = sqrt(1.0 - cost*cost);
  double cosp = cos(phi);
  double sinp = sin(phi);

  double dfr = real->distf(eee, r*vt);
  double dfn = fake->distf(eee, r*vt);
  double rat = dfr/dfn;

  // Deep debug
  //
#ifdef DEBUG
  if (rat <= 0.0) {
    std::cout << "[" << std::setw(3) << myid << "] Bad mass: rat="
	      << std::setw(16) << rat
	      << " df(M)=" << std::setw(16) << dfr
	      << " df(N)=" << std::setw(16) << dfn
	      << " r="     << std::setw(16) << r
	      << " E="     << std::setw(16) << eee
	      << std::endl;
  }
#endif

  out[0] = rat;
  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    
  int ierr = 0;
  for (int i=0; i<7; i++) {
    if (std::isnan(out[i]) || std::isinf(out[i])) ierr = 1;
  }
  
  return {out, ierr};
}

AxiSymModel::PSret
SphericalModelMulti::gen_point(double radius, double theta, double phi)
{
  if (!real->dist_defined || !fake->dist_defined) {
    throw std::runtime_error("SphericalModelMulti: input distribution functions must be defined before realizing!");
  }

  double r, pot, vmax, fmax;
  double vr=0.0, vt=0.0, eee=0.0, vt1=0.0, vt2=0.0;

  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {
    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);
    double dr;

    gen_mass.resize(gen_N);
    gen_rloc.resize(gen_N);
    gen_fmax.resize(gen_N);

    if (rmin_gen <= 1.0e-16) gen_logr = 0;
    
    if (gen_logr)
      dr = (log(rmax_gen) - log(rmin_gen))/(gen_N-1);
    else
      dr = (rmax_gen - rmin_gen)/(gen_N-1);

    for (int i=0; i<gen_N; i++) {

      if (gen_logr) {
	gen_rloc[i] = log(rmin_gen) + dr*i;
	r = exp(gen_rloc[i]);
      }
      else {
	gen_rloc[i] = rmin_gen + dr*i;
	r = gen_rloc[i];
      }

      gen_mass[i] = fake->get_mass(r);

      pot = get_pot(r);
      vmax = sqrt(2.0*fabs(Emax - pot));

      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	double xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  double yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  double zzz = fake->distf(eee, r*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*(1.0 + ftol);

    }

    // Debug
    //
    ofstream test("test_multi.grid");
    if (test) {

      test << "# Rmin=" << rmin_gen
	   << "  Rmax=" << rmax_gen
	   << std::endl;

      for (int i=0; i<gen_N; i++) {
	test << std::setw(15) << gen_rloc[i]
	     << std::setw(15) << gen_mass[i]
	     << std::setw(15) << gen_fmax[i]
	     << std::setw(15) << get_pot(exp(gen_rloc[i]))
	     << std::setw(15) << fake->distf(get_pot(exp(gen_rloc[i])), 0.5)
	     << std::endl;
      }
    }

    gen_firstime = false;
  }

  r    = max<double>(rmin_gen, min<double>(rmax_gen, radius));
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  if (gen_logr) r = exp(r);
  
  pot  = get_pot(r);
  vmax = sqrt(2.0*max<double>(Emax - pot, 0.0));

  // Trial velocity point
  //
				// Diagnostic counters
  int reject=0;
  int it;			// Iteration counter
  for (it=0; it<gen_itmax; it++) {

    double xxx = 2.0*sin(asin(Unit(random_gen))/3.0);
    double yyy = (1.0 - xxx*xxx)*Unit(random_gen);

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (fmax<=0.0) continue;

    if (Unit(random_gen) > fake->distf(eee, r*vt)/fmax ) {
      reject++;
      continue;
    }

    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    // Orientation of tangential velocity vector
    //
    double azi = 2.0*M_PI*Unit(random_gen);
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  Eigen::VectorXd out(7);

  if (it==gen_itmax) {
    if (verbose) {
      static unsigned totcnt = 0;
      std::cerr << "Velocity selection failed [" << std::setw(7) << ++totcnt
		<< "," << std::setw(4) << myid << "]: r="
		<< std::setw(12) << r
		<< std::setw(12) << fmax
		<< " reject="  << std::setw( 6) << reject
		<< std::endl;
    }
    out.setZero();
    return {out, 1};
  }

  if (Unit(random_gen)>=0.5) vr *= -1.0;

  double dfr = real->distf(eee, r*vt);
  double dfn = fake->distf(eee, r*vt);
  double rat = dfr/dfn;

  // Deep debug
  //
#ifdef DEBUG
  if (rat <= 0.0) {
    std::cout << "[" << std::setw(3) << myid << "] Bad mass: rat="
	      << std::setw(16) << rat
	      << " df(M)=" << std::setw(16) << dfr
	      << " df(N)=" << std::setw(16) << dfn
	      << " r="     << std::setw(16) << r
	      << " E="     << std::setw(16) << eee
	      << std::endl;
  }
#endif

  double cost = cos(theta);
  double sint = sin(theta);
  double cosp = cos(phi);
  double sinp = sin(phi);

  out[0] = rat;
  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    
  int ierr = 0;
  for (int i=0; i<7; i++) {
    if (std::isnan(out[i]) || std::isinf(out[i])) ierr = 1;
  }
  
  return {out, ierr};
}

AxiSymModel::PSret SphericalModelMulti::gen_point
(double Emin, double Emax, double Kmin,  double Kmax)
{
  if (!real->dist_defined || !fake->dist_defined) {
    throw std::runtime_error("SphericalModelMulti: input distribution functions must be defined before realizing!");
  }
  
  double r, vr, vt, vt1, vt2, E, K, J, jmax, w1t, pot;
  double phi, sint, cost, sinp, cosp;
  double rmin = max<double>(get_min_radius(), gen_rmin);
  double rmax = get_max_radius();


  if (gen_firstime_E) {
    
    gen_orb = SphericalOrbit(shared_from_this());
    
    Emin_grid = get_pot(rmin)*(1.0 - gen_tolE);
    Emax_grid = get_pot(rmax)*(1.0 + gen_tolE);

    if (verbose) {
      std::cout << "Initial radial bounds: " << rmin
	   << ", " << rmax << std::endl;

      std::cout << "Initial energy bounds: " << Emin_grid 
	   << ", " << Emax_grid << " [" << gen_tolE << "]" << std::endl;
    }

				// Enforce limits
    Emin_grid = max<double>(Emin, Emin_grid);
    Emin_grid = min<double>(Emax, Emin_grid);
    Emax_grid = max<double>(Emin, Emax_grid);
    Emax_grid = min<double>(Emax, Emax_grid);

    dEgrid = (Emax_grid - Emin_grid)/(gen_E-1);
    dKgrid = (1.0 - 2.0*gen_tolK)/(gen_K-1);

    for (int i=0; i<gen_E; i++) Egrid.push_back(Emin_grid + dEgrid*i);
    for (int i=0; i<gen_K; i++) Kgrid.push_back(gen_tolK  + dKgrid*i);

    double ans, K, angles = 2.0*pow(2.0*M_PI, 3.0), tfac;
    vector<double> tmass;
    ANGLE_GRID *agrid;

    for (int i=0; i<gen_E; i++) {
      ans = 0.0;

      vector<WRgrid> wrvec;
      for (int j=0; j<gen_K; j++) {

	// Trapezoidal rule factor
	//
	if (j==0 || j==gen_K) tfac = 0.5;
	else tfac = 1.0;

	K = gen_tolK + dKgrid*j;
	gen_orb.new_orbit(Egrid[i], K);
	ans += K*gen_orb.Jmax()*gen_orb.Jmax() /
	  gen_orb.get_freq(1) * fake->distf(Egrid[i], gen_orb.AngMom()) * tfac;

	WRgrid wr;
	agrid = gen_orb.get_angle_grid();
	
	for (int n=0; n<agrid->num; n++) {
	  wr.w1.push_back(agrid->w1(0, n));
	  wr.r.push_back(agrid->r(0, n));
	}
	wrvec.push_back(wr);

      } // End Kappa loop

      Jmax.push_back(gen_orb.Jmax());
      tmass.push_back(ans*angles*dKgrid*dEgrid);
      Rgrid.push_back(wrvec);

    } // Energy E loop

    EgridMass.push_back(0.0);
    for (int i=1; i<gen_E; i++) 
      EgridMass.push_back(EgridMass[i-1]+tmass[i]);

    // Debug
    
    ofstream test("test.grid");
    if (test) {

      test << "# Emin=" << Emin_grid
	   << "  Emax=" << Emax_grid
	   << std::endl;

      for (int i=0; i<gen_E; i++) {
	test << std::setw(15) << Egrid[i]
	     << std::setw(15) << EgridMass[i]
	     << std::setw(15) << Jmax[i]
	     << std::endl;
      }
    }
    
    gen_firstime_E = false;
  }

  // Enforce limits
  //
  Emin = max<double>(Emin, Emin_grid);
  Emin = min<double>(Emin, Emax_grid);
  Emax = max<double>(Emax, Emin_grid);
  Emax = min<double>(Emax, Emax_grid);

  double Mmin = odd2(Emin, Egrid, EgridMass, 1);
  double Mmax = odd2(Emax, Egrid, EgridMass, 1);
  double kmin = max<double>(Kmin, gen_tolK);
  double kmax = min<double>(Kmax, 1.0 - gen_tolK);

  E = odd2(Mmin + (Mmax-Mmin)*Unit(random_gen), EgridMass, Egrid, 0);
  K = sqrt(kmin*kmin + (kmax*kmax - kmin*kmin)*Unit(random_gen));

  int indxE = int( (E - Emin_grid)/dEgrid );
  int indxK = int( (K - gen_tolK)/dKgrid );
  
  indxE = max<int>(0, min<int>(gen_E-2, indxE));
  indxK = max<int>(0, min<int>(gen_K-2, indxK));

  double cE[2], cK[2];

  cE[1] = (E - (Emin_grid + dEgrid*indxE))/dEgrid;
  cE[0] = 1.0 - cE[1];

  cK[1] = (K - (gen_tolK + dKgrid*indxK))/dKgrid;
  cK[0] = 1.0 - cK[1];


  r = 0.0;
  J = 0.0;
  jmax = 0.0;
  w1t = M_PI*Unit(random_gen);

  for (int ie=0; ie<2; ie++) {
    J += cE[ie]*Jmax[indxE+ie] * K;
    jmax += cE[ie]*Jmax[indxE+ie];
    for (int ik=0; ik<2; ik++) {
      r += cE[ie]*cK[ik]*odd2(w1t, Rgrid[indxE+ie][indxK+ik].w1, 
			     Rgrid[indxE+ie][indxK+ik].r, 0);
    }
  }
      
  pot = get_pot(r);
  vt = J/r;

  // Interpolation check (should be rare). Error condition is set.
  //
  int ierr = 0;
  if (2.0*(E - pot) - vt*vt < 0.0) {
    ierr = 1;
    if (E < pot) E = pot;
    vt = sqrt(E - pot);
  }
  vr = sqrt( 2.0*(E - pot) - vt*vt );

  if (Unit(random_gen)<0.5) vr *= -1.0;
    
  // Orientation of the tangential velocity vector
  //
  double azi = 2.0*M_PI*Unit(random_gen);
  vt1 = vt*cos(azi);
  vt2 = vt*sin(azi);

  Eigen::VectorXd out(7);

  phi = 2.0*M_PI*Unit(random_gen);
  cost = 2.0*(Unit(random_gen) - 0.5);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  out[0] = real->distf(E, r*vt)/fake->distf(E, r*vt);
  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    
  for (int i=0; i<7; i++) {
    if (std::isnan(out[i]) || std::isinf(out[i])) ierr = 1;
  }

  return {out, ierr};
}


