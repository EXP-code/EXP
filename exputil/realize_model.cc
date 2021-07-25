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
#include <cmath>

#include <boost/random/mersenne_twister.hpp>

#include <massmodel.H>
#include <interp.H>

#ifdef DEBUG
#include <orbit.H>
static SphericalOrbit orb;
#endif


bool     AxiSymModel::gen_EJ = true;
int      AxiSymModel::numr = 50;
int      AxiSymModel::numj = 50;
int      AxiSymModel::gen_N = 400;
int      AxiSymModel::gen_E = 400;
int      AxiSymModel::gen_K = 200;
double   AxiSymModel::gen_tolE = 0.01;
double   AxiSymModel::gen_tolK = 0.02;
double   AxiSymModel::gen_rmin = 0.0;
int      AxiSymModel::gen_logr = 1;
double   AxiSymModel::gen_kmin = 0.0;
uint32_t AxiSymModel::gen_seed = 11;
int      AxiSymModel::gen_itmax = 20000;

const bool verbose = true;
const double ftol = 0.01;

extern boost::mt19937 random_gen;

Eigen::VectorXd AxiSymModel::gen_point_2d(int& ierr)
{
  if (!dist_defined) {
    cerr << "AxiSymModel: must define distribution before realizing!\n";
    _exit (-1);
  }

  double r=0.0, pot, vmax, xxx, yyy, zzz, vr=0.0, vt=0.0, eee, fmax, vv;
  double phi=0.0, sinp, cosp;
  double T, w1;
  double tol = 1.0e-5;
  double rmin = max<double>(get_min_radius(), gen_rmin);

  double Emin = get_pot(rmin);
  double Emax = get_pot(get_max_radius());

  int it = 0;

  if (gen_EJ) {

    if (gen_firstime) {

      gen_orb = SphericalOrbit(this);

      double dx = (Emax - Emin - 2.0*tol)/(numr-1);
      double dy = (1.0 - gen_kmin - 2.0*tol)/(numj-1);

      cout << "gen_point_2d[" << ModelID << "]: " << get_max_radius() << endl;
      
      gen_fomax = 0.0;

      for (int j=0; j<numr; j++) {
	xxx = Emin + tol + dx*j;

	for (int k=0; k<numj; k++) {
	  yyy = gen_kmin + tol + dy*k;

	  gen_orb.new_orbit(xxx, yyy);
	    
	  zzz = distf(xxx, gen_orb.AngMom())*gen_orb.Jmax()
	    /gen_orb.get_freq(1);
	  gen_fomax = zzz>gen_fomax ? zzz : gen_fomax;
	}
      }

      gen_firstime = false;
    }

    
				// Trial velocity point
    
    for (it=0; it<gen_itmax; it++) {

      xxx = Emin + tol + (Emax - Emin - 2.0*tol)*Unit(random_gen);
      yyy = gen_kmin + tol + (1.0 - gen_kmin - 2.0*tol)*Unit(random_gen);

      gen_orb.new_orbit(xxx, yyy);

      zzz = distf(xxx, gen_orb.AngMom()) * gen_orb.Jmax()/gen_orb.get_freq(1);

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

      double tol = 1.0e-5;
      double dx = (1.0 - 2.0*tol)/(numr-1);
      double dy = (1.0 - 2.0*tol)/(numj-1);
      double dr;

      gen_mass.resize(gen_N);
      gen_rloc.resize(gen_N);
      gen_fmax.resize(gen_N);

      cout << "gen_point_2d[" << ModelID << "]: " << get_max_radius() << endl;

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
	  xxx = sqrt(tol + dx*j);

	  for (int k=0; k<numj; k++) {
	    yyy = 0.5*M_PI*(tol + dy*k);

	    vr = vmax*xxx*cos(yyy);
	    vt = vmax*xxx*sin(yyy);
	    eee = pot + 0.5*(vr*vr + vt*vt);

	    zzz = distf(eee, gen_rloc[i]*vt);
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
    
    for (it=0; it<gen_itmax; it++) {

      xxx = sqrt(Unit(random_gen));
      yyy = 0.5*M_PI*Unit(random_gen);

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

  if (it==gen_itmax) {
    cerr << "Velocity selection failed, r=" << r << "\n";
    ierr = 1;
    out.setZero();
    return out;
  }

  ierr = 0;
  
  cosp = cos(phi);
  sinp = sin(phi);

  out[0] = 1.0;
  out[1] = r * cosp;
  out[2] = r * sinp;
  out[3] = 0.0;
  out[4] = vr * cosp - vt * sinp;
  out[5] = vr * sinp + vt * cosp;
  out[6] = 0.0;
    
  return out;
}


Eigen::VectorXd AxiSymModel::gen_point_2d(double r, int& ierr)
{
  if (!dist_defined) {
    cerr << "AxiSymModel: must define distribution before realizing!\n";
    _exit (-1);
  }

  int it;
  double pot, vmax, xxx, yyy, zzz, vr=0.0, vt=0.0, eee, fmax;
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

    cout << "gen_point_2d[" << ModelID << "]: " << rmin
	 << ", " << get_max_radius() << endl;

    double dr = (get_max_radius() - rmin)/gen_N;

    for (int i=0; i<gen_N; i++) {
      gen_rloc[i] = rmin + dr*i;
      gen_mass[i] = get_mass(gen_rloc[i]);

      pot = get_pot(gen_rloc[i]);
      vmax = sqrt(2.0*fabs(Emax - pot));

      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	xxx = sqrt(tol + dx*j);

	for (int k=0; k<numj; k++) {
	  yyy = 0.5*M_PI*(tol + dy*k);

	  vr = vmax*xxx*cos(yyy);
	  vt = vmax*xxx*sin(yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  zzz = distf(eee, gen_rloc[i]*vt);
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
    
  for (it=0; it<gen_itmax; it++) {

    xxx = sqrt(Unit(random_gen));
    yyy = 0.5*M_PI*Unit(random_gen);

    vr = vmax*xxx*cos(yyy);
    vt = vmax*xxx*sin(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (Unit(random_gen) > distf(eee, r*vt)/fmax ) continue;
    
    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    phi = 2.0*M_PI*Unit(random_gen);

    break;
  }

  
  Eigen::VectorXd out(7);

  if (it==gen_itmax) {
    cerr << "Velocity selection failed, r=" << r << "\n";
    ierr = 1;
    out.setZero();
    return out;
  }

  ierr = 0;
  
  cosp = cos(phi);
  sinp = sin(phi);

  out[0] = 1.0;
  out[1] = r * cosp;
  out[2] = r * sinp;
  out[3] = 0.0;
  out[4] = vr * cosp - vt * sinp;
  out[5] = vr * sinp + vt * cosp;
  out[6] = 0.0;
    
  return out;
}


Eigen::VectorXd AxiSymModel::gen_point_3d(int& ierr)
{
  if (!dist_defined) {
    cerr << "AxiSymModel: must define distribution before realizing!\n";
    _exit (-1);
  }

#ifdef DEBUG
  static ofstream tout("gen3d.ktest");
#endif

  double r, pot, vmax, xxx, yyy, vr=0.0, vt, eee, vt1=0.0, vt2=0.0, fmax;
  double phi, sint, cost, sinp, cosp, azi;

  double rmin = max<double>(get_min_radius(), gen_rmin);
  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {
    double zzz;
    // const int numr = 200;
    // const int numj = 200;

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
	xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  zzz = distf(eee, r*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*(1.0 + ftol);

    }

    // Debug
    
    std::ofstream test("test.grid");
    if (test) {

      test << "# Rmin=" << rmin
	   << "  Rmax=" << get_max_radius()
	   << endl;

      for (int i=0; i<gen_N; i++) {
	test << setw(15) << gen_rloc[i]
	     << setw(15) << gen_mass[i]
	     << setw(15) << gen_fmax[i]
	     << endl;
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
    
  int it;
  double angmom[3];

  for (it=0; it<gen_itmax; it++) {

    xxx = -2.0*cos(acos(Unit(random_gen))/3.0 - 2.0*M_PI/3.0);
    yyy = (1.0 - xxx*xxx)*Unit(random_gen);

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    // Debug
    /*
    if (sqrt(vr*vr + vt*vt)>0.99*vmax) {
      cout << "Check: df val = " << distf(eee, r*vt)/fmax 
	   << "  v/vmax = " << sqrt(vr*vr + vt*vt)/vmax 
	   << "  eee = " << eee
	   << endl;
    }
    */

    if (Unit(random_gen) > distf(eee, r*vt)/fmax ) continue;

    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    azi = 2.0*M_PI*Unit(random_gen);
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  Eigen::VectorXd out(7);

  if (it==gen_itmax) {
    std::cerr << "Velocity selection failed, r=" << r << std::endl;
    ierr = 1;
    out.setZero();
    return out;
  }

  ierr = 0;
  
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
  tout << setw(15) << eee
       << setw(15) << sqrt(angmom[0]*angmom[0]+
			   angmom[1]*angmom[1]+
			   angmom[2]*angmom[2])/orb.Jmax()
       << setw(15) << r
       << setw(15) << sqrt(out[4]*out[4]+out[5]*out[5]+out[6]*out[6])
       << "\n";
#endif


  return out;
}


Eigen::VectorXd AxiSymModel::gen_point_3d(double Emin, double Emax, 
					  double Kmin, double Kmax, int& ierr)
{
  if (!dist_defined) {
    cerr << "AxiSymModel: must define distribution before realizing!\n";
    _exit (-1);
  }

#ifdef DEBUG
  double angmom[3];
  static ofstream tout("gen3d.Etest");
#endif

  double r, vr, vt, vt1, vt2, E, K, J, jmax, w1t, eee, pot;
  double phi, sint, cost, sinp, cosp, azi;
  double rmin = max<double>(get_min_radius(), gen_rmin);
  double rmax = get_max_radius();

  if (gen_firstime_E) {

#ifdef DEBUG
    orb = SphericalOrbit(this);
#endif
    gen_orb = SphericalOrbit(this);

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
	  wr.w1.push_back(agrid->w1(1, n));
	  wr.r.push_back(agrid->r(1, n));
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
	   << endl;

      for (int i=0; i<gen_E; i++) {
	test << setw(15) << Egrid[i]
	     << setw(15) << EgridMass[i]
	     << setw(15) << Jmax[i]
	     << endl;
      }
    }

    gen_firstime_E = false;

#ifdef DEBUG
    tout << "# Emin=" << Emin_grid << "  Emax=" << Emax_grid
	 << "  Mass frac=" << 
      ( odd2(Emax, Egrid, EgridMass, 1) - 
	odd2(Emin, Egrid, EgridMass, 1) )/EgridMass[gen_E-1] << endl;
#endif    
  }

				// Enforce limits
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
  ierr = 0;
				// Interpolation check (should be rare)
				// Error condition is set
  if (2.0*(E - pot) - vt*vt < 0.0) {
    ierr = 1;
    if (E < pot) E = pot;
    vt = sqrt(E - pot);
  }
  vr = sqrt( 2.0*(E - pot) - vt*vt );

  if (Unit(random_gen)<0.5) vr *= -1.0;
    
  azi = 2.0*M_PI*Unit(random_gen);
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
  if (ierr) return out;

  eee = pot + 0.5*(out[4]*out[4]+out[5]*out[5]+out[6]*out[6]);
  orb.new_orbit(eee, 0.5);
  angmom[0] = out[2]*out[6] - out[3]*out[5];
  angmom[1] = out[3]*out[4] - out[1]*out[6];
  angmom[2] = out[1]*out[5] - out[2]*out[4];
  tout << setw(15) << eee
       << setw(15) << E
       << setw(15) << sqrt(angmom[0]*angmom[0]+
			   angmom[1]*angmom[1]+
			   angmom[2]*angmom[2])/orb.Jmax()
       << setw(15) << K
       << setw(15) << r
       << setw(15) << sqrt(out[4]*out[4]+out[5]*out[5]+out[6]*out[6])
       << setw(15) << orb.Jmax()
       << setw(15) << jmax
       << "\n";
#endif

  return out;
}



Eigen::VectorXd AxiSymModel::gen_point_jeans_3d(int& ierr)
{
  double r, d, xxx, yyy, vr, vt, vt1, vt2, vv, vtot;
  double phi, sint, cost, sinp, cosp, azi;
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
    
    ofstream test("test.grid");
    if (test) {

      test << "# [Jeans] Rmin=" << rmin
	   << "  Rmax=" << get_max_radius()
	   << endl;

      for (int i=0; i<gen_N; i++) {
	if (gen_logr) r = exp(gen_rloc[i]);
	else r = gen_rloc[i];

	d = get_density(r);
	if (d>0.0 && gen_fmax[i]>=0.0)
	  vtot = sqrt(gen_fmax[i]/d);
	else
	  vtot = 0.0;

	test << setw(15) << gen_rloc[i]
	     << setw(15) << gen_mass[i]
	     << setw(15) << gen_fmax[i]
	     << setw(15) << vtot
	     << endl;
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
    
  
  xxx = -2.0*cos(acos(Unit(random_gen))/3.0 - 2.0*M_PI/3.0);
  yyy = (1.0 - xxx*xxx)*Unit(random_gen);

  vr = vtot*xxx;
  vt = vtot*sqrt(yyy);

  azi = 2.0*M_PI*Unit(random_gen);
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
    
  ierr = 0;

  return out;
}

void AxiSymModel::gen_velocity(double* pos, double* vel, int& ierr)
{
  if (dof()!=3)
      bomb( "AxiSymModel: gen_velocity only implemented for 3d model!" );

  if (!dist_defined) {
    cerr << "AxiSymModel: must define distribution before realizing!\n";
    _exit (-1);
  }

  double r, pot, vmax, xxx, yyy, vr=0.0, vt, eee, vt1=0.0, vt2=0.0, fmax;
  double phi, sint, cost, sinp, cosp, azi;
  double rmin = max<double>(get_min_radius(), gen_rmin);

  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {
    double zzz;
    // const int numr = 200;
    // const int numj = 200;

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
	xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  zzz = distf(eee, r*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*(1.0 + ftol);

    }
    gen_firstime = false;

#ifdef DEBUG
    ofstream tout("test.dfgrid");
    tout << "# Rmin=" << rmin 
	 << "  Rmax=" << get_max_radius() << endl;
    for (int i=0; i<gen_N; i++)
      tout << setw(18) << gen_rloc[i]
	   << setw(18) << gen_mass[i]
	   << setw(18) << gen_fmax[i]
	   << endl;
#endif
  }

  r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
  if (gen_logr) r = log(r);
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  if (gen_logr) r = exp(r);
  
  pot = get_pot(r);
  vmax = sqrt(2.0*fabs(Emax - pot));

                                // Trial velocity point
    
  int it;

  for (it=0; it<gen_itmax; it++) {

    xxx = -2.0*cos(acos(Unit(random_gen))/3.0 - 2.0*M_PI/3.0);
    yyy = (1.0 - xxx*xxx)*Unit(random_gen);

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (Unit(random_gen) > distf(eee, r*vt)/fmax ) continue;

    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    azi = 2.0*M_PI*Unit(random_gen);
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  if (it==gen_itmax) {
    cerr << "Velocity selection failed, r=" << r << "\n";
    ierr = 1;
    return;
  }

  ierr = 0;
  
  if (Unit(random_gen)>=0.5) vr *= -1.0;

  phi = atan2(pos[1], pos[0]);
  cost = pos[2]/(r+1.0e-18);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  vel[0] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  vel[1] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  vel[2] = vr * cost      - vt1 * sint;
}

Eigen::VectorXd SphericalModelMulti::gen_point(int& ierr)
{
  if (!real->dist_defined || !fake->dist_defined) {
    cerr << "SphericalModelMulti: input distribution functions must be defined before realizing!\n";
    exit (-1);
  }

  double r, pot, vmax, xxx, yyy;
  double vr=0.0, vt=0.0, eee=0.0, vt1=0.0, vt2=0.0, fmax, emax;
  double mass, phi, sint, cost, sinp, cosp, azi;

  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {
    double zzz;
    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);
    double dr;

    gen_mass.resize(gen_N);
    gen_rloc.resize(gen_N);
    gen_fmax.resize(gen_N);
    vector<double> gen_emax, gen_vmax;

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

      emax = pot;
      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  zzz = fake->distf(eee, r*vt);
	  if (zzz>fmax) {
	    emax = eee;
	  }
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_emax.push_back(emax);
      gen_vmax.push_back(vmax);
      gen_fmax[i] = fmax*(1.0 + ftol);

    }

    // Debug
    
    ofstream test("test_multi.grid");
    if (test) {

      test << "# Rmin=" << rmin_gen
	   << "  Rmax=" << rmax_gen
	   << endl;

      test << left 
	   << endl << setfill('-') // Separator
	   << setw(15) << "#"	
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << endl << setfill(' ') // Labels
	   << setw(15) << "# radius"
	   << setw(15) << "+ mass"
	   << setw(15) << "+ Emax"
	   << setw(15) << "+ Vmax"
	   << setw(15) << "+ Fmax"
	   << setw(15) << "+ Phi(r)"
	   << setw(15) << "+ F_real(Phi)"
	   << setw(15) << "+ F_fake(Phi)"
	   << setw(15) << "+ Ratio"
	   << endl
	   << setw(15) << "# [1]" // Column number
	   << setw(15) << "+ [2]"
	   << setw(15) << "+ [3]"
	   << setw(15) << "+ [4]"
	   << setw(15) << "+ [5]"
	   << setw(15) << "+ [6]"
	   << setw(15) << "+ [7]"
	   << setw(15) << "+ [8]"
	   << setw(15) << "+ [9]"
	   << endl << setfill('-') // Separator
	   << setw(15) << "#"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << setw(15) << "+"
	   << endl << setfill(' ');

      for (int i=0; i<gen_N; i++) {
	double r = exp(gen_rloc[i]);
	double p = get_pot(r);
	test << setw(15) << gen_rloc[i]
	     << setw(15) << gen_mass[i]
	     << setw(15) << gen_emax[i]
	     << setw(15) << gen_vmax[i]
	     << setw(15) << gen_fmax[i]
	     << setw(15) << p
	     << setw(15) << real->distf(p, 0.5)
	     << setw(15) << fake->distf(p, 0.5)
	     << setw(15) << real->get_density(r)/fake->get_density(r)
	     << endl;
      }
    }

    gen_firstime = false;
  }

#if 1
  mass = gen_mass[0] + Unit(random_gen)*(gen_mass[gen_N-1]-gen_mass[0]);
  r = odd2(mass, gen_mass, gen_rloc, 0);
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  if (gen_logr) r = exp(r);
  
  pot = get_pot(r);
  vmax = sqrt(2.0*max<double>(Emax - pot, 0.0));

                                // Trial velocity point
    
  int it;
  for (it=0; it<gen_itmax; it++) {
#else
  int it;
  for (it=0; it<gen_itmax; it++) {

  mass = gen_mass[0] + Unit(random_gen)*(gen_mass[gen_N-1]-gen_mass[0]);
  r = odd2(mass, gen_mass, gen_rloc, 0);
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  if (gen_logr) r = exp(r);
  
  pot = get_pot(r);
  vmax = sqrt(2.0*max<double>(Emax - pot, 0.0));

#endif

    xxx = 2.0*sin(asin(Unit(random_gen))/3.0);
    yyy = (1.0 - xxx*xxx)*Unit(random_gen);

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (fmax<=0.0) continue;
    if (Unit(random_gen) > fake->distf(eee, r*vt)/fmax ) continue;

    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    azi = 2.0*M_PI*Unit(random_gen);
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  Eigen::VectorXd out(7);

  if (it==gen_itmax) {
    if (verbose) cerr << "Velocity selection failed, r=" << r << "\n";
    ierr = 1;
    out.setZero();
    return out;
  }

  ierr = 0;
  
  if (Unit(random_gen)>=0.5) vr *= -1.0;

  phi = 2.0*M_PI*Unit(random_gen);
  cost = 2.0*(Unit(random_gen) - 0.5);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  eee = std::min<double>(eee, gen_ecut);

  out[0] = real->distf(eee, r*vt)/fake->distf(eee, r*vt);
  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    
  for (int i=0; i<7; i++) {
    if (std::isnan(out[i]) || std::isinf(out[i])) ierr = 1;
  }

  /*
  std::cout << "Mass[0]: "
	    << std::setw(18) << eee
	    << std::setw(18) << real->distf(eee, r*vt)
	    << std::setw(18) << fake->distf(eee, r*vt) << std::endl;
  */

  return out;
}


Eigen::VectorXd SphericalModelMulti::gen_point(double radius, int& ierr)
{
  if (!real->dist_defined || !fake->dist_defined) {
    cerr << "SphericalModelMulti: input distribution functions must be defined before realizing!\n";
    exit (-1);
  }

  double r, pot, vmax, xxx, yyy;
  double vr=0.0, vt=0.0, eee=0.0, vt1=0.0, vt2=0.0, fmax;
  double phi, sint, cost, sinp, cosp, azi;

  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {
    double zzz;
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
	xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  zzz = fake->distf(eee, r*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*(1.0 + ftol);

    }

    // Debug
    
    ofstream test("test_multi.grid");
    if (test) {

      test << "# Rmin=" << rmin_gen
	   << "  Rmax=" << rmax_gen
	   << endl;

      for (int i=0; i<gen_N; i++) {
	test << setw(15) << gen_rloc[i]
	     << setw(15) << gen_mass[i]
	     << setw(15) << gen_fmax[i]
	     << setw(15) << get_pot(exp(gen_rloc[i]))
	     << setw(15) << fake->distf(get_pot(exp(gen_rloc[i])), 0.5)
	     << endl;
      }
    }

    gen_firstime = false;
  }

  r = max<double>(rmin_gen, min<double>(rmax_gen, radius));
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  if (gen_logr) r = exp(r);
  
  pot = get_pot(r);
  vmax = sqrt(2.0*max<double>(Emax - pot, 0.0));

                                // Trial velocity point
    
  int it;
  for (it=0; it<gen_itmax; it++) {

    xxx = 2.0*sin(asin(Unit(random_gen))/3.0);
    yyy = (1.0 - xxx*xxx)*Unit(random_gen);

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if (fmax<=0.0) continue;
    if (Unit(random_gen) > fake->distf(eee, r*vt)/fmax ) continue;

    if (Unit(random_gen)<0.5) vr *= -1.0;
    
    azi = 2.0*M_PI*Unit(random_gen);
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  Eigen::VectorXd out(7);

  if (it==gen_itmax) {
    if (verbose) cerr << "Velocity selection failed, r=" << r << "\n";
    ierr = 1;
    out.setZero();
    return out;
  }

  ierr = 0;
  
  if (Unit(random_gen)>=0.5) vr *= -1.0;

  phi = 2.0*M_PI*Unit(random_gen);
  cost = 2.0*(Unit(random_gen) - 0.5);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  out[0] = real->distf(eee, r*vt)/fake->distf(eee, r*vt);
  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    
  for (int i=0; i<7; i++) {
    if (std::isnan(out[i]) || std::isinf(out[i])) ierr = 1;
  }

  /*
  std::cout << "Mass[1]: "
	    << std::setw(18) << eee
	    << std::setw(18) << real->distf(eee, r*vt)
	    << std::setw(18) << fake->distf(eee, r*vt) << std::endl;
  */

  return out;
}


Eigen::VectorXd
  SphericalModelMulti::gen_point(double Emin, double Emax, double Kmin,
				 double Kmax, int& ierr)
{
  if (!real->dist_defined || !fake->dist_defined) {
    cerr << "SphericalModelMulti: input distribution functions must be defined before realizing!\n";
    exit (-1);
  }
  
  double r, vr, vt, vt1, vt2, E, K, J, jmax, w1t, pot;
  double phi, sint, cost, sinp, cosp, azi;
  double rmin = max<double>(get_min_radius(), gen_rmin);
  double rmax = get_max_radius();


  if (gen_firstime_E) {
    
    gen_orb = SphericalOrbit(this);
    
    Emin_grid = get_pot(rmin)*(1.0 - gen_tolE);
    Emax_grid = get_pot(rmax)*(1.0 + gen_tolE);

    if (verbose) {
      cout << "Initial radial bounds: " << rmin
	   << ", " << rmax << endl;

      cout << "Initial energy bounds: " << Emin_grid 
	   << ", " << Emax_grid << " [" << gen_tolE << "]" << endl;
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
	if (j==0 || j==gen_K) tfac = 0.5;
	else tfac = 1.0;

	K = gen_tolK + dKgrid*j;
	gen_orb.new_orbit(Egrid[i], K);
	ans += K*gen_orb.Jmax()*gen_orb.Jmax() /
	  gen_orb.get_freq(1) * fake->distf(Egrid[i], gen_orb.AngMom()) * tfac;

	WRgrid wr;
	agrid = gen_orb.get_angle_grid();
	
	for (int n=0; n<agrid->num; n++) {
	  wr.w1.push_back(agrid->w1(1, n));
	  wr.r.push_back(agrid->r(1, n));
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
	   << endl;

      for (int i=0; i<gen_E; i++) {
	test << setw(15) << Egrid[i]
	     << setw(15) << EgridMass[i]
	     << setw(15) << Jmax[i]
	     << endl;
      }
    }
    
    gen_firstime_E = false;
  }

				// Enforce limits
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
  ierr = 0;
				// Interpolation check (should be rare)
				// Error condition is set
  if (2.0*(E - pot) - vt*vt < 0.0) {
    ierr = 1;
    if (E < pot) E = pot;
    vt = sqrt(E - pot);
  }
  vr = sqrt( 2.0*(E - pot) - vt*vt );

  if (Unit(random_gen)<0.5) vr *= -1.0;
    
  azi = 2.0*M_PI*Unit(random_gen);
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

  return out;
}


