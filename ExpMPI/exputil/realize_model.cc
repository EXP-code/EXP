// This may look like C code, but it is really -*- C++ -*-

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

#include <unistd.h>
#include <math.h>

#include <iostream>
#include <string>

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>

#include <Vector.h>
#include <interp.h>
#include <massmodel.h>


bool AxiSymModel::gen_EJ = true;
int AxiSymModel::numr = 50;
int AxiSymModel::numj = 50;
int AxiSymModel::gen_N = 400;
double AxiSymModel::gen_kmin = 0.0;
unsigned int AxiSymModel::gen_seed = 11;
int AxiSymModel::gen_itmax = 4000;


Vector AxiSymModel::gen_point_2d(int& ierr)
{
  if (!dist_defined) {
    cerr << "AxiSymModel: must define distribution before realizing!\n";
    _exit (-1);
  }

  int it;
  double r, pot, vmax, xxx, yyy, zzz, vr, vt, eee, fmax, vv;
  double phi, sinp, cosp;
  double T, w1;
  double tol = 1.0e-5;

  double Emin = get_pot(get_min_radius());
  double Emax = get_pot(get_max_radius());


  if (gen_EJ) {

    if (gen_firstime) {

      gen_orb = SphericalOrbit(this);

      double dx = (Emax - Emin - 2.0*tol)/(numr-1);
      double dy = (1.0 - gen_kmin - 2.0*tol)/(numj-1);

      gen = new ACG(gen_seed, 20);
      Unit = new Uniform(0.0, 1.0, gen);

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

      xxx = Emin + tol + (Emax - Emin - 2.0*tol)*(*Unit)();
      yyy = gen_kmin + tol + (1.0 - gen_kmin - 2.0*tol)*(*Unit)();

      gen_orb.new_orbit(xxx, yyy);

      zzz = distf(xxx, gen_orb.AngMom()) * gen_orb.Jmax()/gen_orb.get_freq(1);

      if ((*Unit)() > zzz/gen_fomax ) continue;

      w1 = 2.0*M_PI*(*Unit)();
      T = w1/gen_orb.get_freq(1);
      
      r = gen_orb.get_angle(6, T);
      phi = 2.0*M_PI*(*Unit)() - gen_orb.get_angle(5, T);

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

      gen = new ACG(gen_seed, 20);
      Unit = new Uniform(0.0, 1.0, gen);
      
      gen_mass.setsize(1, gen_N);
      gen_rloc.setsize(1, gen_N);
      gen_fmax.setsize(1, gen_N);

      cout << "gen_point_2d[" << ModelID << "]: " << get_max_radius() << endl;

      double dr = (get_max_radius() - get_min_radius())/gen_N;

      for (int i=1; i<=gen_N; i++) {
	gen_rloc[i] = get_min_radius() + dr*((double)i-0.5);
	gen_mass[i] = get_mass(gen_rloc[i]);

	pot = get_pot(gen_rloc[i]);
	vmax = sqrt(2.0*(Emax - pot));

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
	gen_fmax[i] = fmax*1.1;

      }
      gen_firstime = false;
    }

    r = odd2((*Unit)()*gen_mass[gen_N], gen_mass, gen_rloc, 0);
    fmax = odd2(r, gen_rloc, gen_fmax, 1);
    pot = get_pot(r);
    vmax = sqrt(2.0*(Emax - pot));

                                // Trial velocity point
    
    for (it=0; it<gen_itmax; it++) {

      xxx = sqrt((*Unit)());
      yyy = 0.5*M_PI*(*Unit)();

      vr = vmax*xxx*cos(yyy);
      vt = vmax*xxx*sin(yyy);
      eee = pot + 0.5*(vr*vr + vt*vt);

      if ((*Unit)() > distf(eee, r*vt)/fmax ) continue;

      if ((*Unit)()<0.5) vr *= -1.0;
    
      phi = 2.0*M_PI*(*Unit)();

      break;
    }

  }

                
  Vector out(1, 6);

  if (it==gen_itmax) {
    cerr << "Velocity selection failed, r=" << r << "\n";
    ierr = 1;
    out.zero();
    return out;
  }

  ierr = 0;
  
  cosp = cos(phi);
  sinp = sin(phi);

  out[1] = r * cosp;
  out[2] = r * sinp;
  out[3] = 0.0;
  out[4] = vr * cosp - vt * sinp;
  out[5] = vr * sinp + vt * cosp;
  out[6] = 0.0;
    
  return out;
}


Vector AxiSymModel::gen_point_2d(double r, int& ierr)
{
  if (!dist_defined) {
    cerr << "AxiSymModel: must define distribution before realizing!\n";
    _exit (-1);
  }

  int it;
  double pot, vmax, xxx, yyy, zzz, vr, vt, eee, fmax;
  double phi, sinp, cosp;

  double Emin = get_pot(get_min_radius());
  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {

    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);

    gen = new ACG(gen_seed, 20);
    Unit = new Uniform(0.0, 1.0, gen);
      
    gen_mass.setsize(1, gen_N);
    gen_rloc.setsize(1, gen_N);
    gen_fmax.setsize(1, gen_N);

    cout << "gen_point_2d[" << ModelID << "]: " << get_min_radius()
	 << ", " << get_max_radius() << endl;

    double dr = (get_max_radius() - get_min_radius())/gen_N;

    for (int i=1; i<=gen_N; i++) {
      gen_rloc[i] = get_min_radius() + dr*((double)i-0.5);
      gen_mass[i] = get_mass(gen_rloc[i]);

      pot = get_pot(gen_rloc[i]);
      vmax = sqrt(2.0*(Emax - pot));

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
      gen_fmax[i] = fmax*1.1;

    }
    gen_firstime = false;
  }

  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  pot = get_pot(r);
  vmax = sqrt(2.0*(Emax - pot));
  
                                // Trial velocity point
    
  for (it=0; it<gen_itmax; it++) {

    xxx = sqrt((*Unit)());
    yyy = 0.5*M_PI*(*Unit)();

    vr = vmax*xxx*cos(yyy);
    vt = vmax*xxx*sin(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if ((*Unit)() > distf(eee, r*vt)/fmax ) continue;
    
    if ((*Unit)()<0.5) vr *= -1.0;
    
    phi = 2.0*M_PI*(*Unit)();

    break;
  }

  
  Vector out(1, 6);

  if (it==gen_itmax) {
    cerr << "Velocity selection failed, r=" << r << "\n";
    ierr = 1;
    out.zero();
    return out;
  }

  ierr = 0;
  
  cosp = cos(phi);
  sinp = sin(phi);

  out[1] = r * cosp;
  out[2] = r * sinp;
  out[3] = 0.0;
  out[4] = vr * cosp - vt * sinp;
  out[5] = vr * sinp + vt * cosp;
  out[6] = 0.0;
    
  return out;
}


Vector AxiSymModel::gen_point_3d(int& ierr)
{
  if (!dist_defined) {
    cerr << "AxiSymModel: must define distribution before realizing!\n";
    _exit (-1);
  }

  double r, pot, vmax, xxx, yyy, vr, vt, eee, vt1, vt2, fmax;
  double phi, sint, cost, sinp, cosp, azi;

  double Emax = get_pot(get_max_radius());

  if (gen_firstime) {
    double zzz;
    // const int numr = 200;
    // const int numj = 200;

    double tol = 1.0e-5;
    double dx = (1.0 - 2.0*tol)/(numr-1);
    double dy = (1.0 - 2.0*tol)/(numj-1);

    gen = new ACG(gen_seed, 20);
    Unit = new Uniform(0.0, 1.0, gen);

    gen_mass.setsize(1, gen_N);
    gen_rloc.setsize(1, gen_N);
    gen_fmax.setsize(1, gen_N);

    double dr = (get_max_radius() - get_min_radius())/gen_N;

    for (int i=1; i<=gen_N; i++) {
      gen_rloc[i] = get_min_radius() + dr*((double)i-0.5);
      gen_mass[i] = get_mass(gen_rloc[i]);

      pot = get_pot(gen_rloc[i]);
      vmax = sqrt(2.0*(Emax - pot));

      fmax = 0.0;
      for (int j=0; j<numr; j++) {
	xxx = tol + dx*j;

	for (int k=0; k<numj; k++) {
	  yyy = tol + dy*k;

	  vr = vmax*xxx;
	  vt = vmax*sqrt((1.0 - xxx*xxx)*yyy);
	  eee = pot + 0.5*(vr*vr + vt*vt);

	  zzz = distf(eee, gen_rloc[i]*vt);
	  fmax = zzz>fmax ? zzz : fmax;
	}
      }
      gen_fmax[i] = fmax*1.1;

    }
    gen_firstime = false;
  }

  r = odd2((*Unit)()*gen_mass[gen_N], gen_mass, gen_rloc, 0);
  fmax = odd2(r, gen_rloc, gen_fmax, 1);
  pot = get_pot(r);
  vmax = sqrt(2.0*(Emax - pot));

                                // Trial velocity point
    
  int it;
  for (it=0; it<gen_itmax; it++) {

    xxx = 2.0*sin(asin((*Unit)())/3.0);
    yyy = (1.0 - xxx*xxx)*(*Unit)();

    vr = vmax*xxx;
    vt = vmax*sqrt(yyy);
    eee = pot + 0.5*(vr*vr + vt*vt);

    if ((*Unit)() > distf(eee, r*vt)/fmax ) continue;

    if ((*Unit)()<0.5) vr *= -1.0;
    
    azi = 2.0*M_PI*(*Unit)();
    vt1 = vt*cos(azi);
    vt2 = vt*sin(azi);

    break;
  }
                
  Vector out(1, 6);

  if (it==gen_itmax) {
    cerr << "Velocity selection failed, r=" << r << "\n";
    ierr = 1;
    out.zero();
    return out;
  }

  ierr = 0;
  
  if ((*Unit)()>=0.5) vr *= -1.0;

  phi = 2.0*M_PI*(*Unit)();
  cost = 2.0*((*Unit)() - 0.5);
  sint = sqrt(1.0 - cost*cost);
  cosp = cos(phi);
  sinp = sin(phi);

  out[1] = r * sint*cosp;
  out[2] = r * sint*sinp;
  out[3] = r * cost;
  out[4] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
  out[5] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
  out[6] = vr * cost      - vt1 * sint;
    
  return out;
}

