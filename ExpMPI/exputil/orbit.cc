
#pragma implementation

const char rcsid[] = "$Id$";

#include <math.h>
#include <string>
#include <iostream.h>
#include <strstream.h>
#include <stdlib.h>

#include <Vector.h>
#include <orbit.h>
#include <interp.h>

SphericalOrbit::SphericalOrbit(void) {
  dof = 3;
  OrbitID = "SphericalOrbit";
  Gkn = NULL;

  model_defined = false;
  action_defined = false;
  freq_defined = false;
  angle_defined = false;
  biorth_defined = false;

  set_numerical_params();
}

SphericalOrbit::SphericalOrbit(AxiSymModel *Model)
{
  dof = Model->dof();
  OrbitID = "SphericalOrbit";
  Gkn = NULL;

  model = Model;
  model_defined = true;

  action_defined = false;
  freq_defined = false;
  angle_defined = false;
  biorth_defined = false;

  set_numerical_params();
}

SphericalOrbit::SphericalOrbit(AxiSymModel *Model,
			       double Energy, double Kappa, double Beta)
{
  dof = Model->dof();
  OrbitID = "SphericalOrbit";
  Gkn = NULL;

  model = Model;
  model_defined = true;

  set_numerical_params();

  new_orbit(Energy, Kappa, Beta);
}


SphericalOrbit::~SphericalOrbit()
{
  if (Gkn) {
    delete Gkn;
  }
}

/*
SphericalOrbit::SphericalOrbit(const SphericalOrbit &t)
{
  energy = t.energy;
  kappa = t.kappa;
  beta = t.beta;

  jmax = t.jmax;
  r_peri = t.r_peri;
  r_apo = t.r_apo;
  r_circ = t.r_circ;

  model = t.model;

  model_defined = t.model_defined;
  action_defined = t.action_defined;
  freq_defined = t.freq_defined;
  angle_defined = t.angle_defined;
}
*/


SphericalOrbit &SphericalOrbit::operator=(const SphericalOrbit &t)
{
  energy = t.energy;
  kappa = t.kappa;
  beta = t.beta;

  jmax = t.jmax;
  r_peri = t.r_peri;
  r_apo = t.r_apo;
  r_circ = t.r_circ;

  model = t.model;

  model_defined = t.model_defined;
  action_defined = t.action_defined;
  freq_defined = t.freq_defined;
  angle_defined = t.angle_defined;

  if (t.Gkn) {
    Gkn = new LegeQuad(t.Gkn->get_n());
  }

  return *this;
}

void SphericalOrbit::new_orbit(double Energy, double Kappa, double Beta)
{

  if (Kappa < 0.0 || Kappa > 1.0) {
    char* p = new char[128];
    ostrstream ost(p,128);
    ost << "illegal value of Kappa: " << Kappa << '\0';
    bomb(p);
  }

//  if (Energy < model->get_pot(model->get_min_radius()) || 
//      Energy > model->get_pot(model->get_max_radius())) {
  if (Energy < model->get_pot(model->get_min_radius())) {
    char* p = new char[128];
    ostrstream ost(p,128);
    ost << "illegal value of Energy: " << Energy;
    ost << "  Emin[r=" << model->get_min_radius() << "] = " 
	<< model->get_pot(model->get_min_radius()) << '\0';    
    bomb(p);
  }


  energy = Energy;
  kappa = Kappa;
  beta = Beta;

  action_defined = false;
  freq_defined = false;
  angle_defined = false;
  biorth_defined = false;
}

double SphericalOrbit::get_angle(const int i, const double time)
{
  
  if (!freq_defined) compute_freq();
  if (!angle_defined) compute_angles();

  /*
     i     val
     -     ---
     1      w1
     2      w2
     3      t
     4      dw1dt
     5      f
     6      r
     7      phi
   */

  double w1 = freq[1]*time;
  double w2 = freq[2]*time;
  double ans = 0.0;
  int branch=0;

  if (i==1) return w1;
  if (i==2) return w2;

  if (w1 < 0.0) w1 += 2.0*M_PI*(1.0 + abs((int)(w1/(2.0*M_PI))));
  w1 -= 2.0*M_PI*(int)(w1/(2.0*M_PI));
  if (w1 > M_PI) {
    branch=1;
    w1 = 2.0*M_PI - w1;
  }

  switch (i) {

  case 3:
#ifdef SPLINE
    Splint1(angle_grid.w1[1], angle_grid.t[1], angle_grid.t[2], w1, ans);
#else
    ans = odd2(w1, angle_grid.w1[1], angle_grid.t[1]);
#endif
    if (branch) ans = M_PI - ans;
    break;

  case 4:
#ifdef SPLINE    
    Splint1(angle_grid.w1[1], angle_grid.dw1dt[1], angle_grid.dw1dt[2], w1, ans);
#else
    ans = odd2(w1, angle_grid.w1[1], angle_grid.dw1dt[1]);
#endif
    break;

  case 5:
#ifdef SPLINE
    Splint1(angle_grid.w1[1], angle_grid.f[1], angle_grid.f[2], w1, ans);
#else
    ans = odd2(w1, angle_grid.w1[1], angle_grid.f[1]);
#endif
    if (branch) ans *= -1.0;
    break;

  case 6:
#ifdef SPLINE
    Splint1(angle_grid.w1[1], angle_grid.r[1], angle_grid.r[2], w1, ans);
#else
    ans = odd2(w1, angle_grid.w1[1], angle_grid.r[1]);
#endif
    break;

  case 7:
#ifdef SPLINE
    Splint1(angle_grid.w1[1], angle_grid.f[1], angle_grid.f[2], w1, ans);
#else
    ans = odd2(w1, angle_grid.w1[1], angle_grid.f[1]);
#endif
    if (branch) ans *= -1.0;
    ans = w2 - ans;
    break;

  }

  return ans;

}

