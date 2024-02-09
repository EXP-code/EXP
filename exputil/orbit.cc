
#pragma implementation

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>

#include <massmodel.H>
#include <interp.H>
#include <orbit.H>

void RegularOrbit::bomb(const char *s)
{
  std::ostringstream msg;
  msg << "ERROR from " << OrbitID << ": " << s;
  throw std::runtime_error(msg.str());
}

void RegularOrbit::bomb(const string& s) {
  std::ostringstream msg;
  msg << "ERROR from " << OrbitID << ": " << s;
  throw std::runtime_error(msg.str());
}

void RegularOrbit::warn(const string& fct, const string& msg)
{
  cerr << fct << ": from " << OrbitID << ": " << msg << endl;
}

bool SphericalOrbit::guard = true;

SphericalOrbit::SphericalOrbit(void) {
  dof = 3;
  OrbitID = "SphericalOrbit";

  energy = 0.0;
  kappa = 0.0;
  beta = 0.0;

  jmax = 0.0;
  r_peri = 0.0;
  r_apo = 0.0;
  r_circ = 0.0;

  model = 0;

  model_defined  = false;
  action_defined = false;
  freq_defined   = false;
  angle_defined  = false;
  biorth_defined = false;

  set_numerical_params();
}

SphericalOrbit::SphericalOrbit(AxiSymModPtr Model)
{
  dof = Model->dof();
  OrbitID = "SphericalOrbit";

  energy = 0.0;
  kappa = 0.0;
  beta = 0.0;

  jmax = 0.0;
  r_peri = 0.0;
  r_apo = 0.0;
  r_circ = 0.0;

  model = Model;
  model_defined  = true;

  action_defined = false;
  freq_defined   = false;
  angle_defined  = false;
  biorth_defined = false;

  set_numerical_params();
}

SphericalOrbit::SphericalOrbit(AxiSymModPtr Model,
			       double Energy, double Kappa, double Beta)
{
  dof = Model->dof();
  OrbitID = "SphericalOrbit";

  energy = 0.0;
  kappa = 0.0;
  beta = 0.0;

  jmax = 0.0;
  r_peri = 0.0;
  r_apo = 0.0;
  r_circ = 0.0;

  model = Model;
  model_defined = true;

  set_numerical_params();

  new_orbit(Energy, Kappa, Beta);
}


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

  model_defined  = t.model_defined;
  action_defined = t.action_defined;
  freq_defined   = t.freq_defined;
  angle_defined  = t.angle_defined;
}

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

  return *this;
}

void SphericalOrbit::new_orbit(double Energy, double Kappa, double Beta)
{

  if (Kappa < 0.0 || Kappa > 1.0) {
    ostringstream ost;
    ost << "[new_orbit] illegal value of Kappa: " << Kappa;
    bomb(ost.str().c_str());
  }

  //  if (Energy < model->get_pot(model->get_min_radius()) || 
  //      Energy > model->get_pot(model->get_max_radius())) {

  if (guard && Energy < model->get_pot(model->get_min_radius())) {
    ostringstream ost;
    ost << "[new orbit] illegal value of Energy: " << Energy;
    ost << "  Emin[r=" << model->get_min_radius() << "] = " 
	<< model->get_pot(model->get_min_radius());
    bomb(ost.str().c_str());
  }


  energy = Energy;
  kappa  = Kappa;
  beta   = Beta;

  action_defined = false;
  freq_defined   = false;
  angle_defined  = false;
  biorth_defined = false;
}

double SphericalOrbit::get_angle(const int i, const double time)
{
  
  if (!freq_defined ) compute_freq();
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

  double w1 = freq[0]*time;
  double w2 = freq[1]*time;
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

  Eigen::VectorXd X(angle_grid.w1.row(0));
  
  switch (i) {

  case 3:
#ifdef SPLINE
    Splint1(X,
	    Eigen::VectorXd(angle_grid.t.row(0)),
	    Eigen::VectorXd(angle_grid.t.row(1)), w1, ans);
#else
    ans = odd2(w1, X, Eigen::VectorXd(angle_grid.t.row(0)));
#endif
    if (branch) ans = M_PI - ans;
    break;

  case 4:
#ifdef SPLINE    
    Splint1(X,
	    Eigen::VectorXd(angle_grid.dw1dt.row(0)),
	    Eigen::VectorXd(angle_grid.dw1dt.row(1)), w1, ans);
#else
    ans = odd2(w1, X, Eigen::VectorXd(angle_grid.dw1dt.row(0)));
#endif
    break;

  case 5:
#ifdef SPLINE
    Splint1(X,
	    Eigen::VectorXd(angle_grid.f.row(0)),
	    Eigen::VectorXd(angle_grid.f.row(1)), w1, ans);
#else
    ans = odd2(w1, X, Eigen::VectorXd(angle_grid.f.row(0)));
#endif
    if (branch) ans *= -1.0;
    break;

  case 6:
#ifdef SPLINE
    Splint1(X,
	    Eigen::VectorXd(angle_grid.r.row(0)),
	    Eigen::VectorXd(angle_grid.r.row(1)), w1, ans);
#else
    ans = odd2(w1, X, Eigen::VectorXd(angle_grid.r.row(0)));
#endif
    break;

  case 7:
#ifdef SPLINE
    Splint1(X,
	    Eigen::VectorXd(angle_grid.f.row(0)),
	    Eigen::VectorXd(angle_grid.f.row(1)), w1, ans);
#else
    ans = odd2(w1, X, Eigen::VectorXd(angle_grid.f.row(0)));
#endif
    if (branch) ans *= -1.0;
    ans = w2 - ans;
    break;

  }

  return ans;
}

// Get radial angle for a given radius in the range [0, 2*pi]
//
double SphericalOrbit::get_w1(double r, double vr)
{
  
  if (!freq_defined)  compute_freq();
  if (!angle_defined) compute_angles();

  // Check bounds
  if (r<r_peri or r>r_apo) return std::numeric_limits<double>::infinity();

  // End points as special cases
  if (r<r_peri) return 0.0;
  if (r>r_apo ) return M_PI;

  // Interpolate for the angle for peri-to-apo branch
  double ang = odd2(r,
		    Eigen::VectorXd(angle_grid.r. row(0)),
		    Eigen::VectorXd(angle_grid.w1.row(0)) );

  // Use apo-to-peri branch?
  if (vr<0.0) ang = 2.0*M_PI - ang;

  return ang;
}
