// This routine computes orbit information for 1-d slab

// #define NO_SPLINE 1
#undef NO_SPLINE

#include <functional>

#include "numerical.H"
#include "gaussQ.H"
#include "interp.H"
#include "biorth1d.H"
// #include "massmodel1d.H"
#include "Orbit1d.H"
#include "Model1d.H"

using namespace std::literals::complex_literals;

class fct_zmax : public std::function<double(double)>
{
private:

  std::shared_ptr<OneDModel> mdl;
  double E;

public:

  fct_zmax(std::shared_ptr<OneDModel> p, double e) : mdl(p), E(e)  {}

  double operator()(double z)
  {
    return E - mdl->get_pot(z);  
  }
};


OneDOrbit::OneDOrbit(void)
{
  model_defined   = false;
  action_defined  = false;
  freq_defined    = false;
  angle_defined   = false;
  biorth_defined  = false;

  NGRID = 100;
  NINT  = 20;
  NFREQ = 64;
  NREC  = 64;

  dof = 1;
  OrbitID = "OneDimensionalOrbit";
}

OneDOrbit::OneDOrbit(std::shared_ptr<OneDModel> MODEL)
{
  model = MODEL;

  model_defined   = true;
  action_defined  = false;
  freq_defined    = false;
  angle_defined   = false;
  biorth_defined  = false;

  NGRID = 100;
  NINT  = 20;
  NFREQ = 64;
  NREC  = 64;

  dof = 1;
  OrbitID = "OneDimensionalOrbit";
}

OneDOrbit::OneDOrbit(std::shared_ptr<OneDModel> MODEL, double ENERGY, double VPARA)
{
  model = MODEL;
  energy = ENERGY;
  vp = VPARA;

  model_defined   = true;
  action_defined  = false;
  freq_defined    = false;
  angle_defined   = false;
  biorth_defined  = false;

  NGRID = 100;
  NINT  = 20;
  NFREQ = 64;
  NREC  = 64;

  dof = 1;
  OrbitID = "OneDimensionalOrbit";
}

OneDOrbit::OneDOrbit(OneDOrbit& p)
{
  model  = p.model;
  energy = p.energy;
  vp     = p.vp;

  zmax   = p.zmax;
  freq   = p.freq;
  action = p.action;

  model  = p.model;
  biorth = p.biorth;

  model_defined   = p.model_defined;
  action_defined  = p.action_defined;
  freq_defined    = p.freq_defined;
  angle_defined   = p.angle_defined;
  biorth_defined  = p.biorth_defined;

  NGRID = p.NGRID;
  NINT  = p.NINT;
  NFREQ = p.NFREQ;
  NREC  = p.NREC;

  lkw   = p.lkw;
  pot_cache = p.pot_cache;
  pot_eval = p.pot_eval;

  dof = 1;
  OrbitID = "OneDimensionalOrbit";
}


void OneDOrbit::new_orbit(double ENERGY, double VPARA)
{
  energy = ENERGY;
  vp = VPARA;

  freq_defined  = false;
  angle_defined = false;
  if (pot_cache) pot_eval.clear();
}


// Perform integral using centered rectangles

void OneDOrbit::compute_freq(void)
{
  const double E = energy;
  auto fct = [mdl=model, E](double z) { return E - mdl->get_pot(z); };

  const double mmin = model->get_min_radius();
  const double mmax = model->get_max_radius();

  zmax = zbrent(fct, mmin, 4.0*mmax, TOL);

  if (energy - model->get_pot(mmin) < TOL || zmax < TOL) {
    freq = std::sqrt(4.0*M_PI*model->get_density(mmin));
    action = 0.0;
    freq_defined = true;
    action_defined = true;
    return;
  }

  double dt = M_PI/(2.0*NFREQ);

  double ansf = 0.0;
  double ansi = 0.0;
  for (int j=0; j<NFREQ; j++) {
    double t = dt*(0.5+(double)j);
    double z = zmax*sin(t);
    ansf += cos(t) / sqrt( 2.0*(E-model->get_pot(z)) );
    ansi += cos(t) * sqrt( 2.0*(E-model->get_pot(z)) );
  }

  freq = 0.5*M_PI/(ansf*zmax*dt);
  action = 4.0*ansi*zmax*dt;

  freq_defined = true;
}
  


void OneDOrbit::compute_angles(void)
{
  if (!freq_defined) compute_freq();

  double E = energy;

  // Gaussian integration (prevent recomputation)
  //
  if (lkw == nullptr or lkw->get_n() != NINT)
    lkw = std::make_shared<LegeQuad>(NINT);

  angle_grid.size = NGRID;
  angle_grid.time.resize(NGRID+1);
  angle_grid.angl.resize(2, NGRID+1);
  angle_grid.hght.resize(2, NGRID+1);

  double z, dz = zmax/NGRID;
 
  angle_grid.time(0) = 0.0;
  angle_grid.angl(0, 0) = 0.0;
  angle_grid.hght(0, 0) = 0.0;

  double ans, tmin, tmax, t;

  for (int i=1; i<=NGRID; i++) {

    z = angle_grid.hght(0, i) = dz*i;

    tmin = asin(dz*(i-1)/zmax);
    if (dz*i/zmax>=1.0)
      tmax = 0.5*M_PI;
    else
      tmax = asin(dz*i/zmax);

    ans = 0.0;
    for (int j=0; j<NINT; j++) {
      t = tmin + (tmax - tmin) * lkw->knot(j);
      double denom = 2.0*(E-model->get_pot(zmax*sin(t)));
      if (denom<=0.0) {
	double dt = 0.5*M_PI - t;
	denom = model->get_dpot(zmax) * zmax * dt*dt;
      }
      if (denom > 0.0)
	ans += lkw->weight(j) * cos(t) / sqrt(denom);
      else
	std::cerr << "Warning: negative denominator in compute_angles: "
		  << "E=" << E << " pot=" << model->get_pot(zmax*sin(t))
		  << " denom=" << denom << std::endl;
    }
    ans *= zmax*(tmax - tmin);

    angle_grid.time(i) = angle_grid.time(i-1) + ans;
    angle_grid.angl(0, i) = freq * angle_grid.time(i);
  }

#ifndef NO_SPLINE
  Eigen::VectorXd work(NGRID+1);
  Spline(angle_grid.time, angle_grid.angl.row(0), -1.0e30, -1.0e30, work);
  angle_grid.angl.row(1) = work;

  Spline(angle_grid.time, angle_grid.hght.row(0), -1.0e30, -1.0e30, work);
  angle_grid.hght.row(1) =  work;
#endif

  angle_defined = true;
}
  

double OneDOrbit::get_angle(const Type f, const double T)
{
  switch (f) {
  case Type::Angle:
    return get_angle(1, T);
  case Type::Height:
    return get_angle(2, T);
  case Type::Velocity:
    return get_angle(3, T);
  default:
    throw "OneDOrbit::get_angle: invalid angle type";
  }
}

double OneDOrbit::get_angle(const int i, const double T)
{
  if (!angle_defined) compute_angles();

  // Compute the time modulo the quarter period and identify the
  // quadrant
  //
  double time, qperiod = 0.5*M_PI/freq;

  if (T<0.0)
    time = T + 4.0*qperiod * (1.0 + (int)(0.25*fabs(T)/qperiod));
  else
    time = T - 4.0*qperiod * (int)(0.25*fabs(T)/qperiod);

  int phase = (int)(time/qperiod);
  double ptime = time - phase*qperiod;

  phase = phase % 4;

  double zsign = 1.0;
  if (phase > 1) zsign = -1.0;
  if (phase == 1 || phase == 3) ptime = qperiod - ptime;

  double ans, vel;

  switch (i) {
  case 3:
#ifndef NO_SPLINE
    Splint1(angle_grid.time, angle_grid.hght.row(0), angle_grid.hght.row(1), ptime, ans);
#else
    ans = odd2(ptime, angle_grid.time, angle_grid.hght);
#endif
    vel = std::sqrt(2.0*std::fabs((energy - model->get_pot(ans))));
    if (phase == 1 || phase == 2) vel = -vel;
    /*
    if (std::fabs(energy - model->get_pot(ans) - 0.5*vel*vel) > 1.0e-6) {
      std::cout << "Velocity error: energy=" << energy
		<< " pot=" << model->get_pot(ans)
		<< " kinetic=" << 0.5*vel*vel
		<< " error=" << energy - model->get_pot(ans) - 0.5*vel*vel
		<< std::endl;return vel;
    }
    */
    return vel;
  case 2:
#ifndef NO_SPLINE
    Splint1(angle_grid.time, angle_grid.hght.row(0), angle_grid.hght.row(1), ptime, ans);
#else
    ans = odd2(ptime, angle_grid.time, angle_grid.hght);
#endif
    ans *= zsign;
    break;
  case 1:
  default:
#ifndef NO_SPLINE
    Splint1(angle_grid.time, angle_grid.angl.row(0), angle_grid.angl.row(1), ptime, ans);
#else
    ans = odd2(ptime, angle_grid.time, angle_grid.angl);
#endif
    switch (phase) {
    case 1:
      ans = M_PI - ans;
      break;
    case 2:
      ans = M_PI + ans;
      break;
    case 3:
      ans = 2.0*M_PI - ans;
      break;
    }
    break;
  }

  return ans;
}

double OneDOrbit::pot_trans(const int k, double (*func)(double z))
{
  double freq = get_freq();
  double qperiod = 0.5*M_PI/freq;

				// Use centered rectangles
  double t, dt = 2.0*qperiod/NREC;
  double ans = 0.0;
  for (int i=0; i<NREC; i++) {
    t = -qperiod + dt*( (double)i + 0.5 );
    ans += cos(t*freq*k) * (*func)(get_angle(2, t));
  }

  ans *= dt*freq/M_PI;
  
  return ans;
}


double OneDOrbit::pot_trans(const int k, const int n)
{
  if (!biorth_defined) 
    bomb ("Must define biorthogonal set before calling pot_trans(int, int)");

  double freq = get_freq();
  double qperiod = 0.5*M_PI/freq;

				// Use centered rectangles
  double t, dt = 2.0*qperiod/NREC;
  double ans = 0.0;

  for (int i=0; i<NREC; i++) {
    t = -qperiod + dt*( (double)i + 0.5 );
    ans += cos(t*freq*k) * biorth->potl(n, 0, get_angle(2, t));
  }

  ans *= dt*freq/M_PI;
  
  return ans;
}

void OneDOrbit::pot_trans(const int k, Eigen::VectorXd& t)
{
  if (!biorth_defined) 
    bomb ("Must define biorthogonal set before calling pot_trans(int, int)");

  double freq = get_freq();
  double qperiod = 0.5*M_PI/freq;

				// Use centered rectangles
  double T, dT = 2.0*qperiod/NREC;
  Eigen::VectorXd p(t.size()); t.setZero();

  for (int i=0; i<NREC; i++) {
    T = -qperiod + dT*( (double)i + 0.5 );
    // biorth->potl_vec(0, 0, get_angle(2, T), p);
    biorth->potl(0, 0, get_angle(2, T), p);
    t += cos(T*freq*k) * p;
  }

  t *= dT*freq/M_PI;
}

std::complex<double> OneDOrbit::pot_trans_complex
(const int k, std::function<double(double)> func)
{
  double freq = get_freq();
  double hperiod = M_PI/freq;

				// Use centered rectangles
  double t, dt = 2.0*hperiod/NREC;
  Complex ans = 0.0;

  for (int i=0; i<NREC; i++) {
    t = -hperiod + dt*( (double)i + 0.5 );
    ans += exp(-1i*t*freq*static_cast<double>(k)) * func(get_angle(2, t));
  }

  ans *= 0.5*dt*freq/M_PI;
  
  return ans;
}


std::complex<double> OneDOrbit::pot_trans_complex(const int k, const int n)
{
  if (!biorth_defined) 
    bomb ("Must define biorthogonal set before calling "
	  "pot_trans_complex(int, int)");

  double freq = get_freq();
  double hperiod = M_PI/freq;

  // Result
  //
  Complex ans = 0.0;

  // Use centered rectangles
  //
  double dt = 2.0*hperiod/NREC;
  for (int i=0; i<NREC; i++) {
    double t = -hperiod + dt*( (double)i + 0.5 );
    ans += exp(-1i*t*freq*static_cast<double>(k)) *
      biorth->potl(n, 0, get_angle(2, t));
  }

  return ans * 0.5*dt*freq/M_PI;
}

std::complex<double> OneDOrbit::pot_trans_complex_check(const int k, const int n)
{
  if (!biorth_defined) 
    bomb ("Must define biorthogonal set before calling "
	  "pot_trans_complex_check(int, int)");

  double freq = get_freq();

  // Result
  Complex ans = 0.0;

  // Use centered rectangles
  double dw = 2.0*M_PI/NREC;
  
  for (int i=0; i<NREC; i++) {
    double w = dw*(0.5+(double)i);
    ans += exp(-1i*w*static_cast<double>(k)) *
      biorth->potl(n, 0, get_angle(2, w/freq));
  }

  return ans/double(NREC);
}

std::complex<double> OneDOrbit::pot_trans_complex_check
(const int k, const int nx, const int n)
{
  if (!biorth_defined) 
    bomb ("Must define biorthogonal set before calling "
	  "pot_trans_complex_check(int, int)");

  double freq = get_freq();

  // Result
  Complex ans = 0.0;

  // Use centered rectangles
  double dw = 2.0*M_PI/NREC;
  
  for (int i=0; i<NREC; i++) {
    double w = dw*(0.5+(double)i);
    ans += exp(-1i*w*static_cast<double>(k)) *
      biorth->potl(n, nx, get_angle(OneDOrbit::Type::Height, w/freq));
  }

  return ans/double(NREC);
}

void OneDOrbit::set_pot_cache(bool val, int norder)
{
  pot_eval.clear();
  pot_cache = true;
  if (val == false) {
    pot_cache = false;
    return;
  }

  if (!biorth_defined) 
    bomb ("Must define biorthogonal set before calling pot_trans(int, int)");

  double freq = get_freq();
  double hperiod = M_PI/freq;

				// Use centered rectangles
  double T, dT = 2.0*hperiod/NREC;
  Eigen::VectorXd p(norder);

  pot_eval.resize(NREC);
  for (int i=0; i<NREC; i++) {
    T = -hperiod + dT*( (double)i + 0.5 );
    biorth->potl(0, 0, get_angle(2, T), p);
    pot_eval[i] = p;
  }
}

void OneDOrbit::pot_trans_complex(const int k, Eigen::VectorXcd& t)
{
  if (!biorth_defined) 
    bomb ("Must define biorthogonal set before calling pot_trans(int, int)");

  double freq = get_freq();
  double hperiod = M_PI/freq;

				// Use centered rectangles
  double dT = 2.0*hperiod/NREC;

  Eigen::VectorXd p(t.size());	// Evaluation vector for potential transform
  t.setZero();			// Zero the return vector

  for (int i=0; i<NREC; i++) {
    double T = -hperiod + dT*( (double)i + 0.5 );
    Complex factor =  exp(-1i*T*freq*static_cast<double>(k));
    if (pot_cache) {
      t += factor * pot_eval[i].cast<std::complex<double> >();
    } else {
      biorth->potl(0, 0, get_angle(2, T), p);
      t += factor * p.cast<std::complex<double> >();
    }
  }

  t *= 0.5*dT*freq/M_PI;
}

void OneDOrbit::pot_trans_complex(const int nx, const int ny, const int nz,
				  Eigen::VectorXcd& t)
{
  if (!biorth_defined) 
    bomb ("Must define biorthogonal set before calling pot_trans(int, int)");

  double freq    = get_freq();
  double hperiod = M_PI/freq;
				// Use centered rectangles
  double dT = 2.0*hperiod/NREC;
  
  t.setZero();			// Zero the return vector

  Eigen::VectorXd p(t.size());

  for (int i=0; i<NREC; i++) {
    double T = -hperiod + dT*( (double)i + 0.5 );
    biorth->potl(nx, ny, get_angle(2, T), p);
    t += exp(-1i*T*freq*static_cast<double>(nz)) *
      p.cast<std::complex<double> >();
  }

  t *= 0.5*dT*freq/M_PI;
}

// This version is for checking the results of the above routine
// 
void OneDOrbit::pot_trans_complex_check
(const int nx, const int ny, const int nz, Eigen::VectorXcd& t)
{
  if (!biorth_defined) 
    bomb ("Must define biorthogonal set before calling pot_trans(int, int)");

  double freq    = get_freq();
  double hperiod = M_PI/freq;

  // Zero the return vector
  //
  t.setZero();

  // Use centered rectangles
  //
  Eigen::VectorXd p(t.size());
  double dw = 2.0*M_PI/NREC;
  for (int i=0; i<NREC; i++) {
    double w = dw*(0.5+(double)i);
    double z = get_angle(OneDOrbit::Type::Height, w/freq);
    biorth->potl(nx, ny, z, p);
    t += exp(-1i*w*static_cast<double>(nz)) * p.cast<std::complex<double>>();
  }

  t /= NREC;
}

