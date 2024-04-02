/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  These routines computes massmodels for 1-d slab
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *  x        as above
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/88
 *
 ***************************************************************************/


#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include <gaussQ.H>
#include <interp.H>

#include <massmodel1d.H>
#include <RK4.H>

double Sech2::HMAX=1.0e6;
double Sech2mu::HMAX=1.0e6;
double Sech2Halo::HMAX=3.0e1;

OneDModelTable::OneDModelTable(std::string filename, int PARM)
{
  int		i, imodel;
  char		line[144];
  double	radius;

  std::istringstream iline(line);

  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Error opening: " << filename << " . . . quitting"
	      << std::endl;
    exit(-1);
  }
  
				// Read header

  in.getline((char *)line, 144);
  while (string(line).find_first_of("!#") != string::npos) 
    in.getline((char *)line, 144);
  

				// Assign space for model
  iline >> imodel;
  
  std::vector<double> z(imodel), d(imodel), m(imodel), p(imodel);

  for (int i=0; i<imodel; i++) {

    in.getline(line, 144);
    istringstream iline(line);

    iline >> z[i];
    iline >> d[i];
    iline >> m[i];
    iline >> p[i];
  }

  if (PARM) {
    in.getline(line, 144);

    istringstream ins(line);
    string word;
    while (1) {
      ins >> word;
      if (ins.good()) params.push_back(atof(word.c_str()));
    }
  }

				// Compute splines
  dens = Spline1d(z, d, 0.0, 0.0);
  mass = Spline1d(z, m, 0.0, 0.0);
  pot  = Spline1d(z, p, 0.0, 2.0*M_PI*m[imodel-1]);

  even = 0;
  if (fabs((z[1]-z[0]) - (z[2] - z[1])) < 1.0e-6) even = 1;

}

double OneDModelTable::get_mass(const double z)
{
  double zz = fabs(z);
  if (zz>mass.xhi()) return mass.eval(mass.xhi());

  return mass.eval(zz);
}

double OneDModelTable::get_density(const double z)
{
  double zz = fabs(z);
  if (z>dens.xhi()) return 0.0;

  return dens.eval(zz);
}

double OneDModelTable::get_pot(const double z)
{
  double zz = fabs(z);

  if (zz>pot.xhi()) {
    double z0 = pot.xlo(), z1 = pot.xhi();
    double dz = (z1 - z0)*0.01;
    double p1 = pot.eval(z1), pm = pot.eval(z1-dz);
    return p1 + (p1 - pm)/dz * (zz - z1);
  }

  return pot.eval(zz);
}

double OneDModelTable::get_dpot(const double z)
{
  double ans, dum;
  double zz = fabs(z);

  if (zz>pot.xhi()) {
    double z0 = pot.xlo(), z1 = pot.xhi();
    double dz = (z1 - z0)*0.01;
    double p1 = pot.eval(z1), pm = pot.eval(z1-dz);
    ans = (p1 - pm)/dz;
  } else {
    ans = pot.deriv(zz);
  }
  
  return ans*z/(zz+1.0e-18);
}

std::tuple<double, double> OneDModelTable::get_pot_dpot(const double z)
{
  double zz = fabs(z), ur, dur;

  if (zz>pot.xhi()) {
    double z0 = pot.xlo(), z1 = pot.xhi();
    double dz = (z1 - z0)*0.01;
    double p1 = pot.eval(z1), pm = pot.eval(z1-dz);
    ur  = p1 + (p1 - pm)/dz * (zz - z1);
    dur = (p1 - pm)/dz;
  } else {
    ur  = pot.eval(zz);
    dur = pot.deriv(zz);
  }

  dur *= z/(zz+1.0e-18);

  return {ur, dur};
}

double OneDModelTable::get_dpot2(const double z)
{
  double ans;
  double zz = fabs(z);

  if (zz>pot.xhi()) 
    ans = 0.0;
  else
    ans = dens.eval(zz);

  return 4.0*M_PI*ans;
}

void LowIso::setup_model(void)
{
  w0 = params[0];
  Bfac = params[1];
  betak = 1.0/params[2];
  gammak = params[3] / (8.0*sqrt(2.0)*M_PI);

  normx = 1.0/sqrt(2.0*M_PI*dispx);
}

double LowIso::distf(const double E, const double V)
{
  if (E>0.0)
    return 0.0;
  else
    return gammak*(exp(-betak*E) - 1.0) * normx * exp(-0.5*V*V/dispx);
}

double LowIso::dfde(const double E, const double V)
{
  if (E>0.0)
    return 0.0;
  else
    return -betak*gammak*exp(-betak*E) * normx * exp(-0.5*V*V/dispx);;
}

double LowIso::dfdv(const double E, const double V)
{
  if (E>0.0)
    return 0.0;
  else
    return -normx*exp(-0.5*V*V)*V/dispx * gammak*(exp(-betak*E) - 1.0);
}


double LowIso::get_pot(const double z)
{
  double zmax = get_max_radius();

  if (fabs(z) > zmax)
    return OneDModelTable::get_pot(zmax) + 2.0*M_PI*get_mass(zmax)*(z-zmax)
    + Bfac*(z*z - zmax*zmax);
  else
    return OneDModelTable::get_pot(z);
}

double LowIso::get_dpot(const double z)
{
  double zmax = get_max_radius();

  if (fabs(z) > zmax)
    return  2.0*M_PI*get_mass(zmax)*z/(fabs(z)+1.0e-16) + 2.0*Bfac*z;
  else
    return OneDModelTable::get_dpot(z);
}

double LowIso::get_dpot2(const double z)
{
  return OneDModelTable::get_dpot2(z) + 2.0*Bfac;
}

std::tuple<double, double> LowIso::get_pot_dpot(const double z)
{
  double zmax = get_max_radius();

  if (fabs(z) > zmax) {
    double ur = OneDModelTable::get_pot(zmax) + 2.0*M_PI*get_mass(zmax)*(z-zmax)
      + Bfac*(z*z - zmax*zmax);
    double dur = 2.0*M_PI*get_mass(zmax)*z/(fabs(z)+1.0e-16) + 2.0*Bfac*z;
    return {ur, dur};
  }
  else
    return OneDModelTable::get_pot_dpot(z);
}

static double halo_fac, halo_height;

bool Sech2Halo::MU = true;
int Sech2Halo::NTABLE = 400;
double Sech2Halo::OFFSET = 1.0e-4;

Sech2Halo::Sech2Halo(const double DISPZ, const double DRATIO, 
		     const double HRATIO, const double DISPX) 
{
  dispz = DISPZ;
  dispx = DISPX;
  dratio = DRATIO;
  hratio = HRATIO;

  h = 1.0;
  rho0 = 1.0;
  hh = hratio * h;
  rho0h = dratio * rho0;
    
  dist_defined = false;

  reset();
}
  

void Sech2Halo::reset()
{
  even = 1;
  hmax = HMAX*h;
  halo_height = hh/h;
  halo_fac = halo_height*halo_height * rho0h/rho0;
  
  std::vector<double> x(NTABLE), d(NTABLE), m(NTABLE), p(NTABLE);

  double w1, dz = hmax/(NTABLE-1);
  double step = dz/100.0;
  
  x[0] = 0.0;
  d[0] = rho0;
  m[0] = p[0] = 0.0;

  double z, lastz=OFFSET, tol=1.0e-6;

  w1 = 2.0*halo_fac*log(cosh(lastz/halo_height));
  Eigen::VectorXd y(2);
  y(0) = exp(-w1)*lastz*lastz*(1.0 - exp(-w1)*lastz*lastz/6.0);
  y(1) = 2.0*exp(-w1)*lastz*(1.0 - exp(-w1)*lastz*lastz/2.0);

  Eigen::VectorXd dy(2);

  auto halo_derivs = [&](double x, ODE::RK4::Vref y)
  {
    dy(0) = y(1);
    double w1 = 2.0*halo_fac*log(cosh(x/halo_height));
    dy(1) = 2.0*exp(-w1-y(0));
    return ODE::RK4::Vref(dy);
  };

  ODE::RK4 rk4(halo_derivs);

  for (int i=1; i<NTABLE; i++) {

    x[i] = z = dz*i;

    for (int j=0; j<100; j++) {
      lastz = rk4.step(y, lastz, std::min(step, z-lastz));
      if (lastz >= z) break;
    }

    p[i] = y[0] * dispz;
    w1 = 2.0*halo_fac*log(cosh(z/halo_height));
    d[i] = rho0*exp(-w1-y[0]);
    m[i] = y[1];

    lastz = z;
  }

				// Scaling
  if (MU) {
    double mfac = m[NTABLE-1];
    for (auto & v : d) v *= 2.0*M_PI/(dispz*mfac*mfac)/rho0;
    for (auto & m : d) m /= mfac;
    rho0 = 2.0*M_PI/(dispz*mfac*mfac);
  }
  else {
    for (auto & v : d) v /= rho0;
    for (auto & v : m) v *= sqrt(0.5*dispz/M_PI);

    rho0 = 1.0;
  }

  h = sqrt(dispz/(2.0*M_PI*rho0));
  rho0h = rho0 * dratio;
  hh = h * hratio;

  for (auto & v : x) v *= h;

  dens = Spline1d(x, d, 0.0, 0.0);
  mass = Spline1d(x, m, 0.0, 0.0);
  pot  = Spline1d(x, p, 0.0, p[NTABLE-1]/x[NTABLE-1]);

  norm = d[0]/( sqrt(2.0*M_PI*dispz) );
  model_computed = true;
  dist_defined = true;
}
