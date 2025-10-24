#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "orbit.H"
#include "massmodel.H"

#include "SimAnn.H"
#include "FindOrb.H"


void FindOrb::mapvars(std::vector<double>& x, double& ee, double& kk)
{
  ee = Emin + (Emax-Emin)*(atan(x[0])/M_PI + 0.5);
  kk = Kmin + (Kmax-Kmin)*(atan(x[1])/M_PI + 0.5);
}


double FindOrb::KMIN=0.005;
double FindOrb::KMAX=0.995;
int    FindOrb::MAXIT=500;
bool   FindOrb::MELT=true;
double FindOrb::RATE=0.25;
double FindOrb::EFAC=1.0;
double FindOrb::T0=1.0;

FindOrb::FindOrb(std::shared_ptr<AxiSymModel> mod, double PERI, double APO)
{
  halo_model = mod;
  peri = PERI;
  apo =  APO;

  orb = std::make_shared<SphericalOrbit>(halo_model);
  Emin = halo_model->get_pot(halo_model->get_min_radius());
  Emax = halo_model->get_pot(EFAC*halo_model->get_max_radius());
  Kmin = KMIN;
  Kmax = KMAX;
}

FindOrb::~FindOrb()
{
  // Nothing
}

double FindOrb::operator()(std::vector<double>& ek)
{
  double ee, kk;

  mapvars(ek, ee, kk);

  orb->new_orbit(ee, kk);

  double APO = orb->apo() - apo;
  double PERI = orb->peri() - peri;

  double ret =  APO*APO + PERI*PERI;

  return ret;
}


OrbValues FindOrb::Anneal()
{
  std::vector<double> x(2);
				// Guesses
  x[0] = 0.5*(Emax-Emin);
  x[1] = 0.5*(Kmax-Kmin);

  // SimAnn takes a std::function object as its first argument
  //
  SimAnn sa (*this, 2);

  sa.initial(x);          // set initial condition
  sa.learning_rate(RATE);

  double t0;

  if (MELT) {
    sa.melt();			// melt the system
				// make it a bit warmer than melting temperature
    t0 = 1.2 * sa.temperature();
  } else {
    t0 = T0;
  }

  sa.temperature(t0);
        
  double t = sa.anneal(MAXIT);
  x = sa.optimum();

  double ee, kk;
  mapvars(x, ee, kk);
  orb->new_orbit(ee, kk);

  OrbValues ret;

  ret.Boltzmann = sa.Boltzmann();
  ret.rate      = sa.learning_rate();
  ret.t0        = t0;
  ret.tf        = t;
  ret.energy    = ee;
  ret.kappa     = kk;
  ret.value     = (*this)(x);
  ret.peri      = orb->peri();
  ret.apo       = orb->apo();
  ret.radial_period  = 2.0*M_PI/orb->get_freq(0);
  ret.azimuthal_period  = 2.0*M_PI/orb->get_freq(1);

  return ret;
}
