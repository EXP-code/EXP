#ifndef _FIND_ORB_H
#define _FIND_ORB_H

#include <Func1d.H>

struct OrbValues {
  double Boltzmann, rate, t0, tf;
  double energy, kappa, value;
  double peri, apo;
  double radial_period, azimuthal_period;
};

class FindOrb : public Func1d
{
private:
  SphericalOrbit* orb;
  AxiSymModel *halo_model;

  double Kmin, Kmax, Emin, Emax;
  double apo, peri;
  double rate, melt;
  int maxit;

  void mapvars(double *x, double& ee, double& kk);
  double CostFunction(double* ek);

public:

  /// Defaults
  //@{
  static double KMIN;
  static double KMAX;
  static int    MAXIT;
  static bool   MELT;
  static double RATE;
  static double EFAC;
  static double T0;
  //@}
  
  /// Constructor
  FindOrb(AxiSymModel *mod, double PERI, double APO);

  /// Destructor
  ~FindOrb();

  /// Do the work!
  OrbValues Anneal();

  /// Return the orbit instance (presumably now at minimum)
  SphericalOrbit& Orb() { return *orb; }
};

#endif
