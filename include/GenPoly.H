#ifndef _GenPoly_H
#define _GenPoly_H

#include <cmath>
#include "massmodel.H"

class GeneralizedPolytrope : public AxiSymModel
{
private:
  RUN mass;
  RUN dens;
  RUN pot;

  double n, m, KF;

public:

  GeneralizedPolytrope(void);
  GeneralizedPolytrope(int num, double n, double m,
		       double eps0=1.0e-5, double step=1.0e-5);

  // Required member functions

  double get_mass(const double);
  double get_density(const double);
  double get_pot(const double);
  double get_dpot(const double);
  double get_dpot2(const double);
  void get_pot_dpot(const double, double&, double&);
  
  // Additional member functions

  double get_min_radius(void) { return mass.x[1]; }
  double get_max_radius(void) { return mass.x[mass.num]; }

  double distf(double E, double L);
  double dfde(double E, double L);
  double dfdl(double E, double L);
  double d2fde2(double E, double L);

  // Dump model grid
  void dump_model(ostream* out);
};

#endif
