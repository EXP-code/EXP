// This may look like C code, but it is really -*- C++ -*-

//
// Two dimensional disk with surrounding halo
//

#ifndef _DiskWithHalo_h
#define _DiskWithHalo_h

#include "massmodel.h"

class DiskWithHalo : public AxiSymModel
{
private:
  
  AxiSymModel *d;
  AxiSymModel *h;

public:

  DiskWithHalo(AxiSymModel *D, AxiSymModel* H) {
    d       = D;
    h       = H;
    dim     = 2;
    ModelID = "DiskAndHalo [" + D->ModelID + " + " + H->ModelID + "]";
    dist_defined = false;
  }
      

  // Required member functions

  double get_mass(const double r) { return d->get_mass(r); }

  double get_density(const double r) { return d->get_density(r); }

  double get_pot(const double r) {
    return d->get_pot(r) + h->get_pot(r);
  }
							

  double get_dpot(const double r) {
    return d->get_dpot(r) + h->get_dpot(r);
  }

  double get_dpot2(const double r) {
    return d->get_dpot2(r) + h->get_dpot2(r);
  }

  void get_pot_dpot(const double r, double &ur, double &dur) {
    double ur1, dur1, ur2, dur2;
    d->get_pot_dpot(r, ur1, dur1);
    h->get_pot_dpot(r, ur2, dur2);
    ur  = ur1  + ur2;
    dur = dur1 + dur2;
  }
  
  // Addiional member functions

  double get_min_radius(void) 
  { 
    return max<double>(d->get_min_radius(), h->get_min_radius());
  }
  double get_max_radius(void) 
  {
    return min<double>(d->get_max_radius(), h->get_max_radius());
  }

  double distf(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }

  double dfde(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }

  double dfdl(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }
  
  double d2fde2(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }

};

#endif
