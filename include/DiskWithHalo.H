//
// Two dimensional disk with surrounding halo
//

#ifndef _DiskWithHalo_H
#define _DiskWithHalo_H

#include <massmodel.H>

class DiskWithHalo : public AxiSymModel
{
private:
  
  std::shared_ptr<AxiSymModel> d;
  std::shared_ptr<AxiSymModel> h;

public:

  // Null constructor
  DiskWithHalo() : d(0), h(0) {}

  // Useful constructor
  DiskWithHalo(std::shared_ptr<AxiSymModel> D, std::shared_ptr<AxiSymModel> H) 
  {
    d       = D;
    h       = H;
    dim     = 2;
    ModelID = "DiskAndHalo [" + D->ModelID + " + " + H->ModelID + "]";
    dist_defined = false;
  }
      
  // Copy constructor
  DiskWithHalo(const DiskWithHalo& t) 
  {
    d       = t.d;
    h       = t.h;
    dim     = 2;
    ModelID = t.ModelID;
    dist_defined = false;
  }

  // Required member functions

  virtual double get_mass(const double r) override
  { return d->get_mass(r); }

  virtual double get_density(const double r) override
  { return d->get_density(r); }

  virtual double get_pot(const double r) override
  { return d->get_pot(r) + h->get_pot(r); }
							

  virtual double get_dpot(const double r) override
  { return d->get_dpot(r) + h->get_dpot(r); }

  virtual double get_dpot2(const double r) override
  { return d->get_dpot2(r) + h->get_dpot2(r); }

  virtual void get_pot_dpot(const double r, double &ur, double &dur) override
  {
    double ur1, dur1, ur2, dur2;
    d->get_pot_dpot(r, ur1, dur1);
    h->get_pot_dpot(r, ur2, dur2);
    ur  = ur1  + ur2;
    dur = dur1 + dur2;
  }
  
  double get_mass(const double x1, const double x2, const double x3)
  { return get_mass(sqrt(x1*x1 + x2*x2 + x3*x3)); }
  
  double get_density(const double x1, const double x2, const double x3)
  { return get_density(sqrt(x1*x1 + x2*x2 + x3*x3)); }
  
  double get_pot(const double x1, const double x2, const double x3)
  { return get_pot(sqrt(x1*x1 + x2*x2 + x3*x3)); }

  // Addiional member functions

  double get_min_radius(void) override
  { 
    return max<double>(d->get_min_radius(), h->get_min_radius());
  }
  double get_max_radius(void) override
  {
    return min<double>(d->get_max_radius(), h->get_max_radius());
  }

  double distf(double E, double L) override {
    bomb("Dist fct not defined!");
    return 0.0;
  }

  double dfde(double E, double L) override {
    bomb("Dist fct not defined!");
    return 0.0;
  }

  double dfdl(double E, double L) override {
    bomb("Dist fct not defined!");
    return 0.0;
  }
  
  double d2fde2(double E, double L) override {
    bomb("Dist fct not defined!");
    return 0.0;
  }

};

#endif
