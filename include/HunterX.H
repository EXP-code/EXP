//
// Hunter's Legendre Models with Maclauring subtracted
// as defined in Hunter & Toomre
//

#ifndef _HunterX.H
#define _HunterX.H

#include "massmodel.H"
#include <hunter.h>

class HunterDiskX : public AxiSymModel
{
private:
  
  HunterDisk *diskN, *disk1;
  double rmax;
  double factor;

public:

				// Constructor

  HunterDiskX(int N, double RMAX=1.0, double MASS=1.0);

				// Destructor
  ~HunterDiskX();
      
  void set_params(const double KMAX=20.0, const double RCUT=8.0, 
		  const int NINT=200, const int MGRID=200, 
		  const int NHANK=400);


  // Required member functions

  double get_mass(const double r);

  double get_density(const double r);

  double get_pot(const double r);

  double get_dpot(const double r);

  double get_dpot2(const double r);

  void get_pot_dpot(const double r, double &ur, double &dur);
  
  // Addiional member functions

  double get_min_radius(void) { return 0.0; }
  // double get_max_radius(void) { return rmax; }
  double get_max_radius(void) { return 1e10; }

  double distf(double E, double L) {
    bomb("Function <distf> not implemented");
  }

  double dfde(double E, double L) {
    bomb("Function <dfde> not implemented");
  }

  double dfdl(double E, double L) {
    bomb("Function <dfdl> not implemented");
  }

  double d2fde2(double E, double L) {
    bomb("Function <d2fde2> not implemented");
  }

  double get_pot(const double r, const double z);

  double get_dpot(const double r, const double z);

  double get_dpot2(const double r, const double z);

  void get_pot_dpot(const double r, const double z, 
		    double &ur, double &dur);

};



#endif
