//
// Hunter's Legendre Models
// 1964, MN, 129, 321
//

#ifndef _hunter_h
#define _hunter_h

#include "massmodel.H"

class HunterDisk : public AxiSymModel
{
private:
  
  int n;
  double rmax, mass, rcut;
  double sfac, omfac, pfac, kfac;

  Vector gam, cvec;

  Logic model_setup;
  void setup_model(void);
  int nint;

  Vector vplgndr(double), dvplgndr(double), d2vplgndr(double);
  
  Logic mass_setup;
  void setup_mass(void);
  int mgrid;
  Vector r_mass, m_mass;

  Logic hankel_setup;
  Vector S_K, K;
  double kmax;
  int  nhank;

  void tabulate_hankel(void);

  double offset_p;
  double offset_dp;
  double offset_dp2;

  const static int defNINT=200;
  const static int defMGRID=200;
  const static int defNHANK=400;
  const static double defKMAX=20.0;
  const static double defRCUT=8.0;

public:

				// Constructors

  HunterDisk(int N, double RMAX=1.0, double MASS=1.0);
      
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
  double get_max_radius(void) { return 1.0e10; }

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
