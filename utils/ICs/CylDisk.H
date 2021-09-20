// This is -*- C++ -*-

//
// Two dimensional disk from Cylindrical expansion
//

#ifndef _CylDisk_h
#define _CylDisk_h

#include "massmodel.H"
#include "EmpCylSL.H"
#include "interp.H"

class CylDisk : public AxiSymModel
{
private:
  
  EmpCylSL *d;
  void setup_model(double rmin, double rmax);

  bool rlog;
  vector<double> rr, dd, mm, pp;
  double dr, rmin, rmax;

public:

				// Number of radial points in grid
  static int NUMR;		// Default: 200
				// Linear (false) or logarithmic (true) mesh
  static bool RLOG;		// Default: true

  CylDisk(EmpCylSL *D, double Rmin, double Rmax) 
  {
    d       = D;
    dim     = 2;

    ModelID = "CylDisk";
    dist_defined = false;

    setup_model(Rmin, Rmax);
  }
      

  // Required member functions

  double get_mass(const double r) 
  {
    double R = r;
    if (rlog) R = log(r);
    return odd2(R, rr, mm);
  }

  double get_density(const double r) 
  {
    double R = r;
    if (rlog) R = log(r);
    return odd2(R, rr, dd);
  }

  double get_pot(const double r) 
  {
    double R = r;
    if (rlog) R = log(r);
    return odd2(R, rr, pp);
  }
							

  double get_dpot(const double r) 
  {
    double R = r;
    if (rlog) R = log(r);
    return drv2(R, rr, pp);
  }

  double get_dpot2(const double r) 
  {
    double R = r;
    if (rlog) R = log(r);
    int indx = Vlocate(R, rr);
    if (indx >= NUMR-2) indx = NUMR-4;
    else if (indx>0) indx--;
    return ( (pp[indx+3] - 3.0*pp[indx+2]+ 3.0*pp[indx+1] - pp[indx])*(R - rr[indx])
	     - (pp[indx+3] - 4.0*pp[indx+2]+ 5.0*pp[indx+1] - 2.0*pp[indx])*dr ) /
      (dr*dr*dr);
  }

  void get_pot_dpot(const double r, double &ur, double &dur) 
  {
    double R = r;
    if (rlog) R = log(r);
    ur  = drv2(R, rr, pp);

    int indx = Vlocate(R, rr);
    if (indx >= NUMR-2) indx = NUMR-4;
    else if (indx>0) indx--;
    dur = ( (pp[indx+3] - 3.0*pp[indx+2]+ 3.0*pp[indx+1] - pp[indx])*(R - rr[indx])
	    - (pp[indx+3] - 4.0*pp[indx+2]+ 5.0*pp[indx+1] - 2.0*pp[indx])*dr ) /
      (dr*dr*dr);
  }
  
  // Addiional member functions

  double get_min_radius(void) 
  { 
    return rmin;
  }
  double get_max_radius(void) 
  {
    return rmax;
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
