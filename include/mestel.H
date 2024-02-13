#ifndef _mestel_H
#define _mestel_H

#include <massmodel.H>

//! Infiinte Mestel Disk
class MestelDisk : public AxiSymModel
{
protected:
  
  double vrot;
  double rot;
  double rmin;
  double rmax;
  double sig2;
  double q;
  double F;

public:

  //! Constructor
  MestelDisk(double VROT = 1.0, double RMIN = 1.0e-6, double RMAX = 1.0e6);

  //@{
  //! Required member functions
  virtual double get_mass(const double R);

  virtual double get_density(const double R);

  double get_pot(const double R);

  double get_dpot(const double R);

  double get_dpot2(const double R);
  
  void get_pot_dpot(const double R, double &ur, double &dur);
  //@}
  
  //@{
  //! Addiional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }

  virtual void setup_df(double sigma);

  virtual double distf(double E, double L);

  virtual double dfde(double E, double L);

  virtual double dfdl(double E, double L);

  virtual double d2fde2(double E, double L);
  //@}

};


//! Doubly tapered Mestel Disk
class TaperedMestelDisk : public MestelDisk
{
private:

  std::shared_ptr<Linear1d> interp;

protected:
  //! Taper exponents
  double nu, mu;

  //! Inner taper radius
  double Ri;

  //! Outer taper radius
  double Ro;

  //! Taper factors
  double Tifac, Tofac;

  //! Inner taper function
  double Tinner(double Jp);

  //! Outer taper function
  double Touter(double Jp);

  //! Deriv of inner taper function
  double dTinner(double Jp);

  //! Deriv of outer taper function
  double dTouter(double Jp);

public:

  //! Constructor
  TaperedMestelDisk(double nu, double mu, double Ri, double Ro,
		    double VROT = 1.0, 
		    double RMIN = 1.0e-6, double RMAX = 1.0e6);
      
  //! Overloaded density function
  double get_density(const double R);

  //! Overloaded cumulative mass function
  double get_mass(double R);

  //@{
  //! Overloaded DF member functions

  double distf(double E, double L);

  double dfde(double E, double L);

  double dfdl(double E, double L);
  
  double d2fde2(double E, double L);
  //@}

};

#endif
