#ifndef _DISK_MODELS_H
#define _DISK_MODELS_H

#include <EmpCylSL.H>


//! A Ferrers Ellipsoid + Evacuated Exponential Disc (semi-realistic bar+disk model)
class FEED : public EmpCylSL::AxiDisk
{
private:
  double a, b, c, r0, z0, fbar;
  /*
  a    : the major axis length for the Ferrers bar
  b    : the intermediate axis length for Ferrers bar
  c    : the minor axis length for Ferrers bar
  r0   : the scale length for the exponential disc
  z0   : the scale height for the exponential disc
  fbar : the fraction of the total mass in the Ferrers component

  CAUTION: some parameter combinations result in aphysical densities, 
  which this routine will simply zero out. 

  This almost certainly is not the desired behavior; it is recommended that you check the density structure of the input model
  before proceeding.

   */

public:

  FEED(double a, double b, double c, double r0, double z0, double fbar, double M=1) :
    a(a), b(b), c(c), r0(r0), z0(z0), fbar(fbar), EmpCylSL::AxiDisk(M, "FEED") {}

  double operator()(double R, double z, double phi=0.)
  {
     double x = R*cos(phi);
     double y = R*sin(phi);
     double sh = 1.0/cosh(z/z0);

     double fm, bardens, ferrersdisk, basedisk;

     // set up the Ferrers component
     fm = sqrt( (x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c) );

     // only allowing for n=2 at present
     double n = 2;
     double rho0 = fbar*(105/(32*M_PI))/(a*b*c);

     if (fm>1) {
       bardens = 0.;
     } else {
       bardens = rho0 * pow((1-fm*fm),n);
     }
     
     // set up the (to be subtracted) Ferrers exponential disk
     double fdiscr0 = a/3.; // scaled to the Ferrers bar major axis length
     ferrersdisk = fbar*0.25/(M_PI*fdiscr0*fdiscr0*z0) * exp(-R/fdiscr0) * sh * sh;
     
     // set up the base exponential disk
     basedisk = 0.25/(M_PI*r0*r0*z0) * exp(-R/r0) * sh * sh;

     return M * (basedisk - ferrersdisk + bardens);
     
  }

};



//! An n=2 Ferrers Ellipsoid
class Ferrers : public EmpCylSL::AxiDisk
{
private:
  double a, b, c;

public:

  Ferrers(double a, double b, double c, double M=1) :
    a(a), b(b), c(c), EmpCylSL::AxiDisk(M, "Ferrers") {}

  double operator()(double R, double z, double phi=0.)
  {
     double x = R*cos(phi);
     double y = R*sin(phi);
     double fm;

     fm = sqrt( (x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c) );

     // only allowing for n=2 at present
     double n = 2;
     double rho0 = (105/(32*M_PI))/(a*b*c);

     if (fm>1) {
       return 0.;
     } else {
       return rho0 * pow((1-fm*fm),n);
     }
  }

};




//! The usual exponential disk
class Exponential : public EmpCylSL::AxiDisk
{
private:
  double a, h;

public:

  Exponential(double a, double h, double M=1) :
    a(a), h(h), EmpCylSL::AxiDisk(M, "Exponential") {}

  double operator()(double R, double z, double phi=0.)
  {
    double sh = 1.0/cosh(z/h);
    return 0.25*M/(M_PI*a*a*h) * exp(-R/a) * sh * sh;
  }

};




//! Miyamoto-Nagai disk
class MNdisk : public EmpCylSL::AxiDisk
{
private:
  double a, h;

public:

  MNdisk(double a, double h, double M=1) :
    a(a), h(h), EmpCylSL::AxiDisk(M, "Miyamoto-Nagai") {}

  double operator()(double R, double z, double phi=0.)
  {
    double Z2 = z*z + h*h;
    double Z  = sqrt(Z2);
    double Q2 = (a + Z)*(a + Z);
    return 0.25*h*h*M/M_PI*(a*R*R + (a + 3.0*Z)*Q2)/( pow(R*R + Q2, 2.5) * Z*Z2 );
  }
};

//! Truncate a AxiDisk
class Truncated : public EmpCylSL::AxiDisk
{
  typedef boost::shared_ptr<EmpCylSL::AxiDisk> AxiDiskPtr;
  
private:
  
  AxiDiskPtr model;
  double rt, rw;

public:

  Truncated(double R, double W, AxiDiskPtr p) :
    rt(R), rw(W), model(p)
  {
    ID = p->getID() + " [truncated]";
  }

  double operator()(double R, double z, double phi=0.)
  {
    return (*model)(R, z, phi) * erf( (rt - sqrt(R*R + z*z))/rw );
  }

};

#endif

