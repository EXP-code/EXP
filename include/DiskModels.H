#ifndef _DISK_MODELS_H
#define _DISK_MODELS_H

#include <EmpCylSL.H>
#include <DiskDensityFunc.H>

//! A Ferrers Ellipsoid + Evacuated Exponential Disc (semi-realistic
//! bar+disk model)
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

  This almost certainly is not the desired behavior; it is recommended
  that you check the density structure of the input model before
  proceeding.

   */

public:

  FEED(double a, double b, double c, double r0, double z0, double fbar, double M=1) :
    a(a), b(b), c(c), r0(r0), z0(z0), fbar(fbar), EmpCylSL::AxiDisk(M, "FEED")
  {
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);
    params.push_back(r0);
    params.push_back(z0);
    params.push_back(fbar);
  }

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
    a(a), b(b), c(c), EmpCylSL::AxiDisk(M, "Ferrers")
  {
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);
  }

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
    a(a), h(h), EmpCylSL::AxiDisk(M, "Exponential")
  {
    params.push_back(a);
    params.push_back(h);
  }

  double operator()(double R, double z, double phi=0.)
  {
    double sh = 1.0/cosh(z/h);
    return 0.25*M/(M_PI*a*a*h) * exp(-R/a) * sh * sh;
  }

};

//! A mixture of two exponential disks
class DoubleExponential : public EmpCylSL::AxiDisk
{
private:
  double a, h, a2, h2, Mfac, d1, d2;

public:

  DoubleExponential(double a, double h, double a2, double h2, double Mfac, double M=1) :
    a(a), h(h), a2(a2), h2(h2), Mfac(Mfac), EmpCylSL::AxiDisk(M, "DoubleExponential")
  {
    params.push_back(a);
    params.push_back(h);
    params.push_back(a2);
    params.push_back(h2);
    params.push_back(Mfac);
  }

  double operator()(double R, double z, double phi=0.)
  {
    double sh = 1.0/cosh(z/h);
    d1 = 0.25*M*Mfac/(M_PI*a*a*h) * exp(-R/a) * sh * sh;
    
    sh = 1.0/cosh(z/h2);
    d2 = 0.25*M*(1-Mfac)/(M_PI*a2*a2*h2) * exp(-R/a2) * sh * sh;

    return d1 + d2;
  }

};



//! Miyamoto-Nagai disk
class MNdisk : public EmpCylSL::AxiDisk
{
private:
  double a, h;

public:

  MNdisk(double a, double h, double M=1) :
    a(a), h(h), EmpCylSL::AxiDisk(M, "Miyamoto-Nagai")
  {
    params.push_back(a);
    params.push_back(h);
  }

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
  typedef std::shared_ptr<EmpCylSL::AxiDisk> AxiDiskPtr;
  
private:
  
  AxiDiskPtr model;
  double rt, rw;

public:

  Truncated(double R, double W, AxiDiskPtr p) :
    rt(R), rw(W), model(p)
  {
    ID = p->getID() + " [truncated]";
    params.push_back(rt);
    params.push_back(rw);
  }

  double operator()(double R, double z, double phi=0.)
  {
    return (*model)(R, z, phi) * erf( (rt - sqrt(R*R + z*z))/rw );
  }

};


//! The usual exponential disk + a Hernquist bulge
class diskbulge : public EmpCylSL::AxiDisk
{
private:
  double a, h, as, Mfac, d1, d2;

public:

  diskbulge(double a, double h, double as, double Mfac, double M=1) :
    a(a), h(h), as(as), Mfac(Mfac), EmpCylSL::AxiDisk(M, "diskbulge")
  {
    params.push_back(a);
    params.push_back(h);
    params.push_back(as);
    params.push_back(Mfac);
  }

  double operator()(double R, double z, double phi=0.)
  {
    double sh = 1.0/cosh(z/h);
    d1 = 0.25*M*Mfac/(M_PI*a*a*h) * exp(-R/a) * sh * sh;

    double rr = pow(pow(R, 2) + pow(z,2), 0.5);
    d2 = M*(1-Mfac)*pow(as, 4)/(2.0*M_PI*rr)*pow(rr+as,-3.0);

    return d1 + d2;

  }

};


//! A user-defined Python model
class PyModel : public EmpCylSL::AxiDisk
{
private:

  std::shared_ptr<DiskDensityFunc> pyDens;

public:

  PyModel(std::string& pyname)
  {
    pyDens = std::make_shared<DiskDensityFunc>(pyname);
  }

  double operator()(double R, double z, double phi=0.)
  {
    return (*pyDens)(R, z, phi);
  }

};

//! Enforce an axsisymmetric model by using the x-axis only from a
//! user-defined Python model
class AxiSymPyModel : public EmpCylSL::AxiDisk
{
private:

  std::shared_ptr<DiskDensityFunc> pyDens;

public:

  AxiSymPyModel(std::string& pyname)
  {
    pyDens = std::make_shared<DiskDensityFunc>(pyname);
  }

  double operator()(double R, double z, double phi=0.)
  {
    return (*pyDens)(R, z, 0.0);
  }

};

#endif
