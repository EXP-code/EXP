#ifndef _UserMW_H
#define _UserMW_H

#include "Particle.H"
#include "AxisymmetricBasis.H"
#include "ExternalCollection.H"

/** Potential and force of the Milky Way, a la Gala.

    !!!ADD ALL PARAMETERS YOU WANT TO DEFINE, using this format (for automatic documentation)!!! Done
    @param G is the gravitational constant

    @param M_halo is the total mass of the Navarro-Frank-White (NFW) halo component
    @param rs_halo is the scale radius of the Navarro-Frank-White (NFW) halo component

    @param M_disk is the total mass of the Miyamoto-Nagai (MN) disk component
    @param a_disk is the scale length of the Miyamoto-Nagai (MN) disk component
    @param b_disk is the scale height of the Miyamoto-Nagai (MN) disk component

    @param M_nucl is the total mass of the Hernquist (HN) nucleus component
    @param c_nucl is the concentration of the Hernquist (HN) nucleus component

    @param M_bulge is the total mass of the Hernquist (HN) bulge component
    @param c_bulge is the concentration of the Hernquist (HN) bulge component

    @param Ton is the time at the center of the <code>erf</code> turn on
    @param Toff is the time at the center of the <code>erf</code> turn off
    @param DeltaT is the spread of the turn on
    @param ctrname defines the component that defines the center (otherwise the center is fixed to the origin)
*/
class UserMW : public ExternalForce
{
private:

  string ctr_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();


  // !!!define all variable names needed here with their types !!! Done
  double G, M_halo, rs_halo, M_disk, a_disk, b_disk, M_nucl, c_nucl, M_bulge, c_bulge;


  double Ton, Toff, DeltaT;

  string name;

  void userinfo();

  // !!! define necessary force definitions here: they will be accessible in UserMW.cc in this case !!! Done

  double NFW_pot (double r)
  {
      double phi;
      phi =  -(G*M_halo)/r * log( 1 + (r/rs_halo) );
      return phi;
  }

  void NFW_dphi_dr (double x, double y, double z, double *ax1, double *ay1, double *az1) // Work from here
  {
      double dphi_dr, r;

      r = sqrt( pow(x,2) + pow(y,2) );
      dphi_dr = (G*M_halo)/r * (log( 1 + (r/rs_halo) )/r - (1/(rs_halo + r)));

      *ax1 = -dphi_dr * (x/r);
      *ay1 = -dphi_dr * (y/r);
      *az1 = -dphi_dr * (z/r);
      // return dphi_dr, ax, ay, az;
  }

  double HN_nucl_pot (double r)
  {
      double phi;
      phi = -G*M_nucl/(r + c_nucl);
      return phi;
  }

  void HN_nucl_dphi_dr (double x, double y, double z, double *ax2, double *ay2, double *az2)
  {
      double dphi_dr, r;

      r = sqrt( pow(x,2) + pow(y,2) );
      dphi_dr = G*M_nucl/pow ( (r + c_nucl), 2 );

      *ax2 = -dphi_dr * (x/r);
      *ay2 = -dphi_dr * (y/r);
      *az2 = -dphi_dr * (z/r);
      // return dphi_dr, ax, ay, az;
  }

  double HN_bulge_pot (double r)
  {
      double phi;
      phi = -G*M_bulge/(r + c_bulge);
      return phi;
  }

  void HN_bulge_dphi_dr (double x, double y, double z, double *ax3, double *ay3, double *az3)
  {
      double dphi_dr, r;

      r = sqrt( pow(x,2) + pow(y,2) );
      dphi_dr = G*M_bulge/pow ( (r + c_bulge), 2 );

      *ax3 = -dphi_dr * (x/r);
      *ay3 = -dphi_dr * (y/r);
      *az3 = -dphi_dr * (z/r);
      // return dphi_dr, ax, ay, az;
  }

  double MN_pot (double z, double R)
  {
      double phi;
      phi = (-G*M_disk)/(sqrt(pow( R, 2 ) + pow( (sqrt(pow( z, 2)  + pow( b_disk, 2 )) + a_disk), 2 )));
      return phi;
  }

 void MN_dphi_dR_dz (double x, double y, double z, double *ax4, double *ay4, double *az4)
  {
      double dphi_dR, dphi_dz, R;

      R = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
      double n = 3.0/2.0;

      double Rstep3 = sqrt( pow(b_disk,2) + pow(z,2) );
      double Rstep2 = (a_disk + Rstep3);
      double Rstep1 = pow(R,2) + pow(Rstep2, 2);
      dphi_dR = (G*M_disk*R)/pow(Rstep1, n);

      double zstep3 = sqrt( pow(b_disk,2) + pow(z,2) );
      double zstep2 = (a_disk + zstep3);
      double zstep1 = pow(R,2) + pow(zstep2, 2);
      dphi_dz = ((G*M_disk*z)*zstep2)/(zstep3*pow(zstep1, n));

      *ax4 = -dphi_dR * (x/R);
      *ay4 = -dphi_dR * (y/R);
      *az4 = -dphi_dz;
      // return dphi_dR, dphi_dz, ax, ay, az;
  }

  // double MN_dphi_dz (double z, double R)
  // {
  //     double dphi_dz;
  //     double n = 3.0/2.0;
  //     double step3 = sqrt( pow(b_disk,2) + pow(z,2) );
  //     double step2 = (a_disk + step3);
  //     double step1 = pow(R,2) + pow(step2, 2);
  //     dphi_dz = ((G*M_disk*z)*step2)/(step3*pow(step1, n));
  //     return dphi_dz;
  // }

public:

  //! Constructor
  UserMW(const YAML::Node& conf);

  //! Destructor
  ~UserMW();

};

#endif
