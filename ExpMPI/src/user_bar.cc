
/*
  void user_perturbation()

  Compute acceleration and potential for some desired
  disturbance.

  Turned on by setting flag "user".

*/

static char rcsid[] = "$Id";

#include <math.h>
#include "expand.h"
#undef RNUM			// Conflict
#include <SatelliteOrbit.h>

/*  U1 (0.5)	Bar length		*/
/*  U2 (0.3)	Bar amplitude		*/
/*  U3 (-20.0)	Turn on start time	*/
/*  U4 (1.0)	Turn on duration	*/
/*  U5 (0.0)	Corotation factor	*/

class UserID {
public:
  UserID() {
    cout << "\n";
    cout << "User routine id:\n";
    cout << "Bar perturbation\n";
    cout << "\n";
  }
} blab;

				/* Function defs */

extern "C" {

  void determine_fields_at_point_bes(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);

  void determine_fields_at_point_CB(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);

  void determine_fields_at_point_CBDisk(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);

  void determine_fields_at_point_HERNQ(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);

  void determine_fields_at_point_Cyl(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);

  void determine_fields_at_point_SLsphere(double r, double theta, double phi,
   double *dens, double *potl, double *potr, double *pott, double *potp);
}

extern "C" void user_perturbation()
{
  static const double numfac = 3.86274202023190e-01;
  static bool firstime = true;
  static double posang, lastomega, lasttime;

  int i;
   

				// Get frequency
  double R=U1*U5, theta=1.57079632679490, phi;
  double dens, potl, potr, pott, potp;
  double avg = 0.0;

  for (int n=0; n<8; n++) {
    phi = 2.0*M_PI/8.0 * n;

    if (bessel_sph)
      determine_fields_at_point_bes(R, theta, phi, 
                                    &dens, &potl, &potr, &pott, &potp);
    else if (c_brock)
      determine_fields_at_point_CB(R, theta, phi, 
                                    &dens, &potl, &potr, &pott, &potp);
    else if (c_brock_disk)
      determine_fields_at_point_CBDisk(R, theta, phi, 
                                    &dens, &potl, &potr, &pott, &potp);
    else if (hernq)
      determine_fields_at_point_HERNQ(R, theta, phi, 
                                    &dens, &potl, &potr, &pott, &potp);
    else if (sphereSL)
      determine_fields_at_point_SLsphere(R, theta, phi, 
				      &dens, &potl, &potr, &pott, &potp);
    else
      dens = potl = potr = pott = potp = 0.0;

    avg += potr/8.0;
  }

  double omega = sqrt(avg/R);

  if (firstime) {
    posang = 0.0;
    lastomega = omega;
    lasttime = tnow;
    firstime = false;
  } else {
    // posang += 0.5*(omega + lastomega)*(tnow - lasttime);
    posang += 0.5*(omega + lastomega)*dtime;
    lastomega = omega;
    lasttime = tnow;
#ifdef DEBUG
    if (myid==0)
      cout << "Time=" << tnow << " Posang=" << posang << endl;
#endif
  }

  double fac, ffac, amp = U2 * 0.5*(1.0 + erf( (tnow - U3)/U4 ));
  double xx, yy, zz, rr, nn; 
  double cos2p = cos(2.0*posang);
  double sin2p = sin(2.0*posang);

  for (i=1; i<=nbodies; i++)
    {
      if (freeze_particle(i)) continue;

      xx = x[i]-com1[0];
      yy = y[i]-com1[1];
      zz = z[i]-com1[2];
      rr = sqrt( xx*xx + yy*yy + zz*zz );

      fac = U1 + rr;
		     
      ffac = -amp*numfac*U1*U1*U1/pow(fac, 6.0);

      nn = (xx*xx - yy*yy)*cos2p + 2.0*xx*yy*sin2p;

      ax[i] += ffac*
	(  2.0*xx*cos2p*fac + 2.0*yy*sin2p*fac - 5.0*nn*xx/(rr+1.0e-16) );

      ay[i] += ffac*
	( -2.0*yy*cos2p*fac + 2.0*xx*sin2p*fac - 5.0*nn*yy/(rr+1.0e-16) );

      az[i] += ffac*
	( -5.0*nn*zz/(rr+1.0e-16) );

      potext[i] += -ffac*nn*fac;

    }
}
