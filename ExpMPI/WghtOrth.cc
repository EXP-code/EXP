#include <iostream.h>
#include <iomanip.h>
#include <String.h>

#include <numerical.h>
#include <gaussQ.h>
#include "WghtOrth.h"

#include "expand.h"


#undef TINY
#define TINY 1.0e-16


CylindricalCB::CylindricalCB(void)
{
  H = 1.0;
  M = 0;
  MAXR = 0;
  NFFT = 0;
  NINT = defNINT;
  ZMAX = defZMAX;
  initialized = FALSE;
  use_tables = FALSE;
  MPIset = FALSE;
  coefs_made = FALSE;
}

CylindricalCB::~CylindricalCB(void)
{
  if (initialized) {
    delete [] ortho;
    delete [] U;
    delete [] UT;
    delete fft;
  }

  if (MPIset) {
    delete [] MPIin;
    delete [] MPIout;
  }
}

CylindricalCB::CylindricalCB(int mu, int nfft, int m, double a, double h) 
{
  MAXR = mu;
  NFFT = nfft;
  M = m;
  A = a;
  H = h;

  NINT = defNINT;
  ZMAX = defZMAX;
  initialized = FALSE;
  use_tables = FALSE;
  MPIset = FALSE;
  coefs_made = FALSE;
}

void CylindricalCB::reset(int mu, int nfft, int m, double a, double h)
{
  MAXR = mu;
  NFFT = nfft;
  M = m;
  A = a;
  H = h;

  if (initialized) {
    delete [] ortho;
    delete [] U;
    delete [] UT;
    delete fft;
  }

  initialize();
}


void CylindricalCB::initialize(void)
{
  if (M<0)  bomb ("M must be >= 0");
  if (MAXR<0)  bomb ("MU must be >= 0");
  if (NFFT<=0)  bomb ("NFFT must be > 0");

  int i, j, k;

  //======================================================================
  //
  // Step One: setup Fourier transform
  //
  //======================================================================

  double pos, neg, z, offset;

  n = 2<<(NFFT-1);
  n2 = n/2;

  dZ = ZMAX/n;
  dk = 2*M_PI/n/dZ;

  output.setsize(0, MAXR, 0, n-1);
  Z.setsize(0, n2-1, 0, MAXR);

  fft = new FFT (n, dZ);


  //======================================================================
  //
  // Step Two: compute normalization integrals
  //
  //======================================================================

  Matrix ortho_offset(1, MAXR+1, 1, MAXR+1);
  Matrix ev_offset(1, MAXR+1, 1, MAXR+1);
  Vector z_offset(1, MAXR+1);
  
  ortho = new Matrix [n2];
  U = new Matrix [n2];
  UT = new Matrix [n2];

  double A2 = A*A, M1 = 1.0+M;
  double A4 = A2*A2; 
  double faci, facj, factor, r, r2, t, k2, wght;

  LegeQuad lwk(NINT);

  for (int iz = 0; iz<n2; iz++) {

    k2 = dk*dk*iz*iz;

    ortho[iz].setsize(0, MAXR, 0, MAXR);
    ortho[iz].zero();
    
    for (k=1; k<=NINT; k++) {
    
      t = lwk.knot(k);
      r = A*t/sqrt(1.0 - t*t);
      r2 = r*r;

      factor = (A2 + r2)*(A2 + r2)*(k2 + M*M/(r2+TINY));

      wght = A*pow(1.0 - t*t, -1.5) * r *
	pow(A2 + r2, -3.0 - M) * lwk.weight(k);

      faci = 1.0;

      for (i=0; i<=MAXR; i++) {	// Potential loop

	facj = 1.0;

	for (j=0; j<=i; j++) {	// Density loop
	  
	  ortho[iz][i][j] += wght * faci * facj * 
	    (A4*j*j/(r2+TINY) - 2.0*((1.0 + j)*M1 + j)*A2  + M1*M1*r2
	      - factor);

	  facj *= t;
	}

	faci *= t;
      }
    }


    for (i=0; i<=MAXR; i++) {
      for (j=i+1; j<=MAXR; j++) {
	ortho[iz][i][j] = ortho[iz][j][i];
      }
    }

    ortho[iz] *= -1.0;


  //======================================================================
  //
  // Step Five: find orthogonal transform for radial component
  //
  //======================================================================


    U[iz].setsize(0, MAXR, 0, MAXR);

    for (i=0; i<=MAXR; i++) {
      for (j=0; j<=MAXR; j++) {
	ortho_offset[i+1][j+1] = ortho[iz][i][j];
      }
    }

    z_offset = ortho_offset.Symmetric_Eigenvalues_GHQL(ev_offset);

    for (i=0; i<=MAXR; i++) {
      for (j=0; j<=MAXR; j++) {
	U[iz][i][j] = ev_offset[i+1][j+1];
      }
      Z[iz][i] = 1.0/sqrt(fabs(z_offset[i+1]));
    }
    UT[iz] = U[iz].Transpose();

    /*

    cout << endl;
    for (i=0; i<=MAXR; i++)
      cout << setw(5) << i << setw(5) << iz
	   << setw(15) << Z[iz][i] << endl;
    cout << endl;

    Matrix check = UT[iz] * U[iz];
    for (i=0; i<=MAXR; i++) {
      for (j=0; j<=MAXR; j++) {
	if (
	    (i==j && fabs(check[i][j]-1.0) > 1.0e-10) ||
	    (i!=j && fabs(check[i][j]) > 1.0e-10)
	    )
	  cout << "Fail: " 
	       << setw(5) << iz 
	       << setw(5) << i 
	       << setw(5) << j
	       << setw(15) << check[i][j] << endl;
      }
    }
  */

  }

  if (use_tables) pt.make_table(NPT, NPTRMAX, this);

  initialized = TRUE;
}


double CylindricalCB::get_pot(int mu, int nn, double r, double z)
{
  if (mu > MAXR) bomb ("Order must be smaller than max, reset");
  if (nn > n-1) bomb ("Order must be smaller than max, reset");
  if (initialized == FALSE) initialize();

  double r2 = r*r;
  double A2 = A*A;
  double t = r/sqrt(A2 + r2);
  double ret = 0.0;
  int k;

  double trig;
  if (nn<n2) {
    k = nn;
    trig = cos(dk*k*z);
  }
  else {
    k = nn-n2;
    trig = sin(dk*k*z);
  }

  ret = 0.0;
  for (int i=MAXR; i>=0; i--) ret = UT[k][mu][i] + t*ret;
  ret *= pow(A2 + r2, -0.5*(1.0 + M)) * Z[k][mu] * trig;
    
  return ret;
}


double CylindricalCB::get_dens(int mu, int nn, double r, double z)
{
  if (mu > MAXR) bomb ("Order must be smaller than max, reset");
  if (nn > n-1) bomb ("Order must be smaller than max, reset");
  if (initialized == FALSE) initialize();

  double A2 = A*A, M1 = 1.0 + M;
  double A4 = A2*A2;
  double t = r/sqrt(A2 + r*r);
  double r2 = r*r;
  double trig, ret = 0.0;
  int k;

  if (nn<n2) {
    k = nn;
    trig = cos(dk*k*z);
  }
  else {
    k = nn-n2;
    trig = sin(dk*k*z);
  }

  double factor = (A2 + r2)*(A2 + r2) * ( M*M/(r2+TINY) + dk*dk*k*k );

  for (int i=MAXR; i>=0; i--)
    ret = UT[k][mu][i] * 
      ( A4*i*i/(r2+TINY) - 
	 2.0*((1.0 + i)*M1 + i)*A2  + M1*M1*r2 - factor) + t*ret;

  ret *= pow(A2 + r2, -0.5*(5.0 + M)) * Z[k][mu] * trig;

  return 0.25*ret/M_PI;
}


double CylindricalCB::norm(int mu, int nn)
{
  if (mu > MAXR) bomb ("Order must be smaller than max, reset");
  if (nn > n-1) bomb ("Order must be smaller than max, reset");
  if (initialized == FALSE) initialize();

  return Z[mu][nn];
}

Matrix CylindricalCB::inner_product_dens(double (*fct)(double, double))
{
  Matrix coefs(0, n-1, 0, MAXR);
  Vector pos(0, n2-1);
  Vector neg(0, n2-1);
  Vector posin(0, n2-1);
  Vector negin(0, n2-1);
  Vector input(0, n-1);
  Vector hold(0, MAXR);
  coefs.zero();

  LegeQuad lwk(NINT);
  double A2 = A*A, M1 = 1.0+M;
  double A4 = A2*A2; 
  double tz, z, tr, r, r2;
  double wghtr, fctp, fctn, offset;

  int iz, mu, nn;

  for (int kr=1; kr<=NINT; kr++) {

    tr = lwk.knot(kr);
    r = A*tr/sqrt(1.0 - tr*tr);
    r2 = r*r;
    wghtr = A*pow(1.0 - tr*tr, -1.5) * r * lwk.weight(kr) *
      pow(A2 + r2, -0.5*(1.0 + M)) ;

    for (iz=0; iz<n2; iz++) {
      z = dZ*iz;
      input[iz] = fct(r, z);
      if (iz) input[n-iz] = fct(r, -z);
    }
    input[n2] = fct(r, dZ*n2);

    fft->input(input);
    fft->real_half(pos);
    fft->imag_half(neg);

    pos /= n2;
    neg /= n2;

    for (iz=0; iz<n2; iz++) {

      hold.zero();
      for (mu=0; mu<=MAXR; mu++) {
	for (int j=MAXR; j>=0; j--)
	  hold[mu] = UT[iz][mu][j] + tr*hold[mu];

	hold[mu] *= Z[iz][mu];
      }

      coefs[iz] += hold * wghtr * pos[iz];
      coefs[iz+n2] += hold * wghtr * neg[iz];
    }
  }

  return 4.0*M_PI*coefs;
}

Matrix CylindricalCB::inner_product_pot(double (*fct)(double, double))
{
  Matrix coefs(0, n-1, 0, MAXR);
  Vector pos(0, n2-1);
  Vector neg(0, n2-1);
  Vector input(0, n-1);
  Vector hold(0, MAXR);
  coefs.zero();

  LegeQuad lwk(NINT);
  double A2 = A*A, M1 = 1.0+M;
  double A4 = A2*A2; 
  double tz, z, tr, r, r2;
  double wghtr, fctp, fctn, offset, factor;

  int iz, mu, nn;

  for (int kr=1; kr<=NINT; kr++) {

    tr = lwk.knot(kr);
    r = A*tr/sqrt(1.0 - tr*tr);
    r2 = r*r;
    wghtr = A*pow(1.0 - tr*tr, -1.5) * r * lwk.weight(kr) *
      pow(A2 + r2, -0.5*(5.0 + M)) ;

    for (iz=0; iz<n2; iz++) {
      z = dZ*iz;
      input[iz] = fct(r, z);
      if (iz) input[n-iz] = fct(r, -z);
    }
    input[n2] = fct(r, dZ*n2);

    fft->input(input);
    fft->real_half(pos);
    fft->imag_half(neg);

    pos /= n2;
    neg /= n2;

    for (iz=0; iz<n2; iz++) {

      factor = (A2 + r2)*(A2 + r2) * (dk*dk*iz*iz + M*M/(r2+TINY));

      hold.zero();
      for (mu=0; mu<=MAXR; mu++) {
	for (int j=MAXR; j>=0; j--)
	  hold[mu] = UT[iz][mu][j] *
	    ( A4*j*j/(r2+TINY) - 
	      2.0*((1.0 + j)*M1 + j)*A2  + M1*M1*r2 - factor )
	    + tr*hold[mu];
	
	hold[mu] *= Z[iz][mu];
      }

      coefs[iz] += hold * wghtr * pos[iz];
      coefs[iz+n2] += hold * wghtr * neg[iz];
    }
  }

  return coefs;
}


double CylindricalCB::pot_eval(Matrix& coefs, double r, double z)
{
  Vector hold(0, MAXR);

  int mu, nn;

  double A2 = A*A, M1 = 1.0+M;
  double r2=r*r;
  double t = r/sqrt(A2 + r2);
  double ans = 0.0;
  double k;

  double dcos = cos(dk*z);
  double dsin = sin(dk*z);
  double lcos, lsin, cosn=1.0, sinn=0.0;

  for (nn=0; nn<n2; nn++) {

    k = dk*nn;

    hold.zero();
    for (int j=MAXR; j>=0; j--) {

      for (mu=0; mu<=MAXR; mu++)
	hold[mu] = UT[nn][mu][j] + t*hold[mu];
    }

    for (mu=0; mu<=MAXR; mu++) 
      ans += (coefs[nn][mu]*cosn + coefs[nn+n2][mu]*sinn) * hold[mu] *
	Z[nn][mu];
	
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;

  }

  return -ans*pow(A2 + r2, -0.5*(1.0 + M));
}


double CylindricalCB::dens_eval(Matrix& coefs, double r, double z)
{

  Vector hold(0, MAXR);

  int mu, nn;

  double A2 = A*A, M1 = 1.0+M;
  double A4 = A2*A2; 
  double r2 = r*r;
  double t = r/sqrt(A2 + r2);
  double ans = 0.0, factor, dfac;
  double k;

  double dcos = cos(dk*z);
  double dsin = sin(dk*z);
  double lcos, lsin, cosn=1.0, sinn=0.0;

  for (nn=0; nn<n2; nn++) {

    k = dk*nn;

    factor = (A2 + r2)*(A2 + r2) * ( M*M/(r2+TINY) + k*k );
    
    hold.zero();
    for (int j=MAXR; j>=0; j--) {

      dfac = A4*j*j/(r2+TINY) - 
	2.0*((1.0 + j)*M1 + j)*A2  + M1*M1*r2 - factor;

      for (mu=0; mu<=MAXR; mu++)
	hold[mu] = UT[nn][mu][j]*dfac + t*hold[mu];
    }

    for (mu=0; mu<=MAXR; mu++) 
      ans += (coefs[nn][mu]*cosn + coefs[nn+n2][mu]*sinn) * hold[mu] *
	Z[nn][mu];
	
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  return -ans * pow(A2 + r2, -0.5*(5.0 + M)) * 0.25/M_PI;
}



double CylindricalCB::r_force_eval(Matrix& coefs, double r, double z)
{

  Vector hold(0, MAXR);

  int mu, nn;

  double A2 = A*A, M1 = 1.0+M;
  double A4 = A2*A2; 
  double r2 = r*r;
  double t = r/sqrt(A2 + r2);
  double ans = 0.0;
  double k, ffac;

  double dcos = cos(dk*z);
  double dsin = sin(dk*z);
  double lcos, lsin, cosn=1.0, sinn=0.0;

  for (nn=0; nn<n2; nn++) {

    k = dk*nn;

    hold.zero();
    for (int j=MAXR; j>=0; j--) {

      ffac = A2*j/(r+TINY) - M1*r;

      for (mu=0; mu<=MAXR; mu++)
	hold[mu] = UT[nn][mu][j]*ffac + t*hold[mu];
    }

    for (mu=0; mu<=MAXR; mu++) 
      ans += (coefs[nn][mu]*cosn + coefs[nn+n2][mu]*sinn) * hold[mu] *
	Z[nn][mu];
	
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  return ans * pow(A2 + r2, -0.5*(3.0 + M));
}


double CylindricalCB::z_force_eval(Matrix& coefs, double r, double z)
{
  Vector hold(0, MAXR);

  int mu, nn;

  double A2 = A*A, M1 = 1.0+M;
  double r2=r*r;
  double t = r/sqrt(A2 + r2);
  double ans = 0.0;
  double k;

  double dcos = cos(dk*z);
  double dsin = sin(dk*z);
  double lcos, lsin, cosn=1.0, sinn=0.0;
  for (nn=0; nn<n2; nn++) {

    k = dk*nn;

    hold.zero();
    for (int j=MAXR; j>=0; j--) {
      for (mu=0; mu<=MAXR; mu++)
	hold[mu] = UT[nn][mu][j] + t*hold[mu];
    }

    for (mu=0; mu<=MAXR; mu++) 
      ans += (-coefs[nn][mu]*sinn + coefs[nn+n2][mu]*cosn) * k * hold[mu] *
	Z[nn][mu];
	
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;

  }

  return ans*pow(A2 + r2, -0.5*(1.0 + M));
}


void CylindricalCB::force_eval(Matrix& coefs, double r, double z,
				 double& fr, double& fz)
{

  Vector holdr(0, MAXR);
  Vector holdz(0, MAXR);

  int mu, nn, j;

  double A2 = A*A, M1 = 1.0+M;
  double A4 = A2*A2; 
  double r2 = r*r;
  double t = r/sqrt(A2 + r2);
  double k, ffac;

  double dcos = cos(dk*z);
  double dsin = sin(dk*z);
  double lcos, lsin, cosn=1.0, sinn=0.0;

  fr = 0.0;
  fz = 0.0;

  /*

  for (nn=0; nn<n2; nn++) {

    k = dk*nn;

    holdr.zero();
    holdz.zero();
    for (int j=MAXR; j>=0; j--) {

      ffac = A2*j/(r+TINY) - M1*r;

      for (mu=0; mu<=MAXR; mu++) {
	holdr[mu] = UT[nn][mu][j]*ffac + t*holdr[mu];
	holdz[mu] = UT[nn][mu][j] + t*holdz[mu];
      }
    }

    for (mu=0; mu<=MAXR; mu++) {
      fr += (coefs[nn][mu]*cosn + coefs[nn+n2][mu]*sinn) * holdr[mu] *
	Z[nn][mu];
      fz += (-coefs[nn][mu]*sinn + coefs[nn+n2][mu]*cosn) * k * holdz[mu] *
	Z[nn][mu];
	
    }

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  */

  holdr[0] = -M1*r;
  holdz[0] = 1.0;

  for (j=1; j<=MAXR; j++) {
    holdz[j] = holdz[j-1] * t;
    holdr[j] = holdz[j] * (A2*j/(r+TINY) - M1*r);
  }

  double facr, facz;

  for (nn=0; nn<n2; nn++) {

    k = dk*nn;

    for (mu=0; mu<=MAXR; mu++) {

      facr = facz = 0.0;
      for (j=0; j<=MAXR; j++) {
	facr += UT[nn][mu][j] * holdr[j];
	facz += UT[nn][mu][j] * holdz[j];
      }

      fr += (coefs[nn][mu]*cosn + coefs[nn+n2][mu]*sinn) * Z[nn][mu] * facr;
      fz += (-coefs[nn][mu]*sinn + coefs[nn+n2][mu]*cosn) * k * Z[nn][mu] * facz;
    }

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  fr *= pow(A2 + r2, -0.5*(3.0 + M));
  fz *= pow(A2 + r2, -0.5*(1.0 + M));

}

void CylindricalCB::setup_accumulation(void)
{
  accum_cos.setsize(0, n-1, 0, MAXR);
  accum_cos.zero();
  if (M>0) {
    accum_sin.setsize(0, n-1, 0, MAXR);
    accum_sin.zero();
  }
}

void CylindricalCB::accumulate(double r, double z, double phi, double mass)
{
  Vector hold(0, MAXR);

  int mu, nn;

  double A2 = A*A, M1 = 1.0+M;
  double r2=r*r;
  double t = r/sqrt(A2 + r2);
  double ans = 0.0;
  double k;

  double dcos = cos(dk*z);
  double dsin = sin(dk*z);
  double lcos, lsin, cosn=1.0, sinn=0.0;

  double fac, tfac = pow(A2 + r2, -0.5*(1.0 + M)) * mass * 4.0*M_PI;
  double pfac = 1.0/ZMAX;
  double cfac = cos(phi*M) * tfac;
  double sfac = sin(phi*M) * tfac;

  /*

  for (nn=0; nn<n2; nn++) {

    k = dk*nn;
    tfac = pfac;
    if (nn==0) tfac *= 0.5;

    hold.zero();
    for (int j=MAXR; j>=0; j--) {

      for (mu=0; mu<=MAXR; mu++)
	hold[mu] = UT[nn][mu][j] + t*hold[mu];
    }

    for (mu=0; mu<=MAXR; mu++) {
      fac = cfac * tfac * hold[mu]*Z[nn][mu];
      accum_cos[nn][mu] += fac * cosn;
      accum_cos[nn+n2][mu] += fac * sinn;
      if (M>0) {
	fac = sfac * tfac * hold[mu]*Z[nn][mu];
	accum_sin[nn][mu] += fac * cosn;
	accum_sin[nn+n2][mu] += fac * sinn;
      }
    }

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;

  }

  */

  hold[0] = 1.0;
  for (mu=1; mu<=MAXR; mu++) hold[mu] = hold[mu-1] * t;

  for (nn=0; nn<n2; nn++) {

    k = dk*nn;
    tfac = pfac;
    if (nn==0) tfac *= 0.5;

    for (mu=0; mu<=MAXR; mu++) {
      fac = cfac * tfac * hold[mu];
      accum_cos[nn][mu] += fac * cosn;
      accum_cos[nn+n2][mu] += fac * sinn;
      if (M>0) {
	fac = sfac * tfac * hold[mu];
	accum_sin[nn][mu] += fac * cosn;
	accum_sin[nn+n2][mu] += fac * sinn;
      }
    }

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;

  }

  coefs_made = FALSE;
}

Matrix CylindricalCB::accumulated_coefficients(int sincos)
{
  /*
  if (M>0) {
    if (sincos)
      return accum_sin/M_PI;
    else 
      return accum_cos/M_PI;
  }
  else
    return 0.5*accum_cos/M_PI;
    */

  int nn, j, k;
  Matrix ret(0, n-1, 0, MAXR);
  ret.zero();
  double norm;

  if (M>0) {
    if (sincos) {

      for (nn=0; nn<n2; nn++) {
	for (j=0; j<=MAXR; j++) {
	  norm = Z[nn][j]/M_PI;
	  for (k=0; k<=MAXR; k++) ret[nn][j] += 
				    UT[nn][j][k] * accum_sin[nn][k] * norm;
	}
      }

      return ret;
    }
    else {
      for (nn=0; nn<n2; nn++) {
	for (j=0; j<=MAXR; j++) {
	  norm = Z[nn][j]/M_PI;
	  for (k=0; k<=MAXR; k++) ret[nn][j] += 
				    UT[nn][j][k] * accum_cos[nn][k] * norm;
	}
      }

      return ret;
    }
  }
  else {

    for (nn=0; nn<n2; nn++) {
      for (j=0; j<=MAXR; j++) {
	norm = 0.5*Z[nn][j]/M_PI;
	for (k=0; k<=MAXR; k++) ret[nn][j] += 
				  UT[nn][j][k] * accum_cos[nn][k] * norm;
      }
    }

    return ret;
  }
}

void CylindricalCB::make_coefficients(void)
{
  int nn, j, k;
  double norm;

  if (!MPIset) {
    MPIin = new double [n*(MAXR+1)];
    MPIout = new double [n*(MAXR+1)];
    MPIset = TRUE;
  }
  
#ifdef MPE_PROFILE
    MPE_Log_event(7, myid, "b_distrib_c");
#endif

  for (nn=0; nn<n; nn++)
    for (j=0; j<=MAXR; j++) MPIin[nn*(MAXR+1)+j] = accum_cos[nn][j];

  MPI_Allreduce ( MPIin, MPIout, n*(MAXR+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (nn=0; nn<n; nn++)
    for (j=0; j<=MAXR; j++) accum_cos[nn][j] = MPIout[nn*(MAXR+1)+j];

  if (M>0) {

    for (nn=0; nn<n; nn++)
      for (j=0; j<=MAXR; j++) MPIin[nn*(MAXR+1)+j] = accum_sin[nn][j];

    MPI_Allreduce ( MPIin, MPIout, n*(MAXR+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (nn=0; nn<n; nn++)
      for (j=0; j<=MAXR; j++) accum_sin[nn][j] = MPIout[nn*(MAXR+1)+j];

  }

#ifdef MPE_PROFILE
    MPE_Log_event(8, myid, "e_distrib_c");
#endif

  coscoef.setsize(0, n-1, 0, MAXR);
  coscoef.zero();
  if (M>0) {
    sincoef.setsize(0, n-1, 0, MAXR);
    sincoef.zero();
  }

  for (nn=0; nn<n2; nn++) {
    for (j=0; j<=MAXR; j++) {
      norm = Z[nn][j]/M_PI;
      for (k=0; k<=MAXR; k++) {
	coscoef[nn][j] += UT[nn][j][k] * accum_cos[nn][k] * norm;
	coscoef[n2+nn][j] += UT[nn][j][k] * accum_cos[n2+nn][k] * norm;
	if (M>0) {
	  sincoef[nn][j] += UT[nn][j][k] * accum_sin[nn][k] * norm;
	  sincoef[n2+nn][j] += UT[nn][j][k] * accum_sin[n2+nn][k] * norm;
	}
      }
    }
  }

  coefs_made = TRUE;
}

  
void CylindricalCB::accumulated_eval(double r, double z, double phi,
				     double& p, double& fr, double& fz, 
				     double &fp)
{

  if (!coefs_made) make_coefficients();

  if (use_tables && r<=NPTRMAX) {
    pt.feval(r, z, phi, p, fr, fz, fp);
    return;
  }

  Vector holdr(0, MAXR);
  Vector holdz(0, MAXR);

  int mu, nn, j;

  double A2 = A*A, M1 = 1.0+M;
  double A4 = A2*A2; 
  double r2 = r*r;
  double t = r/sqrt(A2 + r2);
  double k, ffac;

  double dcos = cos(dk*z);
  double dsin = sin(dk*z);
  double lcos, lsin, cosn=1.0, sinn=0.0;

  double cosp = cos(phi*M);
  double sinp = sin(phi*M);

  fr = 0.0;
  fz = 0.0;
  fp = 0.0;
  p = 0.0;

  holdr[0] = -M1*r;
  holdz[0] = 1.0;

  for (j=1; j<=MAXR; j++) {
    holdz[j] = holdz[j-1] * t;
    holdr[j] = holdz[j] * (A2*j/(r+TINY) - M1*r);
  }

  double facr, facz;

  for (nn=0; nn<n2; nn++) {

    k = dk*nn;

    for (mu=0; mu<=MAXR; mu++) {

      facr = facz = 0.0;
      for (j=0; j<=MAXR; j++) {
	facr += UT[nn][mu][j] * holdr[j];
	facz += UT[nn][mu][j] * holdz[j];
      }

      if (M==0) {

	fr += (coscoef[nn][mu]*cosn + coscoef[nn+n2][mu]*sinn) * 
	  Z[nn][mu] * facr;
	fz += (-coscoef[nn][mu]*sinn + coscoef[nn+n2][mu]*cosn) * 
	  k * Z[nn][mu] * facz;
	p += (coscoef[nn][mu]*cosn + coscoef[nn+n2][mu]*sinn) * 
	  Z[nn][mu] * facz;

      }
      else {

	fr += ( cosp*(coscoef[nn][mu]*cosn + coscoef[nn+n2][mu]*sinn) +
		sinp*(sincoef[nn][mu]*cosn + sincoef[nn+n2][mu]*sinn) ) *
	  Z[nn][mu] * facr;

	fz += ( cosp*(-coscoef[nn][mu]*sinn + coscoef[nn+n2][mu]*cosn) +
		sinp*(-sincoef[nn][mu]*sinn + sincoef[nn+n2][mu]*cosn) ) *
	  k * Z[nn][mu] * facz;

	fp += ( sinp*(coscoef[nn][mu]*cosn + coscoef[nn+n2][mu]*sinn) -
		cosp*(sincoef[nn][mu]*cosn + sincoef[nn+n2][mu]*sinn) ) * M *
	  Z[nn][mu] * facz;

	p += ( cosp*(coscoef[nn][mu]*cosn + coscoef[nn+n2][mu]*sinn) +
	       sinp*(sincoef[nn][mu]*cosn + sincoef[nn+n2][mu]*sinn) ) *
	  Z[nn][mu] * facz;

      }
    }

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  fr *= pow(A2 + r2, -0.5*(3.0 + M));
  fz *= pow(A2 + r2, -0.5*(1.0 + M));
  fp *= pow(A2 + r2, -0.5*(1.0 + M));
  p *= pow(A2 + r2, -0.5*(1.0 + M));

}

double CylindricalCB::accumulated_dens_eval(double r, double z, double phi)
{
  if (!coefs_made) make_coefficients();

  Vector hold(0, MAXR);

  int mu, nn;

  double A2 = A*A, M1 = 1.0+M;
  double A4 = A2*A2; 
  double r2 = r*r;
  double t = r/sqrt(A2 + r2);
  double ans = 0.0, factor, dfac;
  double k;

  double dcos = cos(dk*z);
  double dsin = sin(dk*z);
  double lcos, lsin, cosn=1.0, sinn=0.0;

  double cosp = cos(phi*M);
  double sinp = sin(phi*M);


  for (nn=0; nn<n2; nn++) {

    k = dk*nn;

    factor = (A2 + r2)*(A2 + r2) * ( M*M/(r2+TINY) + k*k );
    
    hold.zero();
    for (int j=MAXR; j>=0; j--) {

      dfac = A4*j*j/(r2+TINY) - 
	2.0*((1.0 + j)*M1 + j)*A2  + M1*M1*r2 - factor;

      for (mu=0; mu<=MAXR; mu++)
	hold[mu] = UT[nn][mu][j]*dfac + t*hold[mu];
    }


    for (mu=0; mu<=MAXR; mu++) {

      if (M==0) {
	ans += (coscoef[nn][mu]*cosn + coscoef[nn+n2][mu]*sinn) * hold[mu] *
	  Z[nn][mu];
      }
      else {
	ans += ( cosp*(coscoef[nn][mu]*cosn + coscoef[nn+n2][mu]*sinn) +
	       sinp*(sincoef[nn][mu]*cosn + sincoef[nn+n2][mu]*sinn) ) *
	  hold[mu] * Z[nn][mu];

      }

    }
    
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  return -ans * pow(A2 + r2, -0.5*(5.0 + M)) * 0.25/M_PI;
}


void CylindricalCB::PotTable::make_table(int nin, double rmax, 
					  CylindricalCB* p)
{

  int mu, nn, i, j;

  N = nin;
  RMAX = rmax;
  dR = RMAX/N;
  N2 = p->n2;
  NR = p->MAXR;
  ref = p;

  double r;

  fr = new Vector* [N2];
  fz = new Vector* [N2];
 
  for (nn=0; nn<N2; nn++) {
    
    fr[nn] = new Vector [NR+1];
    fz[nn] = new Vector [NR+1];
 
    for (mu=0; mu<=NR; mu++) {
      
      fr[nn][mu].setsize(0, N);
      fz[nn][mu].setsize(0, N);
 
      fr[nn][mu].zero();
      fz[nn][mu].zero();
    }
  }


  Vector holdr(0, NR);
  Vector holdz(0, NR);

  double A2 = p->A*p->A, M1 = 1.0+p->M;
  double A4 = A2*A2; 
  double r2, t;

  double facr, facz, powr, powz;

  for (i=0; i<=N; i++) {
    
    r = dR*i;

    r2 = r*r;
    t = r/sqrt(A2 + r2);

    powr =  pow(A2 + r2, -0.5*(3.0 + p->M));
    powz =  pow(A2 + r2, -0.5*(1.0 + p->M));
    
    holdr[0] = -M1*r;
    holdz[0] = 1.0;

    for (j=1; j<=NR; j++) {
      holdz[j] = holdz[j-1] * t;
      holdr[j] = holdz[j] * (A2*j/(r+TINY) - M1*r);
    }

    for (nn=0; nn<N2; nn++) {

      for (mu=0; mu<=NR; mu++) {

	facr = facz = 0.0;
	for (j=0; j<=NR; j++) {
	  facr += p->UT[nn][mu][j] * holdr[j];
	  facz += p->UT[nn][mu][j] * holdz[j];
	}

	fr[nn][mu][i] = p->Z[nn][mu] * facr * powr;
	fz[nn][mu][i] = p->Z[nn][mu] * facz * powz;

      }

    }

  }
}


CylindricalCB::PotTable::~PotTable(void)
{
  for (int nn=0; nn<N2; nn++) {
    delete [] fr[nn];
    delete [] fz[nn];
  }

  delete [] fr;
  delete [] fz;
}


void CylindricalCB::PotTable::feval(double r, double z, double phi,
				    double& P, double& FR, double& FZ, 
				    double &FP)
{

  if (!N) {
    cerr << "PotTable not initialized\n";
    exit(-3);
  }

  int mu, nn, j;

  int indx = (int)(r/dR);
  if (indx>=N) indx = N-1;
  double fac1 = 1.0 + indx - r/dR;
  double fac2 = 1.0 - fac1;

  double dcos = cos(ref->dk*z);
  double dsin = sin(ref->dk*z);
  double lcos, lsin, cosn=1.0, sinn=0.0;

  int M = ref->M;
  int n2 = ref->n2;

  double cosp = cos(phi*M);
  double sinp = sin(phi*M);
  double k, facr, facz, fac0, fcos, fsin;

  FR = 0.0;
  FZ = 0.0;
  FP = 0.0;
  P = 0.0;

  for (nn=0; nn<ref->n2; nn++) {

    k = ref->dk*nn;

    for (mu=0; mu<=ref->MAXR; mu++) {

      facr = fac1*fr[nn][mu][indx] + fac2*fr[nn][mu][indx+1];
      facz = fac1*fz[nn][mu][indx] + fac2*fz[nn][mu][indx+1];

      if (M==0) {

	fac0 = ref->coscoef[nn][mu]*cosn + ref->coscoef[nn+n2][mu]*sinn;

	FR += fac0 * facr;

	FZ += (-ref->coscoef[nn][mu]*sinn + 
	       ref->coscoef[nn+n2][mu]*cosn) *  k * facz;

	P -= fac0 * facz;

      }
      else {

	fcos =  ref->coscoef[nn][mu]*cosn + 
    	        ref->coscoef[nn+n2][mu]*sinn;

	fsin =  ref->sincoef[nn][mu]*cosn + 
	        ref->sincoef[nn+n2][mu]*sinn;

	fac0 =  cosp*fcos + sinp*fsin;

	FR += fac0 * facr;

	FZ += ( cosp*(-ref->coscoef[nn][mu]*sinn + 
		      ref->coscoef[nn+n2][mu]*cosn) +
		sinp*(-ref->sincoef[nn][mu]*sinn + 
		      ref->sincoef[nn+n2][mu]*cosn) ) * k * facz;

	FP += ( sinp*fcos - cosp*fsin ) * M * facz;

	P -= fac0 * facz;
      }
    }

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

}

