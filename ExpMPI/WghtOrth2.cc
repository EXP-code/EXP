#define DEBUG 1
#define CHECK_OF

#include <iostream.h>
#include <iomanip.h>
#include <String.h>

#include <numerical.h>
#include <gaussQ.h>
#include "WghtOrth2.h"


#undef TINY
#define TINY 1.0e-16


double CylindricalSL::RMIN = 0.001;
double CylindricalSL::RMAX = 100.0;
int CylindricalSL::NUMR = 1000;
int CylindricalSL::NINT = 200;
double CylindricalSL::ZMAX = 5.0;


CylindricalSL::CylindricalSL(void)
{
  MMAX = 0;
  NMAX = 0;
  NFFT = 0;
  initialized = FALSE;
  MPIset = FALSE;
  coefs_made = FALSE;

  accum_cos = 0;
  accum_sin = 0;
  table = 0;
  tablef = 0;
}

CylindricalSL::~CylindricalSL(void)
{
  if (initialized) {
    delete ortho;
    delete fft;
  }

  delete [] accum_cos;
  delete [] accum_sin;
  delete [] table;
  delete [] tablef;

  if (MPIset) {
    delete [] MPIin;
    delete [] MPIout;
  }
}

CylindricalSL::CylindricalSL(int mu, int nfft, int m)
{
  NMAX = mu;
  NFFT = nfft;
  MMAX = m;

  initialized = FALSE;
  MPIset = FALSE;
  coefs_made = FALSE;

  accum_cos = 0;
  accum_sin = 0;
  table = 0;
  tablef = 0;
}

void CylindricalSL::reset(int mu, int nfft, int m)
{
  NMAX = mu;
  NFFT = nfft;
  MMAX = m;

  if (initialized) {
    delete ortho;
    delete fft;
  }

  delete [] accum_cos;
  delete [] accum_sin;
  delete [] table;
  delete [] tablef;

  accum_cos = 0;
  accum_sin = 0;
  table = 0;
  tablef = 0;

  if (MPIset) {
    delete [] MPIin;
    delete [] MPIout;
  }

  initialize();
}


void CylindricalSL::initialize(void)
{
  if (MMAX<0)  bomb ("MMAX must be >= 0");
  if (NMAX<0)  bomb ("NMAX must be >= 0");
  if (NFFT<=0)  bomb ("NFFT must be > 0");

  //
  // Setup Fourier transform
  //

  n = 2<<(NFFT-1);

				// Grid meshes
  dZ = 2.0*ZMAX/n;
  dk = 0.5*M_PI/ZMAX;

  fft = new FFT (2*n, dZ);


  //
  // Initialize grid
  //

#ifdef DEBUG
  cout.form("Process %d: MMAX=%d  NMAX=%d  NUMR=%d  NUMK=%d  RMIN=%f  RMAX=%f  L=%f\n", myid, MMAX, NMAX, NUMR, n, RMIN, RMAX, ZMAX);
#endif


  SLGridCyl::mpi = 1;		// Turn on MPI
  ortho = new SLGridCyl(MMAX, NMAX, NUMR, n, RMIN, RMAX, ZMAX);

  initialized = TRUE;
}


double CylindricalSL::get_pot(int mm, int nn, int kk, double r, double z)
{
  if (nn > NMAX) bomb ("Order must be smaller than max, reset");
  if (kk > n-1) bomb ("Order must be smaller than max, reset");
  if (initialized == FALSE) initialize();

  if (fabs(z)>ZMAX) return 0.0;

  return ortho->get_pot(r, mm, nn, kk) * sin(dk*kk*(z + ZMAX))/sqrt(ZMAX);
}


double CylindricalSL::get_dens(int mm, int kk, int nn, double r, double z)
{
  if (nn > NMAX) bomb ("Order must be smaller than max, reset");
  if (kk > n-1) bomb ("Order must be smaller than max, reset");
  if (initialized == FALSE) initialize();

  if (fabs(z)>ZMAX) return 0.0;

  return 0.25/M_PI * 
    ortho->get_dens(r, mm, nn, kk) * sin(dk*kk*(z + ZMAX))/sqrt(ZMAX);
}


Matrix CylindricalSL::inner_product_dens(int m, double (*fct)(double, double))
{
  if (initialized == FALSE) initialize();

  Vector fk(0, n-1);
  Matrix fr(1, n-1, 1, NMAX);
  Vector fz(0, n-1);

  Matrix coefs(1, n-1, 1, NMAX);
  coefs.zero();

  LegeQuad lwk(NINT);

  double xi, r, wghtr, z;
  int iz;


  for (int kr=1; kr<=NINT; kr++) {

    xi = 2.0*(lwk.knot(kr) - 0.5);
    r = ortho->xi_to_r(xi);
    wghtr = 1.0/ortho->d_xi_to_r(xi) * r * 2.0 * lwk.weight(kr);

    for (iz=0; iz<n; iz++) {
      z = dZ*iz - ZMAX;
      fz[iz] = (*fct)(r, z);
    }

    fft->sin_transform(fk, fz);
    // Test
    /*
    {
      double ans;

      cerr << "r: " << r << endl;

      for (int k=1; k<n; k++) {
	ans = 0.0;
	for (iz=1; iz<n; iz++) {
	  z = dZ*iz - ZMAX;
	  ans += dZ * sin(dk*k*(z+ZMAX))*fz[iz];
	}
	cerr << setw(5) << k << setw(15) << fk[k] << setw(15) << ans
	     << endl;
      }
      cerr << endl;
    }
    */
    // End test
    ortho->get_pot(fr, xi, m, 0);

    for (iz=1; iz<n; iz++) coefs[iz] += fr[iz] * wghtr * fk[iz];
  }

  //
  // Norm----|
  //         |
  //         v
  return -4.0*M_PI / sqrt(ZMAX) * coefs;
}

Matrix CylindricalSL::inner_product_pot(int m, double (*fct)(double, double))
{
  if (initialized == FALSE) initialize();

  Vector fk(0, n);
  Matrix fr(1, n, 1, NMAX);
  Vector fz(0, n);

  Matrix coefs(1, n, 1, NMAX);
  coefs.zero();

  LegeQuad lwk(NINT);

  double xi, r, z, wghtr;
  int iz;


  for (int kr=1; kr<=NINT; kr++) {

    xi = 2.0*(lwk.knot(kr) - 0.5);
    r = ortho->xi_to_r(xi);
    wghtr = 1.0/ortho->d_xi_to_r(xi) * r * 2.0 * lwk.weight(kr);

    for (iz=0; iz<n; iz++) {
      z = dZ*iz - ZMAX;
      fz[iz] = fct(r, z);
    }

    fft->sin_transform(fk, fz);
    ortho->get_dens(fr, xi, m, 0);

    for (iz=1; iz<n; iz++) coefs[iz] += fr[iz] * wghtr * fk[iz];
  }

  //
  //     |---Norm
  //     |
  //     v
  return -coefs / sqrt(ZMAX);
}


double CylindricalSL::pot_eval(Matrix& coefs, int m, double r, double z)
{
  if (initialized == FALSE) initialize();

  double k, ans=0.0;
  int nk;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn=1, sinn=0.0;

  Matrix tab;

  ortho->get_pot(tab, r, m);

  for (nk=1; nk<n; nk++) {

    k = dk*nk;

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;

    ans += (coefs[nk] * tab[nk]) * sinn / sqrt(ZMAX);
    
  }

  return ans;
}


double CylindricalSL::dens_eval(Matrix& coefs, int m, double r, double z)
{
  if (initialized == FALSE) initialize();

  double k, ans=0.0;
  int nk;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn=1, sinn=0.0;

  Matrix tab;

  ortho->get_dens(tab, r, m);

  for (nk=1; nk<n; nk++) {

    k = dk*nk;

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;

    ans += (coefs[nk] * tab[nk])* sinn / sqrt(ZMAX);
    
  }

  return ans * 0.25/M_PI;
}



double CylindricalSL::r_force_eval(Matrix& coefs, int m, double r, double z)
{
  if (initialized == FALSE) initialize();

  double k, ans=0.0;
  int nk;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn=1, sinn=0.0;

  Matrix tab;

  ortho->get_force(tab, r, m);

  for (nk=1; nk<n; nk++) {

    k = dk*nk;

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
    
    ans += (coefs[nk] * tab[nk])* sinn / sqrt(ZMAX);
    
  }

  return -ans;
}


double CylindricalSL::z_force_eval(Matrix& coefs, int m, double r, double z)
{
  if (initialized == FALSE) initialize();

  double k, ans=0.0;
  int nk;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn=1, sinn=0.0;

  Matrix tab;

  ortho->get_pot(tab, r, m);

  for (nk=1; nk<n; nk++) {

    k = dk*nk;

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
    
    ans += (coefs[nk] * tab[nk])* k*cosn / sqrt(ZMAX);
    
  }

  return -ans;
}


void CylindricalSL::force_eval(Matrix& coefs, int m, double r, double z,
				 double& fr, double& fz)
{
  if (initialized == FALSE) initialize();

  double k, ansr=0.0, ansz=0.0;
  int nk;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn=1, sinn=0.0;

  Matrix tabp, tabf;

  ortho->get_pot(tabp, r, m);
  ortho->get_force(tabf, r, m);

  for (nk=1; nk<n; nk++) {

    k = dk*nk;

    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
    
    ansr += (coefs[nk] * tabf[nk])* sinn / sqrt(ZMAX);
    ansz += (coefs[nk] * tabp[nk])* k*cosn / sqrt(ZMAX);
    
  }

  fr = -ansr;
  fz = -ansz;

}

void CylindricalSL::setup_accumulation(void)
{
  if (initialized == FALSE) initialize();

  if (!accum_cos) {
    accum_cos = new Matrix [MMAX+1];
    accum_sin = new Matrix [MMAX+1];
    table = new Matrix [MMAX+1];
    tablef = new Matrix [MMAX+1];
  }

  for (int m=0; m<=MMAX; m++) {
    accum_cos[m].setsize(1, n, 1, NMAX);
    accum_cos[m].zero();
    if (m>0) {
      accum_sin[m].setsize(1, n, 1, NMAX);
      accum_sin[m].zero();
    }
  }

  coefs_made = FALSE;

}

void CylindricalSL::accumulate(double r, double z, double phi, double mass)
{
  if (coefs_made) setup_accumulation();

  double k, msin, mcos;
  int nk, mm;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn, sinn;
  double norm = -4.0*M_PI;

  ortho->get_pot(table, r, 0, MMAX);

  for (mm=0; mm<=MMAX; mm++) {

    mcos = cos(phi*mm);
    msin = sin(phi*mm);

    cosn=1.0;
    sinn=0.0;

    for (nk=1; nk<n; nk++) {

      k = dk*nk;

      lcos = cosn;
      lsin = sinn;
      cosn = lcos*dcos - lsin*dsin;
      sinn = lsin*dcos + lcos*dsin;

      accum_cos[mm][nk] += norm * mass * mcos * table[mm][nk] * sinn ;
#ifdef CHECK_OF
      if (accum_cos[mm][nk]*accum_cos[mm][nk] > 1.0e6) {
	cout.form("Process %d: mm=%d  nk=%d  r=%f  phi=%f  z=%f  accum_cos=%f\n", myid, mm, nk, r, phi, z, table[mm][nk]*table[mm][nk]);
      }
#endif // CHECK_OF
      if (mm>0) {
	accum_sin[mm][nk] += norm * mass * msin * table[mm][nk] * sinn ;
#ifdef CHECK_OF
	if (accum_sin[mm][nk]*accum_sin[mm][nk] > 1.0e6) {
	  cout.form("Process %d: mm=%d  nk=%d  r=%f  phi=%f  z=%f  accum_sin=%f\n", myid, mm, nk, r, phi, z, norm*mass*msin*table[mm][nk]*table[mm][nk]);
	}
#endif // CHECK_OF
      }
    }
  }

}


void CylindricalSL::make_coefficients(void)
{
  int mm, nn, j;

  if (!MPIset) {
    MPIin = new double [(n-1)*NMAX*(MMAX+1)];
    MPIout = new double [(n-1)*NMAX*(MMAX+1)];
    MPIset = TRUE;
  }
  
#ifdef MPE_PROFILE
  MPE_Log_event(7, myid, "b_distrib_c");
#endif

  for (mm=0; mm<=MMAX; mm++)
    for (nn=1; nn<n; nn++)
      for (j=1; j<=NMAX; j++) 
	MPIin[mm*(n-1)*NMAX + (nn-1)*NMAX + (j-1)] = accum_cos[mm][nn][j];
  
  MPI_Allreduce ( MPIin, MPIout, (n-1)*NMAX*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (mm=0; mm<=MMAX; mm++)
    for (nn=1; nn<n; nn++)
      for (j=1; j<=NMAX; j++) 
	accum_cos[mm][nn][j] = MPIout[mm*(n-1)*NMAX + (nn-1)*NMAX + (j-1)];
  



  for (mm=1; mm<=MMAX; mm++)
    for (nn=1; nn<n; nn++)
      for (j=1; j<=NMAX; j++) 
	MPIin[mm*(n-1)*NMAX + (nn-1)*NMAX + (j-1)] = accum_sin[mm][nn][j];
  
  MPI_Allreduce ( MPIin, MPIout, (n-1)*NMAX*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

  for (mm=1; mm<=MMAX; mm++)
    for (nn=1; nn<n; nn++)
      for (j=1; j<=NMAX; j++) 
	accum_sin[mm][nn][j] = MPIout[mm*(n-1)*NMAX + (nn-1)*NMAX + (j-1)];
  
#ifdef MPE_PROFILE
  MPE_Log_event(8, myid, "e_distrib_c");
#endif

  coefs_made = TRUE;
}

  
void CylindricalSL::accumulated_eval(double r, double z, double phi,
				     double& p, double& fr, double& fz, 
				     double &fp)
{
  if (!coefs_made) make_coefficients();

  fr = 0.0;
  fz = 0.0;
  fp = 0.0;
  p = 0.0;


  double k, ccos, ssin, znorm, facpc, facfc, facps, facfs;
  int nk, mm;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn, sinn;

  ortho->get_pot(table, r, 0, MMAX);
  ortho->get_force(tablef, r, 0, MMAX);

  for (mm=0; mm<=MMAX; mm++) {

    ccos = cos(phi*mm);
    ssin = sin(phi*mm);

    if (mm)
      znorm = 1.0/M_PI/ZMAX;
    else
      znorm = 0.5/M_PI/ZMAX;

    cosn = 1.0;
    sinn = 0.0;

    for (nk=1; nk<n; nk++) {

      k = dk*nk;

      lcos = cosn;
      lsin = sinn;
      cosn = lcos*dcos - lsin*dsin;
      sinn = lsin*dcos + lcos*dsin;
    
      facpc = accum_cos[mm][nk] * table[mm][nk]; 
      facfc = accum_cos[mm][nk] * tablef[mm][nk]; 
      if (mm>0) {
	facps = accum_sin[mm][nk] * table[mm][nk]; 
	facfs = accum_sin[mm][nk] * tablef[mm][nk]; 
      }
      else
	facps = facfs = 0.0;

      p +=   (facpc*ccos + facps*ssin) * sinn * znorm;
      fr += -(facfc*ccos + facfs*ssin) * sinn * znorm;
      fz += -(facpc*ccos + facps*ssin) * k*cosn * znorm;
      fp +=  (facpc*ssin - facps*ccos) * sinn * mm * znorm;
    }

  }

}

double CylindricalSL::accumulated_dens_eval(double r, double z, double phi)
{
  if (!coefs_made) make_coefficients();

  double ans = 0.0;

  double k, ccos, ssin, znorm, facpc, facps;
  int nk, mm;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn, sinn;

  ortho->get_pot(table, r, 0, MMAX);

  for (mm=0; mm<=MMAX; mm++) {

    ccos = cos(phi*mm);
    ssin = sin(phi*mm);

    if (mm)
      znorm = 1.0/M_PI/ZMAX;
    else
      znorm = 0.5/M_PI/ZMAX;

    cosn = 1.0;
    sinn = 0.0;

    for (nk=1; nk<n; nk++) {

      k = dk*nk;

      lcos = cosn;
      lsin = sinn;
      cosn = lcos*dcos - lsin*dsin;
      sinn = lsin*dcos + lcos*dsin;
    
      facpc = accum_cos[mm][nk] * table[mm][nk]; 
      if (mm>0)
	facps = accum_sin[mm][nk] * table[mm][nk]; 
      else
	facps = 0.0;

      ans += (facpc*ccos + facps*ssin) * sinn * znorm;
    }

  }

  return 0.25*ans/M_PI;
}


#include <stdio.h>

void CylindricalSL::dump_coefs(FILE *fout)
{
  double znorm;

  for (int mm=0; mm<=MMAX; mm++) {

    if (mm)
      znorm = 1.0/M_PI/ZMAX;
    else
      znorm = 0.5/M_PI/ZMAX;

    for (int nk=1; nk<n; nk++) {

      fprintf(fout, "%4d %4d", mm, nk);

      for (int j=1; j<=NMAX; j++)
	fprintf(fout, " %13.4e", accum_cos[mm][nk][j]*znorm);
      fprintf(fout, "\n");

      if (mm) {

	fprintf(fout, "%4d %4d", mm, nk);

	for (int j=1; j<=NMAX; j++)
	  fprintf(fout, " %13.4e", accum_sin[mm][nk][j]*znorm);
	fprintf(fout, "\n");
      }

    }
  }
}

