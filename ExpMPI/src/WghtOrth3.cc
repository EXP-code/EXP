// #define DEBUG 1
// #define DEBUG_NAN 1

#include <iostream>
#include <iomanip>

#include <numerical.h>
#include <WghtOrth3.h>

#include <gaussQ.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#undef TINY
#define TINY 1.0e-16

#include <expand.h>

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
  initialized = false;
  MPIset = false;
  coefs_made = false;

  accum_cos = 0;
  accum_sin = 0;
  table = 0;
  tablef = 0;
}

CylindricalSL::~CylindricalSL(void)
{
  if (initialized) {
    delete ortho;
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

  initialized = false;
  MPIset = false;
  coefs_made = false;

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
  dk = M_PI/ZMAX;

				// Sin/cos normalization
  zn.setsize(0, n-1);
  zn.zero();
  zn += 1.0;
  zn[0] = 1.0/sqrt(2.0);

  //
  // Initialize grid
  //

#ifdef DEBUG
  cout << "Process " << myid 
       << ": MMAX=" << MMAX
       << "  NMAX=" << NMAX
       << "  NUMR=" << NUMR
       << "  NUMK=" << n
       << "  RMIN=" << RMIN
       << "  RMAX=" << RMAX
       << "  L=" << ZMAX
       << '\n';
#endif


  SLGridCyl::mpi = 1;		// Turn on MPI
  ortho = new SLGridCyl(MMAX, NMAX, NUMR, n, RMIN, RMAX, ZMAX);

  initialized = true;
}


double CylindricalSL::get_pot(int mm, int nn, int kk, double r, double z)
{
  if (nn > NMAX) bomb ("Order must be smaller than max, reset");
  if (kk > 2*n-1) bomb ("Order must be smaller than max, reset");
  if (initialized == false) initialize();

  if (fabs(z)>ZMAX) return 0.0;

  double ans;

  if (kk<n)
    ans = ortho->get_pot(r, mm, nn, kk) * cos(dk*kk*(z + ZMAX)) * zn[kk] /
      sqrt(ZMAX);
  else {
    kk -= n;
    ans = ortho->get_pot(r, mm, nn, kk) * sin(dk*kk*(z + ZMAX)) * zn[kk] /
      sqrt(ZMAX);
  }

  return ans;
}


double CylindricalSL::get_dens(int mm, int kk, int nn, double r, double z)
{
  if (nn > NMAX) bomb ("Order must be smaller than max, reset");
  if (kk > n-1) bomb ("Order must be smaller than max, reset");
  if (initialized == false) initialize();

  if (fabs(z)>ZMAX) return 0.0;

  double ans;

  if (kk<n)
    ans = ortho->get_dens(r, mm, nn, kk) * cos(dk*kk*(z + ZMAX)) * zn[kk] /
      sqrt(ZMAX);
  else {
    kk -= n;
    ans = ortho->get_dens(r, mm, nn, kk) * sin(dk*kk*(z + ZMAX)) * zn[kk] /
      sqrt(ZMAX);
  }

  return 0.25/M_PI * ans;
}


Matrix CylindricalSL::inner_product_dens(int m, double (*fct)(double, double))
{
  if (initialized == false) initialize();

  Vector fkc(0, n-1);
  Vector fks(0, n-1);
  Matrix fr(0, n, 1, NMAX);

  Matrix coefs(0, 2*n-1, 1, NMAX);
  coefs.zero();

  LegeQuad lwk(NINT);

  double xi, r, wghtr, z;
  int iz;


  for (int kr=1; kr<=NINT; kr++) {

    xi = 2.0*(lwk.knot(kr) - 0.5);
    r = ortho->xi_to_r(xi);
    wghtr = 1.0/ortho->d_xi_to_r(xi) * r * 2.0 * lwk.weight(kr);

    double dZ2 = 2.0*ZMAX/NINT;

    fkc.zero();
    fks.zero();

    for (int k=0; k<n; k++) {
      for (iz=0; iz<NINT; iz++) {
	z = dZ2*( (double)iz + 0.5) - ZMAX;
	fkc[k] += dZ2 * cos(dk*k*(z+ZMAX))*(*fct)(r, z);
	fks[k] += dZ2 * sin(dk*k*(z+ZMAX))*(*fct)(r, z);
      }
    }

    ortho->get_pot(fr, xi, m, 0);

    for (iz=0; iz<n; iz++) coefs[iz]   += fr[iz] * wghtr * fkc[iz] * zn[iz];
    for (iz=1; iz<n; iz++) coefs[iz+n] += fr[iz] * wghtr * fks[iz] * zn[iz];
  }

  //
  // Norm----|
  //         |
  //         v
  return -4.0*M_PI / sqrt(ZMAX) * coefs;
}

Matrix CylindricalSL::inner_product_pot(int m, double (*fct)(double, double))
{
  if (initialized == false) initialize();

  Vector fkc(0, n-1);
  Vector fks(0, n-1);
  Matrix fr(0, n, 1, NMAX);

  Matrix coefs(0, 2*n-1, 1, NMAX);
  coefs.zero();

  LegeQuad lwk(NINT);

  double xi, r, z, wghtr;
  int iz;


  for (int kr=1; kr<=NINT; kr++) {

    xi = 2.0*(lwk.knot(kr) - 0.5);
    r = ortho->xi_to_r(xi);
    wghtr = 1.0/ortho->d_xi_to_r(xi) * r * 2.0 * lwk.weight(kr);

    double dZ2 = 2.0*ZMAX/NINT;

    fkc.zero();
    fks.zero();

    for (int k=0; k<n; k++) {
      for (iz=0; iz<NINT; iz++) {
	z = dZ2*( (double)iz + 0.5) - ZMAX;
	fkc[k] += dZ2 * cos(dk*k*(z+ZMAX))*(*fct)(r, z);
	fks[k] += dZ2 * sin(dk*k*(z+ZMAX))*(*fct)(r, z);
      }
    }

    ortho->get_dens(fr, xi, m, 0);

    for (iz=0; iz<n; iz++) coefs[iz]   += fr[iz] * wghtr * fkc[iz] * zn[iz];
    for (iz=1; iz<n; iz++) coefs[iz+n] += fr[iz] * wghtr * fks[iz] * zn[iz];
  }

  //
  //     |---Norm
  //     |
  //     v
  return -coefs / sqrt(ZMAX);
}


double CylindricalSL::pot_eval(Matrix& coefs, int m, double r, double z)
{
  if (initialized == false) initialize();

  double k, ans=0.0;
  int nk;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn=1, sinn=0.0;

  Matrix tab;

  ortho->get_pot(tab, r, m);

  for (nk=0; nk<n; nk++) {

    k = dk*nk;

    ans += (coefs[nk]*cosn + coefs[nk+n]*sinn) * tab[nk] * zn[nk];
    
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  return ans/sqrt(ZMAX);
;
}


double CylindricalSL::dens_eval(Matrix& coefs, int m, double r, double z)
{
  if (initialized == false) initialize();

  double k, ans=0.0;
  int nk;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn=1, sinn=0.0;

  Matrix tab;

  ortho->get_dens(tab, r, m);

  for (nk=0; nk<n; nk++) {

    k = dk*nk;

    ans += (coefs[nk]*cosn + coefs[nk+n]*sinn) * tab[nk] * zn[nk];
    
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  return ans/sqrt(ZMAX) * 0.25/M_PI;
}



double CylindricalSL::r_force_eval(Matrix& coefs, int m, double r, double z)
{
  if (initialized == false) initialize();

  double k, ans=0.0;
  int nk;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn=1, sinn=0.0;

  Matrix tab;

  ortho->get_force(tab, r, m);

  for (nk=0; nk<n; nk++) {

    k = dk*nk;

    ans += (coefs[nk]*cosn + coefs[nk+n]*sinn) * tab[nk] * zn[nk];
    
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
    
  }

  return -ans/sqrt(ZMAX);
;
}


double CylindricalSL::z_force_eval(Matrix& coefs, int m, double r, double z)
{
  if (initialized == false) initialize();

  double k, ans=0.0;
  int nk;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn=1, sinn=0.0;

  Matrix tab;

  ortho->get_pot(tab, r, m);

  for (nk=0; nk<n; nk++) {

    k = dk*nk;

    ans += (-coefs[nk]*sinn + coefs[nk+n]*cosn) * tab[nk] * k * zn[nk];
    
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  return -ans/sqrt(ZMAX);
;
}


void CylindricalSL::force_eval(Matrix& coefs, int m, double r, double z,
				 double& fr, double& fz)
{
  if (initialized == false) initialize();

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

  for (nk=0; nk<n; nk++) {

    k = dk*nk;

    ansr += ( coefs[nk]*cosn + coefs[nk+n]*sinn) * tabf[nk] * zn[nk];
    ansz += (-coefs[nk]*sinn + coefs[nk+n]*cosn) * tabp[nk] * k * zn[nk];
    
    lcos = cosn;
    lsin = sinn;
    cosn = lcos*dcos - lsin*dsin;
    sinn = lsin*dcos + lcos*dsin;
  }

  fr = -ansr/sqrt(ZMAX);
  fz = -ansz/sqrt(ZMAX);

}

void CylindricalSL::setup_accumulation(void)
{
  if (initialized == false) initialize();

  if (!accum_cos) {
    accum_cos = new Matrix [MMAX+1];
    accum_sin = new Matrix [MMAX+1];
    table = new Matrix [MMAX+1];
    tablef = new Matrix [MMAX+1];
#ifdef DEBUG_NAN
    cerr.form("Slave %d: tables allocated, MMAX=%d\n", myid, MMAX);
#endif // DEBUG_NAN
  }

  for (int m=0; m<=MMAX; m++) {
    accum_cos[m].setsize(0, 2*n-1, 1, NMAX);
    accum_cos[m].zero();
    if (m>0) {
      accum_sin[m].setsize(0, 2*n-1, 1, NMAX);
      accum_sin[m].zero();
    }
  }

  coefs_made = false;

}

void check_vector_values(const Vector& v)
{
  for (int i=v.getlow(); i<=v.gethigh(); i++)
    if (isinf(v[i]) || isnan(v[i]))
      {
	cerr << "check_vector: Illegal value\n";
      }
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

    for (nk=0; nk<n; nk++) {

      k = dk*nk;

      accum_cos[mm][nk  ] += norm * mass * mcos * table[mm][nk] * cosn * zn[nk];
      accum_cos[mm][nk+n] += norm * mass * mcos * table[mm][nk] * sinn * zn[nk];
#ifdef CHECK_OF
      if (accum_cos[mm][nk]*accum_cos[mm][nk] > 1.0e6) {
	cout.form("Process %d: mm=%d  nk=%d  r=%f  phi=%f  z=%f  accum_cos=%f\n", myid, mm, nk, r, phi, z, table[mm][nk]*table[mm][nk]);
      }
#endif // CHECK_OF
#ifdef DEBUG_NAN
      check_vector_values(accum_cos[mm][nk]);
      check_vector_values(accum_cos[mm][nk+n]);
#endif // DEBUG_NAN
      if (mm>0) {
	accum_sin[mm][nk  ] += norm * mass * msin * table[mm][nk] * cosn * zn[nk];
	accum_sin[mm][nk+n] += norm * mass * msin * table[mm][nk] * sinn * zn[nk];
#ifdef DEBUG_NAN
	check_vector_values(accum_sin[mm][nk]);
	check_vector_values(accum_sin[mm][nk+n]);
#endif // DEBUG_NAN
#ifdef CHECK_OF
	if (accum_sin[mm][nk]*accum_sin[mm][nk] > 1.0e6) {
	  cout.form("Process %d: mm=%d  nk=%d  r=%f  phi=%f  z=%f  accum_sin=%f\n", myid, mm, nk, r, phi, z, norm*mass*msin*table[mm][nk]*table[mm][nk]);
	}
#endif // CHECK_OF
      }

      lcos = cosn;
      lsin = sinn;
      cosn = lcos*dcos - lsin*dsin;
      sinn = lsin*dcos + lcos*dsin;

    }
  }

}


void CylindricalSL::make_coefficients(void)
{
  int mm, nn, j;

  if (!MPIset) {
    MPIin  = new double [2*n*NMAX*(MMAX+1)];
    MPIout = new double [2*n*NMAX*(MMAX+1)];
    MPIset = true;
  }
  
#ifdef MPE_PROFILE
  MPE_Log_event(7, myid, "b_distrib_c");
#endif

  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<2*n; nn++)
      for (j=1; j<=NMAX; j++) 
	MPIin[mm*2*n*NMAX + nn*NMAX + (j-1)] = accum_cos[mm][nn][j];
  
  MPI_Allreduce ( MPIin, MPIout, 2*n*NMAX*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<2*n; nn++)
      for (j=1; j<=NMAX; j++) 
	accum_cos[mm][nn][j] = MPIout[mm*2*n*NMAX + nn*NMAX + (j-1)];
  



  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<2*n; nn++)
      for (j=1; j<=NMAX; j++) 
	MPIin[mm*2*n*NMAX + nn*NMAX + (j-1)] = accum_sin[mm][nn][j];
  
  MPI_Allreduce ( MPIin, MPIout, 2*n*NMAX*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<2*n; nn++)
      for (j=1; j<=NMAX; j++) 
	accum_sin[mm][nn][j] = MPIout[mm*2*n*NMAX + nn*NMAX + (j-1)];
  
#ifdef MPE_PROFILE
  MPE_Log_event(8, myid, "e_distrib_c");
#endif

  coefs_made = true;
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


  double k, ccos, ssin, znorm;
  double facpcC, facfcC, facpsC, facfsC;
  double facpcS, facfcS, facpsS, facfsS;
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

    for (nk=0; nk<n; nk++) {

      k = dk*nk;

      facpcC = accum_cos[mm][nk  ] * table[mm][nk]; 
      facfcC = accum_cos[mm][nk  ] * tablef[mm][nk]; 
      facpcS = accum_cos[mm][nk+n] * table[mm][nk]; 
      facfcS = accum_cos[mm][nk+n] * tablef[mm][nk]; 
      if (mm>0) {
	facpsC = accum_sin[mm][nk  ] * table[mm][nk]; 
	facfsC = accum_sin[mm][nk  ] * tablef[mm][nk]; 
	facpsS = accum_sin[mm][nk+n] * table[mm][nk]; 
	facfsS = accum_sin[mm][nk+n] * tablef[mm][nk]; 
      }
      else
	facpsC = facfsC = facpsS = facfsS = 0.0;

#ifdef DEBUG_NAN
      if (
	  isnan(facfcC) || isinf(facfcC) ||
	  isnan(facfcS) || isinf(facfcS) ||
	  isnan(facfsC) || isinf(facfsC) ||
	  isnan(facfsS) || isinf(facfsS)
	  )
	{
	  cerr << "accumulated_force: invalid value\n";
	}
#endif


      p +=   (
	       (facpcC*ccos + facpsC*ssin) * cosn
	      +(facpcS*ccos + facpsS*ssin) * sinn
	      ) * znorm * zn[nk];
      fr += (
	     -(facfcC*ccos + facfsC*ssin) * cosn
	     -(facfcS*ccos + facfsS*ssin) * sinn
	     ) * znorm * zn[nk];
      fz += (
	      (facpcC*ccos + facpsC*ssin) * sinn
	     -(facpcS*ccos + facpsS*ssin) * cosn
	      ) * k * znorm * zn[nk];

      fp +=  (
	       (facpcC*ssin - facpsC*ccos) * cosn
	      +(facpcS*ssin - facpsS*ccos) * sinn
	       ) * mm * znorm * zn[nk];

      lcos = cosn;
      lsin = sinn;
      cosn = lcos*dcos - lsin*dsin;
      sinn = lsin*dcos + lcos*dsin;
    }

  }

}

double CylindricalSL::accumulated_dens_eval(double r, double z, double phi)
{
  if (!coefs_made) make_coefficients();

  double ans = 0.0;

  double k, ccos, ssin, znorm, facpcC, facpsC, facpcS, facpsS;
  int nk, mm;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn, sinn;

  ortho->get_dens(table, r, 0, MMAX);

  for (mm=0; mm<=MMAX; mm++) {

    ccos = cos(phi*mm);
    ssin = sin(phi*mm);

    if (mm)
      znorm = 1.0/M_PI/ZMAX;
    else
      znorm = 0.5/M_PI/ZMAX;

    cosn = 1.0;
    sinn = 0.0;

    for (nk=0; nk<n; nk++) {

      k = dk*nk;

      facpcC = accum_cos[mm][nk  ] * table[mm][nk]; 
      facpcS = accum_cos[mm][nk+n] * table[mm][nk]; 
      if (mm>0) {
	facpsC = accum_sin[mm][nk  ] * table[mm][nk]; 
	facpsS = accum_sin[mm][nk+n] * table[mm][nk]; 
      }
      else
	facpsC = facpsS = 0.0;

      ans += (
	       (facpcC*ccos + facpsC*ssin) * cosn
	      +(facpcS*ccos + facpsS*ssin) * sinn
	       ) * znorm * zn[nk];

      lcos = cosn;
      lsin = sinn;
      cosn = lcos*dcos - lsin*dsin;
      sinn = lsin*dcos + lcos*dsin;
    }

  }

  return 0.25*ans/M_PI;
}


double CylindricalSL::accumulated_potential_energy(void)
{
  if (!coefs_made) make_coefficients();

  double znorm, ans = 0.0;

  int nk, mm;

  for (mm=0; mm<=MMAX; mm++) {

    if (mm)
      znorm = 1.0/M_PI/ZMAX;
    else
      znorm = 0.5/M_PI/ZMAX;

    for (nk=0; nk<n; nk++) {
      
      ans += (
	      accum_cos[mm][nk  ]*accum_cos[mm][nk  ] +
	      accum_cos[mm][nk+n]*accum_cos[mm][nk+n]
	      ) * zn[nk]*zn[nk]*znorm;

      if (mm>0) {
	ans += (
		accum_sin[mm][nk  ]*accum_sin[mm][nk  ] +
		accum_sin[mm][nk+n]*accum_sin[mm][nk+n]
		) * zn[nk]*zn[nk]*znorm;
      }

    }

  }

  return -0.25*ans/M_PI * ans;
}


double CylindricalSL::accumulated_vertical_mass_eval(double r, double z, double phi)
{
  if (!coefs_made) make_coefficients();

  double ans = 0.0;

  double k, ccos, ssin, znorm, facpcC, facpsC, facpcS, facpsS;
  int nk, mm;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn, sinn, signnk;

  ortho->get_dens(table, r, 0, MMAX);

  for (mm=0; mm<=MMAX; mm++) {

    ccos = cos(phi*mm);
    ssin = sin(phi*mm);

    if (mm)
      znorm = 1.0/M_PI/ZMAX;
    else
      znorm = 0.5/M_PI/ZMAX;

    cosn = 1.0;
    sinn = 0.0;

    signnk = 1.0;

    for (nk=0; nk<n; nk++) {

      k = dk*nk;

      facpcC = accum_cos[mm][nk  ] * table[mm][nk]; 
      facpcS = accum_cos[mm][nk+n] * table[mm][nk]; 
      if (mm>0) {
	facpsC = accum_sin[mm][nk  ] * table[mm][nk]; 
	facpsS = accum_sin[mm][nk+n] * table[mm][nk]; 
      }
      else
	facpsC = facpsS = 0.0;

      if (nk==0)
	ans += (facpcC*ccos + facpsC*ssin)*(z + ZMAX) * znorm * zn[nk];
      else
	ans += ZMAX/M_PI/nk *
	  (  (facpcC*ccos + facpsC*ssin) * sinn
	    -(facpcS*ccos + facpsS*ssin) * (cosn - signnk)
	     ) * znorm * zn[nk];

      lcos = cosn;
      lsin = sinn;
      cosn = lcos*dcos - lsin*dsin;
      sinn = lsin*dcos + lcos*dsin;

      signnk *= -1.0;
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

    for (int nk=0; nk<n; nk++) {

      fprintf(fout, "%4d %4d\n", mm, nk);

      for (int j=1; j<=NMAX; j++)
	fprintf(fout, " %13.4e", accum_cos[mm][nk  ][j]*znorm*zn[nk]);
      fprintf(fout, "\n");

      for (int j=1; j<=NMAX; j++)
	fprintf(fout, " %13.4e", accum_cos[mm][nk+n][j]*znorm*zn[nk]);
      fprintf(fout, "\n");

      if (mm) {

	fprintf(fout, "%4d %4d\n", mm, nk);

	for (int j=1; j<=NMAX; j++)
	  fprintf(fout, " %13.4e", accum_sin[mm][nk  ][j]*znorm*zn[nk]);
	fprintf(fout, "\n");

	for (int j=1; j<=NMAX; j++)
	  fprintf(fout, " %13.4e", accum_sin[mm][nk+n][j]*znorm*zn[nk]);
	fprintf(fout, "\n");
      }

    }
  }
}


void CylindricalSL::read_coefs(istream *in)
{
  double znorm;
  int mm1, nk1;

  for (int mm=0; mm<=MMAX; mm++) {

    for (int nk=0; nk<n; nk++) {

      *in >> mm1;
      *in >> nk1;

      if (mm != mm1)
	bomb("problem reading coefficient dump . . . m values");
      if (nk != nk1)
	bomb("problem reading coefficient dump . . . k values");

      if (mm)
	znorm = 1.0/M_PI/ZMAX;
      else
	znorm = 0.5/M_PI/ZMAX;

      for (int j=1; j<=NMAX; j++)
	*in >> accum_cos[mm][nk  ][j];
      
      for (int j=1; j<=NMAX; j++)
	*in >> accum_cos[mm][nk+n][j];

      accum_cos[mm][nk  ] /= znorm * zn[nk];
      accum_cos[mm][nk+n] /= znorm * zn[nk];

      if (mm) {

	*in >> mm1;
	*in >> nk1;

	if (mm != mm1)
	  bomb("problem reading coefficient dump . . . m values");
	if (nk != nk1)
	  bomb("problem reading coefficient dump . . . k values");

	for (int j=1; j<=NMAX; j++)
	  *in >> accum_sin[mm][nk  ][j];

	for (int j=1; j<=NMAX; j++)
	  *in >> accum_sin[mm][nk+n][j];

	accum_sin[mm][nk  ] /= znorm * zn[nk];
	accum_sin[mm][nk+n] /= znorm * zn[nk];

      }

    }
  }

  coefs_made = true;

}


void CylindricalSL::dump_coefs_binary(FILE *fout, double time)
{
  double znorm, p;

  coefheader.time = time;
  coefheader.mmax = MMAX;
  coefheader.nord = n;
  coefheader.nmax = NMAX;

  fwrite(&coefheader, sizeof(CoefHeader), 1, fout);

  for (int mm=0; mm<=MMAX; mm++) {

    if (mm)
      znorm = 1.0/M_PI/ZMAX;
    else
      znorm = 0.5/M_PI/ZMAX;

    for (int nk=0; nk<n; nk++) {

      for (int j=1; j<=NMAX; j++)
	fwrite(&(p=accum_cos[mm][nk  ][j]*znorm*zn[nk]), sizeof(double), 1, fout);
      for (int j=1; j<=NMAX; j++)
	fwrite(&(p=accum_cos[mm][nk+n][j]*znorm*zn[nk]), sizeof(double), 1, fout);

      if (mm) {

	for (int j=1; j<=NMAX; j++)
	  fwrite(&(p=accum_sin[mm][nk  ][j]*znorm*zn[nk]), sizeof(double), 1, fout);
	for (int j=1; j<=NMAX; j++)
	  fwrite(&(p=accum_sin[mm][nk+n][j]*znorm*zn[nk]), sizeof(double), 1, fout);
      }

    }
  }
}


void CylindricalSL::read_coefs_binary(istream *in, double* time)
{
  *time = -999.0;		// Default error signal

  double znorm;

  int tst = in->rdstate();
  in->read((char *)&coefheader, sizeof(CoefHeader));
  tst = in->rdstate();
  
  if (MMAX != coefheader.mmax)
    bomb("MMAX mismatch reading binary coefficient file");
  if (n    != coefheader.nord)
    bomb("NFFT mismatch reading binary coefficient file");
  if (NMAX != coefheader.nmax)
    bomb("NMAX mismatch reading binary coefficient file");


  for (int mm=0; mm<=MMAX; mm++) {

    if (in->good()) 
      cout << "Read ok at m=" << mm << endl;
    else {
      int tst = in->rdstate();
      if (tst == (int)ios::goodbit)
	  cout << "Stream is really ok\n";
      if (tst == (int)ios::eofbit)
	cout << "Stream at end of file\n";
      if (tst == (int)ios::failbit)
	cout << "Stream input failed\n";
      if (tst == (int)ios::badbit)
	cout << "Stream is in unusable state\n";
    }

    for (int nk=0; nk<n; nk++) {

      if (mm)
	znorm = 1.0/M_PI/ZMAX;
      else
	znorm = 0.5/M_PI/ZMAX;

      for (int j=1; j<=NMAX; j++)
	in->read((char *)&accum_cos[mm][nk  ][j], sizeof(double));

      for (int j=1; j<=NMAX; j++)
	in->read((char *)&accum_cos[mm][nk+n][j], sizeof(double));

      accum_cos[mm][nk  ] /= znorm * zn[nk];
      accum_cos[mm][nk+n] /= znorm * zn[nk];

      if (mm) {

	for (int j=1; j<=NMAX; j++)
	  in->read((char *)&accum_sin[mm][nk  ][j], sizeof(double));

	for (int j=1; j<=NMAX; j++)
	  in->read((char *)&accum_sin[mm][nk+n][j], sizeof(double));

	accum_sin[mm][nk  ] /= znorm * zn[nk];
	accum_sin[mm][nk+n] /= znorm * zn[nk];

      }

    }
  }

  if (!*in) return;

  coefs_made = true;

  *time = coefheader.time;

}


