// #define DEBUG 1
// #define DEBUG_NAN 1
// #define GHQL 1
// #define STURM 1
// #define EIGEN 1
// #define DEBUG_PCA

#include <iostream.h>
#include <iomanip.h>
#include <strstream.h>

#include <string>
#include <algorithm>

#include <interp.h>

#include "exp_thread.h"

				// Constants from expand.h
const int nthrds = 1;
double tnow = 0.0;
static pthread_mutex_t used_lock;


#include <numerical.h>
#include <gaussQ.h>
#include <EmpOrth9thd.h>

Vector Symmetric_Eigenvalues_MSRCH(Matrix& a, Matrix& ef, int M);


#undef TINY
#define TINY 1.0e-16


bool EmpCylSL::DENS=false;
bool EmpCylSL::SELECT=false;
bool EmpCylSL::CMAP=false;
bool EmpCylSL::logarithmic=false;
EmpCylSL::EmpModel EmpCylSL::mtype = Exponential;
int EmpCylSL::NUMX=64;
int EmpCylSL::NUMY=128;
int EmpCylSL::NOUT=12;
int EmpCylSL::NUMR=1000;
double EmpCylSL::RMIN=0.001;
double EmpCylSL::RMAX=20.0;
string EmpCylSL::CACHEFILE = ".eof.cache.file";
string EmpCylSL::TABLEFILE = ".eof.table.file";


EmpCylSL::EmpCylSL(void)
{
  NORDER=0;
  MPIset = false;
  MPIset_eof = false;
  coefs_made = false;
  eof_made = false;
  eof_recompute = true;

  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  SC = 0;
  SS = 0;

  ortho = 0;
  model = 0;

  accum_cos = 0;
  accum_sin = 0;
  accum_cos2 = 0;
  accum_sin2 = 0;

  accum_cos0 = 0;
  accum_sin0 = 0;

  mpi_double_buf1 = 0;
  mpi_double_buf2 = 0;
  mpi_double_buf3 = 0;
}

EmpCylSL::~EmpCylSL(void)
{
  delete ortho;
  delete model;

  if (SC) {

    for (int m=0; m<=MMAX; m++) {
      delete [] potC[m];
      delete [] rforceC[m];
      delete [] zforceC[m];
      if (DENS) delete [] densC[m];
      if (m) {
	delete [] potS[m];
	delete [] rforceS[m];
	delete [] zforceS[m];
	if (DENS) delete [] densS[m];
      }
    }


    delete [] potC;
    delete [] rforceC;
    delete [] zforceC;
    if (DENS) delete [] densC;

    delete [] potS;
    delete [] rforceS;
    delete [] zforceS;
    if (DENS) delete [] densS;

    delete [] tpot;
    delete [] trforce;
    delete [] tzforce;
    if (DENS) delete [] tdens;

    delete [] mpi_double_buf1;
    delete [] mpi_double_buf2;
    delete [] mpi_double_buf3;

    for (int nth=0; nth<nthrds; nth++) {

      for (int mm=0; mm<=MMAX; mm++) {
      
	for (int j1=1; j1<=NMAX*(LMAX-mm+1); j1++) {

	  delete [] (SC[nth][mm][j1] + 1);
	  if (mm) delete [] (SS[nth][mm][j1] + 1);
	}
	  
	delete [] (SC[nth][mm] + 1);
	if (mm) delete [] (SS[nth][mm] + 1);
      }
	
      delete [] SC[nth];
      delete [] SS[nth];
    }
    
    delete [] SC;
    delete [] SS;

    delete [] vc;
    delete [] vs;
    delete [] hold;
    delete [] var;

    delete [] cosm;
    delete [] sinm;
    delete [] legs;
    delete [] dlegs;

    delete [] table;
    delete [] facC;
    delete [] facS;
  }

  if (accum_cos) {

    for (int nth=0; nth<nthrds; nth++) {
      delete [] accum_cos[nth];
      delete [] accum_sin[nth];
      if (SELECT) {
	delete [] accum_cos2[nth];
	delete [] accum_sin2[nth];
      }
      
      delete [] accum_cos0[nth];
      delete [] accum_sin0[nth];
    }

    delete [] accum_cos;
    delete [] accum_sin;
    delete [] accum_cos2;
    delete [] accum_sin2;
    
    delete [] accum_cos0;
    delete [] accum_sin0;
  }
  
  if (MPIset) {
    delete [] MPIin;
    delete [] MPIout;
  }

  if (MPIset_eof) {
    delete [] MPIin_eof;
    delete [] MPIout_eof;
  }

}


EmpCylSL::EmpCylSL(int nmax, int lmax, int mmax, int nord, 
		   double ascale, double hscale)
{
  NMAX = nmax;
  MMAX = mmax;
  LMAX = lmax;
  NORDER = nord;

  ASCALE = ascale;
  HSCALE = hscale;
  pfac = 1.0/sqrt(ascale);
  ffac = pfac/ascale;
  dfac = ffac/ascale;

  SLGridSph::mpi = 1;		// Turn on MPI
  ortho = new SLGridSph(LMAX, NMAX, NUMR, RMIN, RMAX*0.99, make_sl(), 1, 1.0);

  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  SC = 0;
  SS = 0;

  MPIset = false;
  MPIset_eof = false;
  coefs_made = false;
  eof_made = false;
  eof_recompute = true;

  accum_cos = 0;
  accum_sin = 0;
  accum_cos2 = 0;
  accum_sin2 = 0;

  accum_cos0 = 0;
  accum_sin0 = 0;

  mpi_double_buf1 = 0;
  mpi_double_buf2 = 0;
  mpi_double_buf3 = 0;
}


void EmpCylSL::reset(int numr, int lmax, int mmax, int nord, 
		     double ascale, double hscale)
{
  NMAX = numr;
  MMAX = mmax;
  LMAX = lmax;
  NORDER = nord;

  ASCALE = ascale;
  HSCALE = hscale;
  pfac = 1.0/sqrt(ascale);
  ffac = pfac/ascale;
  dfac = ffac/ascale;

  SLGridSph::mpi = 1;		// Turn on MPI
  ortho = new SLGridSph(LMAX, NMAX, NUMR, RMIN, RMAX*0.99, make_sl(), 1, 1.0);

  SC = 0;
  SS = 0;

  MPIset = false;
  MPIset_eof = false;
  coefs_made = false;
  eof_made = false;
  eof_recompute = true;

  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  accum_cos = 0;
  accum_sin = 0;
  accum_cos2 = 0;
  accum_sin2 = 0;

  accum_cos0 = 0;
  accum_sin0 = 0;

  mpi_double_buf1 = 0;
  mpi_double_buf2 = 0;
  mpi_double_buf3 = 0;
}

/*
  Note that the produced by the following three routines
  are in dimensionless units
*/
double EmpCylSL::massR(double R)
{
  double ans=0.0, fac, arg;

  switch (mtype) {
  case Exponential:
    ans = 1.0 - (1.0 + R)*exp(-R); 
    break;
  case Gaussian:
    arg = 0.5*R*R;
    ans = 1.0 - exp(-arg);
    break;
  case Plummer:
    fac = R/(1.0+R);
    ans = pow(fac, 3.0);
    break;
  }

  return ans;
}

double EmpCylSL::densR(double R)
{
  double ans=0.0, fac, arg;

  switch (mtype) {
  case Exponential:
    ans = exp(-R)/(4.0*M_PI*R);
    break;
  case Gaussian:
    arg = 0.5*R*R;
    ans = exp(-arg)/(4.0*M_PI*R);
    break;
  case Plummer:
    fac = 1.0/(1.0+R);
    ans = 3.0*pow(fac, 4.0)/(4.0*M_PI);
    break;
  }

  return ans;
}

SphericalModelTable* EmpCylSL::make_sl()
{
  const int number = 10000;

  r =  vector<double>(number);
  d =  vector<double>(number);
  m =  vector<double>(number);
  p =  vector<double>(number);

  vector<double> mm(number);
  vector<double> pw(number);

				// ------------------------------------------
				// Make radial, density and mass array
				// ------------------------------------------
  double dr;
  if (logarithmic)
    dr = (log(RMAX) - log(RMIN))/(number - 1);
  else
    dr = (RMAX - RMIN)/(number - 1);

  for (int i=0; i<number; i++) {
    if (logarithmic)
      r[i] = RMIN*exp(dr*i);
    else
      r[i] = RMIN + dr*i;

    m[i] = massR(r[i]);
    d[i] = densR(r[i]);
  }

  mm[0] = 0.0;
  pw[0] = 0.0;
  for (int i=1; i<number; i++) {
    mm[i] = mm[i-1] +
      2.0*M_PI*(r[i-1]*r[i-1]*d[i-1] + r[i]*r[i]*d[i])*
      (r[i] - r[i-1]);
    pw[i] = pw[i-1] +
      2.0*M_PI*(r[i-1]*d[i-1] + r[i]*d[i])*(r[i] - r[i-1]);
  }

  for (int i=0; i<number; i++) 
    p[i] = -mm[i]/(r[i]+1.0e-10) - (pw[number-1] - pw[i]);

#ifdef DEBUG
  ostrstream outf;
  outf << "test_adddisk_sl." << myid << '\0';
  ofstream out(outf.str());
  for (int i=0; i<number; i++) {
    out 
      << setw(15) << r[i] 
      << setw(15) << d[i] 
      << setw(15) << m[i] 
      << setw(15) << p[i] 
      << setw(15) << mm[i] 
      << endl;
  }
  out.close();
#endif

  model = new SphericalModelTable(number, &r[0]-1, &d[0]-1, &m[0]-1, &p[0]-1); 

  return model;
}

void EmpCylSL::send_eof_grid()
{
  double *MPIbuf  = new double [MPIbufsz];

#ifdef DEBUG
  cerr << "Process " << myid << ": size=" << MPIbufsz << "\n";
#endif
				// Sum up over threads

  for (int nth=1; nth<nthrds; nth++) {

    for (int m=0; m<=MMAX; m++) {
      for (int v=0; v<rank3; v++) 
	accum_cos[0][m][v] += accum_cos[nth][m][v];
    }

    for (int m=1; m<=MMAX; m++) {
      for (int v=0; v<rank3; v++) 
	accum_sin[0][m][v] += accum_sin[nth][m][v];
    }
  }

				// Send to slaves

  if (myid==0) {

    for (int m=0; m<=MMAX; m++) {

				// Coefficients for current step

      for (int v=0; v<rank3; v++) mpi_double_buf3[v] = accum_cos[0][m][v];
      MPI_Bcast(mpi_double_buf3, rank3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Grids in X--Y

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = potC[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = rforceC[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = zforceC[m][v][ix][iy];
	
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      MPIbuf[ix*(NUMY+1) + iy] = densC[m][v][ix][iy];

	  MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

      }

    }

    for (int m=1; m<=MMAX; m++) {

				// Coefficients for current step

      for (int v=0; v<rank3; v++) mpi_double_buf3[v] = accum_sin[0][m][v];
      MPI_Bcast(mpi_double_buf3, rank3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      

				// Grids in X--Y

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = potS[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = rforceS[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = zforceS[m][v][ix][iy];
	
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      MPIbuf[ix*(NUMY+1) + iy] = densS[m][v][ix][iy];

	  MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	
      }

    }

  }
  else {
				// Get tables from Master

    for (int m=0; m<=MMAX; m++) {

				// Coefficients for current step

      MPI_Bcast(mpi_double_buf3, rank3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      for (int v=0; v<rank3; v++) accum_cos[0][m][v] = mpi_double_buf3[v];


				// Grids in X--Y

      for (int v=0; v<rank3; v++) {

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    potC[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    rforceC[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    zforceC[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	if (DENS) {

	  MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      densC[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
	}

      }
    }

    for (int m=1; m<=MMAX; m++) {

				// Coefficients for current step

      MPI_Bcast(mpi_double_buf3, rank3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      for (int v=0; v<rank3; v++) accum_sin[0][m][v] = mpi_double_buf3[v];


				// Grids in X--Y

      for (int v=0; v<rank3; v++) {

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    potS[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    rforceS[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    zforceS[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	if (DENS) {

	  MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      densS[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	}

      }
    }

  }

  delete [] MPIbuf;
  
}


int EmpCylSL::read_cache(void)
{
  setup_accumulation();

				// Master tries to read table
  int retcode;
  if (myid==0) retcode = cache_grid(0);
  MPI_Bcast(&retcode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (!retcode) return 0;
				// Send table to slave processes
  send_eof_grid();

  if (myid==0) 
    cerr << "EmpCylSL::read_cache: table forwarded to all processes\n";


  eof_made = true;
  coefs_made = true;
  eof_recompute = false;

  return 1;
}


int EmpCylSL::cache_grid(int readwrite)
{

  if (readwrite) {

    ofstream out(CACHEFILE.c_str());
    if (!out) {
      cerr << "EmpCylSL::cache_grid: error writing file" << endl;
      return 0;
    }

    const int one = 1;
    const int zero = 0;

    out.write(&MMAX, sizeof(int));
    out.write(&NUMX, sizeof(int));
    out.write(&NUMY, sizeof(int));
    out.write(&NMAX, sizeof(int));
    out.write(&NORDER, sizeof(int));
    if (DENS) out.write(&one, sizeof(int));
    else      out.write(&zero, sizeof(int));
    if (CMAP) out.write(&one, sizeof(int));
    else      out.write(&zero, sizeof(int));
    out.write(&RMIN, sizeof(double));
    out.write(&RMAX, sizeof(double));
    out.write(&ASCALE, sizeof(double));
    out.write(&HSCALE, sizeof(double));
    out.write(&tnow, sizeof(double));

				// Write table

    for (int m=0; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write(&potC[m][v][ix][iy], sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write(&rforceC[m][v][ix][iy], sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write(&zforceC[m][v][ix][iy], sizeof(double));
	  
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      out.write(&densC[m][v][ix][iy], sizeof(double));

	}
	
      }

    }

    for (int m=1; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write(&potS[m][v][ix][iy], sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write(&rforceS[m][v][ix][iy], sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write(&zforceS[m][v][ix][iy], sizeof(double));
	
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      out.write(&densS[m][v][ix][iy], sizeof(double));
	}
	
      }

    }

  }
  else {

    ifstream in(CACHEFILE.c_str());
    if (!in) {
      cerr << "EmpCylSL::cache_grid: error opening file" << endl;
      return 0;
    }

    int mmax, numx, numy, nmax, norder, tmp;
    bool cmap=false, dens=false;
    double rmin, rmax, ascl, hscl;

    in.read(&mmax, sizeof(int));
    in.read(&numx, sizeof(int));
    in.read(&numy, sizeof(int));
    in.read(&nmax, sizeof(int));
    in.read(&norder, sizeof(int));
    in.read(&tmp, sizeof(int)); if (tmp) dens = true;
    in.read(&tmp, sizeof(int)); if (tmp) cmap = true;
    in.read(&rmin, sizeof(double));
    in.read(&rmax, sizeof(double));
    in.read(&ascl, sizeof(double));
    in.read(&hscl, sizeof(double));

				// Spot check compatibility
    if ( (MMAX    != mmax   ) |
	 (NUMX    != numx   ) |
	 (NUMY    != numy   ) |
	 (NMAX    != nmax   ) |
	 (NORDER  != norder ) |
	 (DENS    != dens   ) |
	 (CMAP    != cmap   ) |
	 (fabs(rmin-RMIN)>1.0e-12 ) |
	 (fabs(rmax-RMAX)>1.0e-12 ) |
	 (fabs(ascl-ASCALE)>1.0e-12 ) |
	 (fabs(hscl-HSCALE)>1.0e-12 )
	 ) 
      {
	cout << "MMAX=" << MMAX << " mmax=" << mmax << endl;
	cout << "NUMX=" << NUMX << " numx=" << numx << endl;
	cout << "NUMY=" << NUMY << " numy=" << numy << endl;
	cout << "NMAX=" << NMAX << " nmax=" << nmax << endl;
	cout << "NORDER=" << NORDER << " norder=" << norder << endl;
	cout << "DENS=" << DENS << " dens=" << dens << endl;
	cout << "CMAP=" << CMAP << " cmap=" << cmap << endl;
	cout << "RMIN=" << RMIN << " rmin=" << rmin << endl;
	cout << "RMAX=" << RMAX << " rmax=" << rmax << endl;
	cout << "ASCALE=" << ASCALE << " ascale=" << ascl << endl;
	cout << "HSCALE=" << HSCALE << " hscale=" << hscl << endl;
	return 0;
      }
    
    double time;
    in.read(&time, sizeof(double));

				// Read table

    for (int m=0; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read(&potC[m][v][ix][iy], sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read(&rforceC[m][v][ix][iy], sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read(&zforceC[m][v][ix][iy], sizeof(double));
	  
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      in.read(&densC[m][v][ix][iy], sizeof(double));

	}
	
      }

    }

    for (int m=1; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read(&potS[m][v][ix][iy], sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read(&rforceS[m][v][ix][iy], sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read(&zforceS[m][v][ix][iy], sizeof(double));
	
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      in.read(&densS[m][v][ix][iy], sizeof(double));
	}
	
      }

    }

    cerr << "EmpCylSL::cache_grid: file read successfully" << endl;
  }

  return 1;
}

void EmpCylSL::receive_eof(int request_id, int MM)
{

  int type, icnt, off;
  int mm;

#ifdef DEBUG
  cerr << "master listening . . . \n";
#endif

  MPI_Recv(&type, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	   MPI_COMM_WORLD, &status);

  int current_source = status.MPI_SOURCE;

#ifdef DEBUG
  cerr << "master beginning to receive from " << current_source << " . . . \n";
#endif

  MPI_Recv(&mm, 1, MPI_INT, current_source, MPI_ANY_TAG, 
	   MPI_COMM_WORLD, &status);

#ifdef DEBUG	
  cerr << "master receiving from " << current_source << ": type=" << type 
       << "   M=" << mm << "\n";
#endif

				// Receive rest of data

  MPI_Recv(mpi_double_buf1, NORDER, MPI_DOUBLE, current_source,
	    13, MPI_COMM_WORLD, &status);

  for (int n=0; n<NORDER; n++) {
    MPI_Recv(&mpi_double_buf2[MPIbufsz*(MPItable*n+0)], 
	     MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+1, 
	     MPI_COMM_WORLD, &status);
    MPI_Recv(&mpi_double_buf2[MPIbufsz*(MPItable*n+1)], 
	     MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+2, 
	     MPI_COMM_WORLD, &status);
    MPI_Recv(&mpi_double_buf2[MPIbufsz*(MPItable*n+2)], 
	     MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+3, 
	     MPI_COMM_WORLD, &status);
    if (DENS)
      MPI_Recv(&mpi_double_buf2[MPIbufsz*(MPItable*n+3)], 
	       MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+4, 
	       MPI_COMM_WORLD, &status);
  }
  

				// Send slave new orders
  if (request_id >=0) {
    MPI_Send(&request_id, 1, MPI_INT, current_source, 1, MPI_COMM_WORLD);
    MPI_Send(&MM, 1, MPI_INT, current_source, 2, MPI_COMM_WORLD);
  }
  else {
    MPI_Send(&request_id, 1, MPI_INT, current_source, 1, MPI_COMM_WORLD);
  }
  

				// Read from buffers

  for (int n=0; n<NORDER; n++) {
    if (type)
      accum_cos[0][mm][n] = mpi_double_buf1[n];
    else
      accum_sin[0][mm][n] = mpi_double_buf1[n];

  }

  for (int n=0; n<NORDER; n++) {
  
    off = MPIbufsz*(MPItable*n+0);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	potC[mm][n][ix][iy]  = mpi_double_buf2[off+icnt++];
      else
	potS[mm][n][ix][iy]  = mpi_double_buf2[off+icnt++];
	


    off = MPIbufsz*(MPItable*n+1);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	rforceC[mm][n][ix][iy]  = mpi_double_buf2[off+icnt++];
      else
	rforceS[mm][n][ix][iy]  = mpi_double_buf2[off+icnt++];
	

    off = MPIbufsz*(MPItable*n+2);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	zforceC[mm][n][ix][iy]  = mpi_double_buf2[off+icnt++];
      else
	zforceS[mm][n][ix][iy]  = mpi_double_buf2[off+icnt++];
    

    if (DENS) {
      off = MPIbufsz*(MPItable*n+3);
      icnt = 0;
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++)
	  if (type)
	    densC[mm][n][ix][iy]  = mpi_double_buf2[off+icnt++];
	  else
	    densS[mm][n][ix][iy]  = mpi_double_buf2[off+icnt++];
    }
  }
  
#ifdef DEBUG
  cerr << "master finished receiving: type=" << type << "   M=" 
       << mm << "\n";
#endif

  return;

}

void EmpCylSL::compute_eof_grid(int request_id, int m)
{
#ifdef EIGEN
  static Vector sum(ev.getlow(), ev.gethigh());
  {
    cerr << "Node " << myid 
	 << ": m=" << m << " dimen=" << ev.getlength() << "\n";
    sum[ev.getlow()] = ev[ev.getlow()];
    for (int nk=ev.getlow()+1; nk<=ev.gethigh(); nk++) sum[nk] = sum[nk-1] +
							 ev[nk];
    for (int nk=ev.getlow(); nk<=ev.gethigh(); nk++) {
      cerr << "       " << m << "," << nk << ">    "
	   << ev[nk] << "    " << sum[nk] << "\n";
      if (sum[nk]/sum[ev.gethigh()] > 1.0 - 1.0e-5) break;
    }
  }
#endif

  static int mpi_in_progress = 0;

  //  Read in coefficient matrix or
  //  make grid if needed


				// Sin/cos normalization
  double x, y, r, z;
  double costh, fac1, fac2, dens, potl, potr, pott, fac3, fac4;

  int icnt, off;

  for (int v=0; v<NORDER; v++) {
    tpot[v].zero();
    trforce[v].zero();
    tzforce[v].zero();
    if (DENS) tdens[v].zero();
  }

  for (int ix=0; ix<=NUMX; ix++) {

    x = XMIN + dX*ix;
    r = xi_to_r(x);

    for (int iy=0; iy<=NUMY; iy++) {

      y = YMIN + dY*iy;
      z = y_to_z(y);

      double rr = sqrt(r*r + z*z) + 1.0e-18;

      ortho->get_pot(potd, rr/ASCALE);
      ortho->get_force(dpot, rr/ASCALE);
      if (DENS) ortho->get_dens(dend, rr/ASCALE);

      costh = z/rr;
      dlegendre_R(LMAX, costh, legs[0], dlegs[0]);
      
      for (int v=0; v<NORDER; v++) {

	for (int ir=1; ir<=NMAX; ir++) {

	  for (int l=m; l<=LMAX; l++) {

	    fac1 = sqrt((2.0*l+1.0)/(4.0*M_PI));

	    if (m==0) {
	      fac2 = fac1*legs[0][l][m];

	      dens = fac2*dend[l][ir] * dfac;
	      potl = fac2*potd[l][ir] * pfac;
	      potr = fac2*dpot[l][ir] * ffac;
	      pott = fac1*dlegs[0][l][m]*potd[l][ir] * pfac;

	    } else {

	      fac2 = M_SQRT2 * fac1 * exp(0.5*(lgamma(l-m+1) - lgamma(l+m+1)));
	      fac3 = fac2 * legs[0][l][m];
	      fac4 = fac2 * dlegs[0][l][m];
	      
	      dens = fac3*dend[l][ir] * dfac;
	      potl = fac3*potd[l][ir] * pfac;
	      potr = fac3*dpot[l][ir] * ffac;
	      pott = fac4*potd[l][ir] * pfac;
	    }
	    
	    int nn = ir + NMAX*(l-m);

	    tpot[v][ix][iy] +=  ef[v+1][nn] * potl;

	    trforce[v][ix][iy] += 
	      -ef[v+1][nn] * (potr*r/rr - pott*z*r/(rr*rr*rr));

	    tzforce[v][ix][iy] += 
	      -ef[v+1][nn] * (potr*z/rr + pott*r*r/(rr*rr*rr));

	    if (DENS) 
	      tdens[v][ix][iy] +=  ef[v+1][nn] * dens;
	  }
	}
      }
    }
  }

				// Wait for previous delivery?

  mpi_in_progress = 0;

				// Send stuff back to master
      
  MPI_Send(&request_id, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
  MPI_Send(&m, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
      
				// First send back part of accum

  for (int n=0; n<NORDER; n++) {
    mpi_double_buf1[n] = 0.0;
    for (int i=1; i<=NMAX*(LMAX-m+1); i++) {
      if (request_id)
	mpi_double_buf1[n] += ef[n+1][i]*accum_cos0[0][m][i];
      else
	mpi_double_buf1[n] += ef[n+1][i]*accum_sin0[0][m][i];
    }
  }

  MPI_Send(mpi_double_buf1, NORDER, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);

    
  for (int n=0; n<NORDER; n++) {

				// normalization factors
    if (DENS)
      tdens[n] *= 0.25/M_PI;

				// Potential

    off = MPIbufsz*(MPItable*n+0);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = tpot[n][ix][iy];
      
    MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	     13 + MPItable*n+1, MPI_COMM_WORLD);

				// R force

    off = MPIbufsz*(MPItable*n+1);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = trforce[n][ix][iy];
    
    MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	     13 + MPItable*n+2, MPI_COMM_WORLD);

				// Z force

    off = MPIbufsz*(MPItable*n+2);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = tzforce[n][ix][iy];
    
    MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	     13 + MPItable*n+3, MPI_COMM_WORLD);

				// Density

    if (DENS) {
      off = MPIbufsz*(MPItable*n+3);
      icnt = 0;
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++)
	  mpi_double_buf2[off + icnt++] = tdens[n][ix][iy];
    
      MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	       13 + MPItable*n+4, MPI_COMM_WORLD);

    }

  }

  mpi_in_progress = 1;

}


void EmpCylSL::setup_accumulation(void)
{
  if (eof_recompute) setup_eof();


  if (!accum_cos) {

    accum_cos = new Vector* [nthrds];
    accum_sin = new Vector* [nthrds];
    for (int nth=0; nth<nthrds;nth++) {
      accum_cos[nth] = new Vector [MMAX+1];
      accum_sin[nth] = new Vector [MMAX+1];
      if (SELECT) {
	accum_cos2[nth] = new Vector [MMAX+1];
	accum_sin2[nth] = new Vector [MMAX+1];
      }
    }
#ifdef DEBUG_NAN
    cerr << "Slave " << myid << ": tables allocated, MMAX=" << MMAX << "\n";
#endif // DEBUG_NAN
  }


  cylused = cylused1 = 0;
  cylmass = 0.0;

  for (int nth=0; nth<nthrds;nth++) {

    for (int m=0; m<=MMAX; m++) {
      accum_cos[nth][m].setsize(0, rank3-1);
      accum_cos[nth][m].zero();
      if (SELECT) {
	accum_cos2[nth][m].setsize(0, rank3-1);
	accum_cos2[nth][m].zero();
      }
      if (m>0) {
	accum_sin[nth][m].setsize(0, rank3-1);
	accum_sin[nth][m].zero();
	if (SELECT) {
	  accum_sin2[nth][m].setsize(0, rank3-1);
	  accum_sin2[nth][m].zero();
	}
      }
    }
  }
  coefs_made = false;


}

void EmpCylSL::setup_eof()
{
  if (!SC) {

    rank2 = NMAX*(LMAX+1);
    rank3 = NORDER;

    ev.setsize(1, rank2);
    ef.setsize(1, rank2, 1, rank2);

    Rtable = M_SQRT1_2 * RMAX;
    XMIN = r_to_xi(RMIN);
    XMAX = r_to_xi(Rtable);
    dX = (XMAX - XMIN)/NUMX;

    YMIN = z_to_y(-Rtable);
    YMAX = z_to_y( Rtable);
    dY = (YMAX - YMIN)/NUMY;

    potC = new Matrix* [MMAX+1];
    rforceC = new Matrix* [MMAX+1];
    zforceC = new Matrix* [MMAX+1];
    if (DENS) densC = new Matrix* [MMAX+1];

    potS = new Matrix* [MMAX+1];
    rforceS = new Matrix* [MMAX+1];
    zforceS = new Matrix* [MMAX+1];
    if (DENS) densS = new Matrix* [MMAX+1];

    for (int m=0; m<=MMAX; m++) {

      potC[m] = new Matrix [rank3];
      rforceC[m] = new Matrix [rank3];
      zforceC[m] = new Matrix [rank3];
      if (DENS) densC[m] = new Matrix [rank3];

      for (int v=0; v<rank3; v++) {
	potC[m][v].setsize(0, NUMX, 0, NUMY);
	rforceC[m][v].setsize(0, NUMX, 0, NUMY);
	zforceC[m][v].setsize(0, NUMX, 0, NUMY);
	if (DENS) densC[m][v].setsize(0, NUMX, 0, NUMY);
      }

    }


    for (int m=1; m<=MMAX; m++) {

      potS[m] = new Matrix [rank3];
      rforceS[m] = new Matrix [rank3];
      zforceS[m] = new Matrix [rank3];
      if (DENS) densS[m] = new Matrix [rank3];

      for (int v=0; v<rank3; v++) {
	potS[m][v].setsize(0, NUMX, 0, NUMY);
	rforceS[m][v].setsize(0, NUMX, 0, NUMY);
	zforceS[m][v].setsize(0, NUMX, 0, NUMY);
	if (DENS) densS[m][v].setsize(0, NUMX, 0, NUMY);
      }

    }

    tpot = new Matrix [NORDER];
    trforce = new Matrix [NORDER];
    tzforce = new Matrix [NORDER];
    if (DENS) tdens = new Matrix [NORDER];

    for (int n=0; n<NORDER; n++) {
      tpot[n].setsize(0, NUMX, 0, NUMY);
      trforce[n].setsize(0, NUMX, 0, NUMY);
      tzforce[n].setsize(0, NUMX, 0, NUMY);
      if (DENS) tdens[n].setsize(0, NUMX, 0, NUMY);
    }

				// Allocate storage for each thread
    accum_cos0 = new Vector* [nthrds];
    accum_sin0 = new Vector* [nthrds];

    SC = new double*** [nthrds];
    SS = new double*** [nthrds];

    for (int nth=0; nth<nthrds; nth++) {
      accum_cos0[nth] = new Vector [MMAX+1];
      accum_sin0[nth] = new Vector [MMAX+1];

      SC[nth] = new double** [MMAX+1];
      SS[nth] = new double** [MMAX+1];
    }

    vc = new Matrix [nthrds];
    vs = new Matrix [nthrds];
    hold = new Vector[nthrds];
    for (int i=0; i<nthrds; i++) {
      vc[i].setsize(0, max<int>(1,MMAX), 0, rank3-1);
      vs[i].setsize(0, max<int>(1,MMAX), 0, rank3-1);
      hold[i].setsize(0, rank3-1);
    }

    var = new Matrix[MMAX+1];
    for (int m=0; m<=MMAX; m++)
      var[m].setsize(1, NMAX*(LMAX-m+1), 1, NMAX*(LMAX-m+1));
      
    potd.setsize(0, LMAX, 1, NMAX);
    dpot.setsize(0, LMAX, 1, NMAX);
    dend.setsize(0, LMAX, 1, NMAX);

    cosm = new Vector [nthrds];
    sinm = new Vector [nthrds];
    legs = new Matrix [nthrds];
    dlegs = new Matrix [nthrds];
    for (int i=0; i<nthrds; i++) {
      cosm[i].setsize(0, LMAX);
      sinm[i].setsize(0, LMAX);
      legs[i].setsize(0, LMAX, 0, LMAX);
      dlegs[i].setsize(0, LMAX, 0, LMAX);
    }


    for (int nth=0; nth<nthrds; nth++) {

      for (int m=0; m<=MMAX; m++) {

	SC[nth][m] = new double* [NMAX*(LMAX-m+1)] - 1;
	if (m) SS[nth][m] = new double* [NMAX*(LMAX-m+1)] - 1;

	accum_cos0[nth][m].setsize(1, NMAX*(LMAX-m+1));

	if (m) accum_sin0[nth][m].setsize(1, NMAX*(LMAX-m+1));

	for (int i=1; i<=NMAX*(LMAX-m+1); i++) {

	  SC[nth][m][i] = new double [NMAX*(LMAX-m+1)] - 1;
	  if (m) SS[nth][m][i] = new double [NMAX*(LMAX-m+1)] - 1;

	}

      }
    
    }

    table = new Matrix [nthrds];
    facC = new Matrix [nthrds];
    facS = new Matrix [nthrds];
    for (int i=0; i<nthrds; i++) {
      table[i].setsize(0, LMAX, 1, NMAX);
      facC[i].setsize(1, NMAX, 0, LMAX);
      facS[i].setsize(1, NMAX, 0, LMAX);
    }

    MPIbufsz = (NUMX+1)*(NUMY+1);

    mpi_double_buf1 = new double [rank2*rank2];
    mpi_double_buf2 = new double [MPIbufsz*NORDER*MPItable];
    mpi_double_buf3 = new double [rank3];

  }

  for (int nth=0; nth<nthrds; nth++) {
    for (int m=0; m<=MMAX; m++)  {
      accum_cos0[nth][m].zero();
      if (m>0) accum_sin0[nth][m].zero();
      for (int i=1; i<=NMAX*(LMAX-m+1); i++)  {
	for (int j=1; j<=NMAX*(LMAX-m+1); j++)  {
	  SC[nth][m][i][j] = 0.0;
	  if (m>0) SS[nth][m][i][j] = 0.0;
	}
      }
    }
  }

  eof_made = false;
}


void EmpCylSL::accumulate_eof(double r, double z, double phi, double mass, 
			      int id)
{
  if (eof_made) {
#ifdef DEBUG
    cerr << "accumulate_eof: Process " << myid << ", Thread " 
	 << id << " calling setup_eof()\n";
#endif      
    setup_eof();
  }

  double rr = sqrt(r*r + z*z);

  if (rr>Rtable) return;

  double fac0 = 4.0*M_PI, ylm;

  ortho->get_pot(table[id], rr/ASCALE);
  double costh = z/(rr+1.0e-18);
  legendre_R(LMAX, costh, legs[id]);
  sinecosine_R(LMAX, phi, cosm[id], sinm[id]);

  int nn1, nn2;

  // *** m loop
  for (int m=0; m<=MMAX; m++) {

    // *** ir loop
    for (int ir=1; ir<=NMAX; ir++) {

      // *** l loop
      for (int l=m; l<=LMAX; l++) {

	ylm = sqrt((2.0*l+1.0)/(4.0*M_PI)) * pfac *
	  exp(0.5*(lgamma(l-m+1) - lgamma(l+m+1))) * legs[0][l][m];

	if (m==0) {

	  facC[id][ir][l-m] = ylm*table[id][l][ir];

	}
	else {

	  facC[id][ir][l-m] = ylm*table[id][l][ir]*cosm[id][m];
	  facS[id][ir][l-m] = ylm*table[id][l][ir]*sinm[id][m];

	}

      } // *** l loop

    } // *** ir loop

    for (int ir1=1; ir1<=NMAX; ir1++) {
      for (int l1=m; l1<=LMAX; l1++) {
	nn1 = ir1 + NMAX*(l1-m);

	if (m==0) {

	  accum_cos0[id][m][nn1] += facC[id][ir1][l1-m]*mass*fac0;

	  for (int ir2=1; ir2<=NMAX; ir2++) {
	    for (int l2=m; l2<=LMAX; l2++) {
	      nn2 = ir2 + NMAX*(l2-m);

	      SC[id][m][nn1][nn2] += 
		facC[id][ir1][l1-m]*facC[id][ir2][l2-m] * mass;
	    }
	  }

	} else {

	  accum_cos0[id][m][nn1] += facC[id][ir1][l1-m]*mass*fac0;
	  accum_sin0[id][m][nn1] += facS[id][ir1][l1-m]*mass*fac0;
	  
	  for (int ir2=1; ir2<=NMAX; ir2++) {
	    for (int l2=m; l2<=LMAX; l2++) {
	      nn2 = ir2 + NMAX*(l2-m);

	      SC[id][m][nn1][nn2] += 
		facC[id][ir1][l1-m]*facC[id][ir2][l2-m] * mass;

	      SS[id][m][nn1][nn2] += 
		facS[id][ir1][l1-m]*facS[id][ir2][l2-m] * mass;
	    }
	  }
	}

      } // *** l loop

    } // *** ir loop

  } // *** m loop
  
}

void EmpCylSL::make_eof(void)
{
  int icnt;
  double tmp;

  if (!MPIset_eof) {
    MPIin_eof  = new double [rank2*(rank2+1)/2];
    MPIout_eof = new double [rank2*(rank2+1)/2];
    MPIset_eof = true;
  }
  
  //
  //  Sum up over threads
  //

  for (int nth=1; nth<nthrds; nth++) {

    for (int mm=0; mm<=MMAX; mm++) {

				// Mean
      for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	accum_cos0[0][mm][i] += accum_cos0[nth][mm][i];

				// Covariance
      for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	for (int j=i; j<=NMAX*(LMAX-mm+1); j++)
	  SC[0][mm][i][j] += SC[nth][mm][i][j];
  
    }

    for (int mm=1; mm<=MMAX; mm++) {

				// Mean
      for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	accum_sin0[0][mm][i] += accum_sin0[nth][mm][i];

				// Covariance
      for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	for (int j=i; j<=NMAX*(LMAX-mm+1); j++)
	  SS[0][mm][i][j] += SS[nth][mm][i][j];
  
    }

  }

  //
  //  Distribute mean and covariance to all processes
  //

  for (int mm=0; mm<=MMAX; mm++) {

				// Mean
    icnt=0;
    for (int i=1; i<=NMAX*(LMAX-mm+1); i++) 
      MPIin_eof[icnt++] = accum_cos0[0][mm][i];

    MPI_Allreduce ( MPIin_eof, MPIout_eof, NMAX*(LMAX-mm+1),
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      
    icnt=0;
    for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
      accum_cos0[0][mm][i] = MPIout_eof[icnt++];


				// Covariance
    icnt=0;
    for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
      for (int j=i; j<=NMAX*(LMAX-mm+1); j++)
	MPIin_eof[icnt++] = SC[0][mm][i][j];
    
    MPI_Allreduce ( MPIin_eof, MPIout_eof, 
		    NMAX*(LMAX-mm+1)*(NMAX*(LMAX-mm+1)+1)/2,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    icnt=0;
    for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
      for (int j=i; j<=NMAX*(LMAX-mm+1); j++)
	SC[0][mm][i][j] = MPIout_eof[icnt++];
    
  }
  
  for (int mm=1; mm<=MMAX; mm++) {

				// Mean
    icnt=0;
    for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
      MPIin_eof[icnt++] = accum_sin0[0][mm][i];

    MPI_Allreduce ( MPIin_eof, MPIout_eof, NMAX*(LMAX-mm+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    icnt=0;
    for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
      accum_sin0[0][mm][i] = MPIout_eof[icnt++];


				// Covariance
    icnt=0;
    for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
      for (int j=i; j<=NMAX*(LMAX-mm+1); j++)
	MPIin_eof[icnt++] = SS[0][mm][i][j];
  
    MPI_Allreduce ( MPIin_eof, MPIout_eof, 
		    NMAX*(LMAX-mm+1)*(NMAX*(LMAX-mm+1)+1)/2,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    icnt=0;
    for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
      for (int j=i; j<=NMAX*(LMAX-mm+1); j++)
	SS[0][mm][i][j] = MPIout_eof[icnt++];
  }


  if (myid==0) {

    int slave = 0;
    int request_id = 1;		// Begin with cosine case
    int M;

    M = 0;			// Initial counters
    while (M<=MMAX) {
	
      // Send request to slave
      if (slave<numprocs-1) {
	  
	slave++;
	  
	MPI_Send(&request_id, 1, MPI_INT, slave, 1, MPI_COMM_WORLD);
	MPI_Send(&M,  1, MPI_INT, slave, 2, MPI_COMM_WORLD);
	  
	// Increment counters
	request_id++;
	if (request_id>1) {
	  M++;
	  request_id = 0;
	}
	
#ifdef DEBUG
	cerr << "master in make_eof: done waiting on Slave " << slave << "\n";
#endif
      }
	
				// If M>MMAX before processor queue exhausted,
				// exit loop and reap the slave data
      if (M>MMAX) break;

      if (slave == numprocs-1) {
	  
	//
	// <Wait and receive and send new request>
	//
	  
	receive_eof(request_id, M);
	  
	// Increment counters

	request_id++;
	if (request_id>1) {
	  M++;
	  request_id = 0;
	}
      }
    }
    
    //
    // <Wait for all slaves to return and flag to continue>
    //
      
				// Dispatch resting slaves
    if (slave < numprocs-1) {
      request_id = -1;		// request_id < 0 means continue
      for (int s=numprocs-1; s>slave; s--) {
	MPI_Send(&request_id, 1, MPI_INT, s, 1, MPI_COMM_WORLD);
      }
    }

				// Get data from working slaves
    while (slave) {
      receive_eof(-1,0);
      slave--;
    }
      
  } else {

    int M, request_id;

    while (1) {

#ifdef DEBUG
      cerr << "Slave " << myid << ": listening . . . \n";
#endif


      MPI_Recv(&request_id, 1, 
	       MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
      
				// Done!
      if (request_id<0) break;

      MPI_Recv(&M, 1, 
	       MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
	

#ifdef DEBUG
      cerr << "Slave " << myid << ": received orders type="
	   << request_id << "  M=" << M << "\n";
#endif

      if (request_id) {

				// Complete symmetric part
    

	for (int i=1; i<NMAX*(LMAX-M+1); i++) {
	  for (int j=i; j<=NMAX*(LMAX-M+1); j++)
	    var[M][i][j] = SC[0][M][i][j];
	}

	for (int i=1; i<NMAX*(LMAX-M+1); i++) {
	  for (int j=i+1; j<=NMAX*(LMAX-M+1); j++) {
	    var[M][j][i] = SC[0][M][i][j];
	  }
	}
    
	double maxV = 0.0;
	for (int i=1; i<=NMAX*(LMAX-M+1); i++) {
	  for (int j=i; j<=NMAX*(LMAX-M+1); j++) {
	    tmp = fabs(var[M][i][j]);
	    if (tmp > maxV) maxV = tmp;
	  }
	}

	var[M] /= maxV;
    
    /*==========================*/
    /* Solve eigenvalue problem */
    /*==========================*/
    
#if defined(STURM)
	ev = Symmetric_Eigenvalues_MSRCH(var[M], ef, NORDER);

#elif defined(GHQL)
	ev = var[M].Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#else
	ev = var[M].Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#endif
      }
      else {

				// Complete symmetric part
    

	for (int i=1; i<NMAX*(LMAX-M+1); i++) {
	  for (int j=i; j<=NMAX*(LMAX-M+1); j++)
	    var[M][i][j] = SS[0][M][i][j];
	}

	for (int i=1; i<NMAX*(LMAX-M+1); i++) {
	  for (int j=i+1; j<=NMAX*(LMAX-M+1); j++) {
	    var[M][j][i] = SS[0][M][i][j];
	  }
	}
    
	double maxV = 0.0;
	for (int i=1; i<=NMAX*(LMAX-M+1); i++) {
	  for (int j=i; j<=NMAX*(LMAX-M+1); j++) {
	    tmp = fabs(var[M][i][j]);
	    if (tmp > maxV) maxV = tmp;
	  }
	}
	
	var[M] /= maxV;
    
    /*==========================*/
    /* Solve eigenvalue problem */
    /*==========================*/
    
#if defined(STURM)
	ev = Symmetric_Eigenvalues_MSRCH(var[M], ef, NORDER);

#elif defined(GHQL)
	ev = var[M].Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#else
	ev = var[M].Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#endif
      }

      compute_eof_grid(request_id, M);

    }

  }
				// Send grid to all processes
  send_eof_grid();

				// Cache table for restarts
				// (it would be nice to multithread or fork
				//  this call . . . )
  if (myid==0) cache_grid(1);

  eof_made = true;
  coefs_made = true;
  eof_recompute = false;
}


void EmpCylSL::accumulate(vector<Particle>& part)
{
  double r, phi, z, mass;

  setup_accumulation();

  vector<Particle>::iterator p;
  for (p=part.begin(); p!=part.end(); p++) {

    mass = p->mass;
    r = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
    phi = atan2(p->pos[1], p->pos[0]);
    z = p->pos[2];
    
    accumulate(r, z, phi, mass, 0);
  }

}

void EmpCylSL::accumulate(double r, double z, double phi, double mass, int id)
{

  if (coefs_made) {
#ifdef DEBUG
    cerr << "EmpCylSL::accumulate: Process " << myid << ", Thread " << id 
	 << ": calling setup_accumulation from accumulate\n";
#endif
    setup_accumulation();
  }

  if (eof_recompute) {
    accumulate_eof(r, z, phi, mass, id);
    return;
  }

  double rr = sqrt(r*r+z*z);
  if (rr>Rtable) return;

  double msin, mcos;
  int mm;

  double norm = -4.0*M_PI;
  
  if (SELECT) {
    pthread_mutex_lock(&used_lock);
    cylused1++;
    pthread_mutex_unlock(&used_lock);
  }

  get_pot(vc[id], vs[id], r, z);

  for (mm=0; mm<=MMAX; mm++) {

    mcos = cos(phi*mm);
    msin = sin(phi*mm);

    for (int nn=0; nn<rank3; nn++) 
      hold[id][nn] = norm * mass * mcos * vc[id][mm][nn];

    for (int nn=0; nn<rank3; nn++) 
      accum_cos[id][mm][nn] += hold[id][nn];

    if (SELECT) {
      for (int nn=0; nn<rank3; nn++) 
	accum_cos2[id][mm][nn] += hold[id][nn]*hold[id][nn]/mass;
    }
    if (mm>0) {
      for (int nn=0; nn<rank3; nn++) 
	hold[id][nn] = norm * mass * msin * vs[id][mm][nn];
      for (int nn=0; nn<rank3; nn++) 
	accum_sin[id][mm][nn] += hold[id][nn];
      if (SELECT) {
	for (int nn=0; nn<rank3; nn++) 
	  accum_sin2[id][mm][nn] += hold[id][nn]*hold[id][nn]/mass;
      }
    }
  }

  cylmass += mass;

}


void EmpCylSL::make_coefficients(void)
{
  if (eof_recompute) {
    make_eof();
    return;
  }

  int mm, nn;

  if (!MPIset) {
    MPIin  = new double [rank3*(MMAX+1)];
    MPIout = new double [rank3*(MMAX+1)];
    MPIset = true;
  }
  

				// Sum up over threads

  for (int nth=1; nth<nthrds; nth++) {

    for (mm=0; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++) {
	accum_cos[0][mm][nn] += accum_cos[nth][mm][nn];
	if (SELECT)
	  accum_cos2[0][mm][nn] += accum_cos2[nth][mm][nn];
      }


    for (mm=1; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++) {
	accum_sin[0][mm][nn] += accum_sin[nth][mm][nn];
	if (SELECT)
	  accum_sin2[0][mm][nn] += accum_sin2[nth][mm][nn];
      }

  }

				// Begin distribution loop


#ifdef MPE_PROFILE
  MPE_Log_event(7, myid, "b_distrib_c");
#endif

  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      MPIin[mm*rank3 + nn] = accum_cos[0][mm][nn];
  
  MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      accum_cos[0][mm][nn] = MPIout[mm*rank3 + nn];
  

  if (SELECT) {
    for (mm=0; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++)
	MPIin[mm*rank3 + nn] = accum_cos2[0][mm][nn];
  
    MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (mm=0; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++)
	accum_cos2[0][mm][nn] = MPIout[mm*rank3 + nn];
  }  


  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      MPIin[mm*rank3 + nn] = accum_sin[0][mm][nn];
  
  MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      accum_sin[0][mm][nn] = MPIout[mm*rank3 + nn];
  


  if (SELECT) {
    for (mm=1; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++)
	MPIin[mm*rank3 + nn] = accum_sin2[0][mm][nn];
  
    MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

    for (mm=1; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++)
	accum_sin2[0][mm][nn] = MPIout[mm*rank3 + nn];
  }
  
#ifdef MPE_PROFILE
  MPE_Log_event(8, myid, "e_distrib_c");
#endif

  if (SELECT) pca_hall();

  coefs_made = true;
}


void EmpCylSL::pca_hall(void)
{
  double sqr, fac;
  int mm, nn;

#ifdef DEBUG_PCA
  cerr << "Process " << myid << ": made it to pca_hall\n";
#endif

  
				// Need number of particles to compute variance
  MPI_Allreduce ( &cylused1, &cylused, 1,
		  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#ifdef DEBUG_PCA
  cerr << "Process " << myiud << ": using " << cylused << " particles\n";
#endif

  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++) {
      sqr = accum_cos[0][mm][nn]*accum_cos[0][mm][nn];
      //      fac = (accum_cos2[0][mm][nn] - sqr + 1.0e-10)/(sqr + 1.0e-18)/cylused;
      fac = (accum_cos2[0][mm][nn] - sqr + 1.0e-10)/(sqr + 1.0e-18);
#ifdef DEBUG_PCA
      if (myid==0) cerr << mm << ", " << nn << ", C:   "
			<< accum_cos2[0][mm][nn] << "  " 
			<< sqr << "  " << fac << '\n';
#endif
      accum_cos[0][mm][nn] *= 1.0/(1.0 + fac);
    }
  

  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++) {
      sqr = accum_sin[0][mm][nn]*accum_sin[0][mm][nn];
      //      fac = (accum_sin2[mm][nn] - sqr + 1.0e-10)/(sqr + 1.0e-18)/cylused;
      fac = (accum_sin2[0][mm][nn] - sqr + 1.0e-10)/(sqr + 1.0e-18);
#ifdef DEBUG_PCA
      if (myid==0) cerr << mm << ", " << nn << ", S:   "
			<< accum_sin2[0][mm][nn] << "  " 
			<< sqr << "  " << fac << '\n';
#endif
      accum_sin[0][mm][nn] *= 1.0/(1.0 + fac);
    }
  
#ifdef DEBUG_PCA
  cerr << "Process " << myid << ": exiting to pca_hall\n";
#endif
}


void EmpCylSL::accumulated_eval(double r, double z, double phi,
				double& p, double& fr, double& fz, 
				double &fp)
{
  //  if (!coefs_made) make_coefficients();

  fr = 0.0;
  fz = 0.0;
  fp = 0.0;
  p = 0.0;

  double rr = sqrt(r*r + z*z);
  if (rr>Rtable) return;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;
  
  if (ix < 0) ix = 0;
  if (iy < 0) iy = 0;
  
  if (ix >= NUMX) ix = NUMX-1;
  if (iy >= NUMY) iy = NUMY-1;
  
  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;
  
  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;
  
  double ccos, ssin=0.0, fac;
  int n, mm;
  
  for (mm=0; mm<=MMAX; mm++) {
    
    ccos = cos(phi*mm);
    ssin = sin(phi*mm);

    for (n=0; n<rank3; n++) {
      
      fac = accum_cos[0][mm][n] * ccos;
      
      p += fac *
	(
	 potC[mm][n][ix  ][iy  ] * c00 +
	 potC[mm][n][ix+1][iy  ] * c10 +
	 potC[mm][n][ix  ][iy+1] * c01 +
	 potC[mm][n][ix+1][iy+1] * c11 
	 );
      
      fr += fac *
	(
	 rforceC[mm][n][ix  ][iy  ] * c00 +
	 rforceC[mm][n][ix+1][iy  ] * c10 +
	 rforceC[mm][n][ix  ][iy+1] * c01 +
	 rforceC[mm][n][ix+1][iy+1] * c11
	 );
      
      fz += fac *
	(
	 zforceC[mm][n][ix  ][iy  ] * c00 +
	 zforceC[mm][n][ix+1][iy  ] * c10 +
	 zforceC[mm][n][ix  ][iy+1] * c01 +
	 zforceC[mm][n][ix+1][iy+1] * c11 
	 );
      
      fac = accum_cos[0][mm][n] * ssin;
      
      fp += fac * mm *
	(
	 potC[mm][n][ix  ][iy  ] * c00 +
	 potC[mm][n][ix+1][iy  ] * c10 +
	 potC[mm][n][ix  ][iy+1] * c01 +
	 potC[mm][n][ix+1][iy+1] * c11 
	 );
      
      
      if (mm) {
	
	fac = accum_sin[0][mm][n] * ssin;
	
	p += fac *
	  (
	   potS[mm][n][ix  ][iy  ] * c00 +
	   potS[mm][n][ix+1][iy  ] * c10 +
	   potS[mm][n][ix  ][iy+1] * c01 +
	   potS[mm][n][ix+1][iy+1] * c11 
	   );
	
	fr += fac *
	  (
	   rforceS[mm][n][ix  ][iy  ] * c00 +
	   rforceS[mm][n][ix+1][iy  ] * c10 +
	   rforceS[mm][n][ix  ][iy+1] * c01 +
	   rforceS[mm][n][ix+1][iy+1] * c11
	   );
	
	fz += fac *
	  (
	   zforceS[mm][n][ix  ][iy  ] * c00 +
	   zforceS[mm][n][ix+1][iy  ] * c10 +
	   zforceS[mm][n][ix  ][iy+1] * c01 +
	   zforceS[mm][n][ix+1][iy+1] * c11 
	   );
	
	fac = -accum_sin[0][mm][n] * ccos;
	
	fp += fac * mm *
	  (
	   potS[mm][n][ix  ][iy  ] * c00 +
	   potS[mm][n][ix+1][iy  ] * c10 +
	   potS[mm][n][ix  ][iy+1] * c01 +
	   potS[mm][n][ix+1][iy+1] * c11 
	   );
	
      }
      
    }
    
  }
  
}


double EmpCylSL::accumulated_dens_eval(double r, double z, double phi)
{
  if (!DENS) return 0.0;

  if (!coefs_made) make_coefficients();

  double ans = 0.0;

  double rr = sqrt(r*r + z*z);

  if (rr > Rtable) return ans;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;

  if (ix < 0) ix = 0;
  if (iy < 0) iy = 0;
  
  if (ix >= NUMX) ix = NUMX-1;
  if (iy >= NUMY) iy = NUMY-1;

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;

  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;

  double ccos, ssin=0.0, fac;
  int n, mm;

  for (mm=0; mm<=MMAX; mm++) {

    ccos = cos(phi*mm);
    ssin = sin(phi*mm);

    for (n=0; n<rank3; n++) {

      fac = accum_cos[0][mm][n]*ccos;

      ans += fac *
	(
	 densC[mm][n][ix  ][iy  ] * c00 +
	 densC[mm][n][ix+1][iy  ] * c10 +
	 densC[mm][n][ix  ][iy+1] * c01 +
	 densC[mm][n][ix+1][iy+1] * c11 
	 );

      if (mm) {

	fac = accum_sin[0][mm][n]*ssin;

	ans += fac *
	  (
	   densS[mm][n][ix  ][iy  ] * c00 +
	   densS[mm][n][ix+1][iy  ] * c10 +
	   densS[mm][n][ix  ][iy+1] * c01 +
	   densS[mm][n][ix+1][iy+1] * c11 
	   );
      }

    }

  }

  return ans;
}


  
void EmpCylSL::get_pot(Matrix& Vc, Matrix& Vs, double r, double z)
{
  Vc.setsize(0, max(1,MMAX), 0, rank3-1);
  Vs.setsize(0, max(1,MMAX), 0, rank3-1);

  if (z > Rtable) z = Rtable;
  if (z <-Rtable) z = -Rtable;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;

  if (ix < 0) ix = 0;
  if (iy < 0) iy = 0;
  
  if (ix >= NUMX) ix = NUMX-1;
  if (iy >= NUMY) iy = NUMY-1;
    
  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;

  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;

  double fac = 1.0;

  for (int mm=0; mm<=MMAX; mm++) {
    
    for (int n=0; n<rank3; n++) {

      Vc[mm][n] = fac *
	(
	 potC[mm][n][ix  ][iy  ] * c00 +
	 potC[mm][n][ix+1][iy  ] * c10 +
	 potC[mm][n][ix  ][iy+1] * c01 +
	 potC[mm][n][ix+1][iy+1] * c11 
	 );

      if (mm) {

	Vs[mm][n] = fac *
	  (
	   potS[mm][n][ix  ][iy  ] * c00 +
	   potS[mm][n][ix+1][iy  ] * c10 +
	   potS[mm][n][ix  ][iy+1] * c01 +
	   potS[mm][n][ix+1][iy+1] * c11 
	   );
      }

    }

  }

}
    
  
void EmpCylSL::get_all(int mm, int nn, 
		       double r, double z, double phi,
		       double& p, double& d, 
		       double& fr, double& fz, double &fp)
{
  if (!coefs_made) make_coefficients();

  fr = 0.0;
  fz = 0.0;
  fp = 0.0;
  p = 0.0;
  d = 0.0;

  double rr = sqrt(r*r + z*z);

  if (rr>Rtable) {
    p = -cylmass/(rr+1.0e-16);
    fr = p*r/(rr+1.0e-16)/(rr+1.0e-16);
    fz = p*z/(rr+1.0e-16)/(rr+1.0e-16);

    return;
  }

  if (z > Rtable) z = Rtable;
  if (z <-Rtable) z = -Rtable;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;

  if (ix < 0) ix = 0;
  if (iy < 0) iy = 0;
  
  if (ix >= NUMX) ix = NUMX-1;
  if (iy >= NUMY) iy = NUMY-1;

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;

  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;

  double ccos, ssin;

  ccos = cos(phi*mm);
  ssin = sin(phi*mm);

  p += ccos *
    (
     potC[mm][nn][ix  ][iy  ] * c00 +
     potC[mm][nn][ix+1][iy  ] * c10 +
     potC[mm][nn][ix  ][iy+1] * c01 +
     potC[mm][nn][ix+1][iy+1] * c11 
	 );

  fr += ccos *
    (
     rforceC[mm][nn][ix  ][iy  ] * c00 +
     rforceC[mm][nn][ix+1][iy  ] * c10 +
     rforceC[mm][nn][ix  ][iy+1] * c01 +
     rforceC[mm][nn][ix+1][iy+1] * c11
	 );
  
  fz += ccos *
    (
     zforceC[mm][nn][ix  ][iy  ] * c00 +
     zforceC[mm][nn][ix+1][iy  ] * c10 +
     zforceC[mm][nn][ix  ][iy+1] * c01 +
     zforceC[mm][nn][ix+1][iy+1] * c11 
     );
  
  fp += ssin * mm *
    (
     potC[mm][nn][ix  ][iy  ] * c00 +
     potC[mm][nn][ix+1][iy  ] * c10 +
     potC[mm][nn][ix  ][iy+1] * c01 +
     potC[mm][nn][ix+1][iy+1] * c11 
     );
  
  if (DENS)
  d += ccos *
    (
     densC[mm][nn][ix  ][iy  ] * c00 +
     densC[mm][nn][ix+1][iy  ] * c10 +
     densC[mm][nn][ix  ][iy+1] * c01 +
     densC[mm][nn][ix+1][iy+1] * c11 
     );


  if (mm) {
    
    p += ssin *
      (
       potS[mm][nn][ix  ][iy  ] * c00 +
       potS[mm][nn][ix+1][iy  ] * c10 +
       potS[mm][nn][ix  ][iy+1] * c01 +
       potS[mm][nn][ix+1][iy+1] * c11 
       );

    fr += ssin *
      (
       rforceS[mm][nn][ix  ][iy  ] * c00 +
       rforceS[mm][nn][ix+1][iy  ] * c10 +
       rforceS[mm][nn][ix  ][iy+1] * c01 +
       rforceS[mm][nn][ix+1][iy+1] * c11
       );

    fz += ssin *
      (
       zforceS[mm][nn][ix  ][iy  ] * c00 +
       zforceS[mm][nn][ix+1][iy  ] * c10 +
       zforceS[mm][nn][ix  ][iy+1] * c01 +
	   zforceS[mm][nn][ix+1][iy+1] * c11 
	   );

    fp += -ccos * mm *
	  (
	   potS[mm][nn][ix  ][iy  ] * c00 +
	   potS[mm][nn][ix+1][iy  ] * c10 +
	   potS[mm][nn][ix  ][iy+1] * c01 +
	   potS[mm][nn][ix+1][iy+1] * c11 
	   );
      
    if (DENS)
    d += ssin *
      (
       densS[mm][nn][ix  ][iy  ] * c00 +
       densS[mm][nn][ix+1][iy  ] * c10 +
       densS[mm][nn][ix  ][iy+1] * c01 +
       densS[mm][nn][ix+1][iy+1] * c11 
       );

  }
  
}


void EmpCylSL::dump_coefs(ostream& out)
{
  double znorm = 1.0;

  out.setf(ios::scientific);

  for (int mm=0; mm<=MMAX; mm++) {

    for (int l=mm; l<=LMAX; l++) {

      out << setw(4) << mm << setw(4) << l << endl;

      for (int j=1; j<=NMAX; j++)
	out << " " << setw(15) << accum_cos0[0][mm][j + NMAX*(l-mm)]*znorm;
      out << endl;

      if (mm) {

	out << setw(4) << mm << setw(4) << l << endl;

	for (int j=1; j<=NMAX; j++)
	  out << " " << setw(15) << accum_sin0[0][mm][j + NMAX*(l-mm)]*znorm;
	out << endl;

      }

    }
  }
}


#include <coef.H>
static CoefHeader coefheader;
static CoefHeader2 coefheader2;

void EmpCylSL::dump_coefs_binary_last(ostream& out, double time)
{
  double p, znorm = 1.0;

  coefheader.time = time;
  coefheader.mmax = MMAX;
  coefheader.nord = NORDER;
  coefheader.nmax = NMAX;

  out.write(&coefheader, sizeof(CoefHeader));

  for (int mm=0; mm<=MMAX; mm++) {

    for (int l=mm; l<=LMAX; l++) {
      
      for (int j=1; j<=NMAX; j++)
	out.write(&(p=accum_cos0[0][mm][j+NMAX*(l-mm)]*znorm), sizeof(double));
      
      if (mm) {

	for (int j=1; j<=NMAX; j++)
	  out.write(&(p=accum_sin0[0][mm][j+NMAX*(l-mm)]*znorm), sizeof(double));

      }
      
    }
  }
}


void EmpCylSL::dump_coefs_binary_curr(ostream& out, double time)
{
  coefheader2.time = time;
  coefheader2.mmax = MMAX;
  coefheader2.nmax = rank3;

  out.write(&coefheader2, sizeof(CoefHeader2));

  for (int mm=0; mm<=MMAX; mm++) {

    for (int j=0; j<rank3; j++)
      out.write(&accum_cos[0][mm][j], sizeof(double));
    
    if (mm) {

      for (int j=0; j<rank3; j++)
	out.write(&accum_sin[0][mm][j], sizeof(double));
      
    }
  }
}


void EmpCylSL::dump_basis(const string& name, int step)
{
  static string labels [] = {"pot.", "fr.", "fz.", "dens."};
  static int numx = 60;
  static int numy = 60;
  
  double rmax = 0.33*Rtable;
  double r, dr = rmax/(numx-1);
  double z, dz = 2.0*rmax/(numy-1);

  float zz;

  char strbuf[128];
  double fac=1;
  int n, mm;
  ofstream** outC = new ofstream* [MPItable];
  ofstream** outS = new ofstream* [MPItable];
  
  for (mm=0; mm<=MMAX; mm++) {

    for (n=0; n<=min<int>(NOUT, rank3); n++) {

				// Make output streams
      for (int i=0; i<MPItable; i++) {

	ostrstream ins(strbuf, 128);
	ins << name << ".C." << labels[i] << mm << "." << n
	    << "." << step << '\0';
	
	outC[i] = new ofstream(strbuf);
	outC[i]->write(&numx, sizeof(int));
	outC[i]->write(&numy, sizeof(int));
	outC[i]->write(&(zz=  0.0), sizeof(float));
	outC[i]->write(&(zz= rmax), sizeof(float));
	outC[i]->write(&(zz=-rmax), sizeof(float));
	outC[i]->write(&(zz= rmax), sizeof(float));
      }

      if (mm) {

	for (int i=0; i<MPItable; i++) {

	  ostrstream ins(strbuf, 128);
	  ins << name << ".S." << labels[i] << mm << "." << n 
	      << "." << step << '\0';
	
	  outS[i] = new ofstream(strbuf);
	  outS[i]->write(&numx,       sizeof(int));
	  outS[i]->write(&numy,       sizeof(int));
	  outS[i]->write(&(zz=  0.0), sizeof(float));
	  outS[i]->write(&(zz= rmax), sizeof(float));
	  outS[i]->write(&(zz=-rmax), sizeof(float));
	  outS[i]->write(&(zz= rmax), sizeof(float));
	}

      }


				// Ok, write data

      for (int k=0; k<numy; k++) {

	z = -rmax + dz*k;

	for (int j=0; j<numx; j++) {
	  
	  r = dr*j;

	  double X = (r_to_xi(r) - XMIN)/dX;
	  double Y = (z_to_y(z)  - YMIN)/dY;

	  int ix = (int)X;
	  int iy = (int)Y;

	  if (ix < 0) ix = 0;
	  if (iy < 0) iy = 0;
  
	  if (ix >= NUMX) ix = NUMX-1;
	  if (iy >= NUMY) iy = NUMY-1;

	  double delx0 = (double)ix + 1.0 - X;
	  double dely0 = (double)iy + 1.0 - Y;
	  double delx1 = X - (double)ix;
	  double dely1 = Y - (double)iy;

	  double c00 = delx0*dely0;
	  double c10 = delx1*dely0;
	  double c01 = delx0*dely1;
	  double c11 = delx1*dely1;

	  zz = fac*(
	    potC[mm][n][ix  ][iy  ] * c00 +
	    potC[mm][n][ix+1][iy  ] * c10 +
	    potC[mm][n][ix  ][iy+1] * c01 +
	    potC[mm][n][ix+1][iy+1] * c11 );

	  outC[0]->write(&zz, sizeof(float));
	    
	  zz = fac*(
	    rforceC[mm][n][ix  ][iy  ] * c00 +
	    rforceC[mm][n][ix+1][iy  ] * c10 +
	    rforceC[mm][n][ix  ][iy+1] * c01 +
	    rforceC[mm][n][ix+1][iy+1] * c11 );

	  outC[1]->write(&zz, sizeof(float));

	  zz = fac*(
	    zforceC[mm][n][ix  ][iy  ] * c00 +
	    zforceC[mm][n][ix+1][iy  ] * c10 +
	    zforceC[mm][n][ix  ][iy+1] * c01 +
	    zforceC[mm][n][ix+1][iy+1] * c11 );

	  outC[2]->write(&zz, sizeof(float));

	  if (DENS) {
	    zz = fac*(
	      densC[mm][n][ix  ][iy  ] * c00 +
	      densC[mm][n][ix+1][iy  ] * c10 +
	      densC[mm][n][ix  ][iy+1] * c01 +
	      densC[mm][n][ix+1][iy+1] * c11 );

	    outC[3]->write(&zz, sizeof(float));
	  }

	  if (mm) {
      
	    zz = fac*(
	      potS[mm][n][ix  ][iy  ] * c00 +
	      potS[mm][n][ix+1][iy  ] * c10 +
	      potS[mm][n][ix  ][iy+1] * c01 +
	      potS[mm][n][ix+1][iy+1] * c11 );

	    outS[0]->write(&zz, sizeof(float));
	    
	    zz = fac*(
	      rforceS[mm][n][ix  ][iy  ] * c00 +
	      rforceS[mm][n][ix+1][iy  ] * c10 +
	      rforceS[mm][n][ix  ][iy+1] * c01 +
	      rforceS[mm][n][ix+1][iy+1] * c11 );

	    outS[1]->write(&zz, sizeof(float));

	    zz = fac*(
	      zforceS[mm][n][ix  ][iy  ] * c00 +
	      zforceS[mm][n][ix+1][iy  ] * c10 +
	      zforceS[mm][n][ix  ][iy+1] * c01 +
	      zforceS[mm][n][ix+1][iy+1] * c11 );
	    
	    outS[2]->write(&zz, sizeof(float));
	    
	    if (DENS) {
	      zz = fac*(
		densS[mm][n][ix  ][iy  ] * c00 +
		densS[mm][n][ix+1][iy  ] * c10 +
		densS[mm][n][ix  ][iy+1] * c01 +
		densS[mm][n][ix+1][iy+1] * c11 );
	      
	      outS[3]->write(&zz, sizeof(float));
	    }

	  }

	}
      }
      
				// Close and delete output streams
      for (int i=0; i<MPItable; i++) {
	outC[i]->close();
	delete outC[i];

	if (mm) {
	  outS[i]->close();
	  delete outS[i];
	}
      }

    }
  }

  delete [] outC;
  delete [] outS;


}

double EmpCylSL::r_to_xi(double r)
{
  if (CMAP) {
    if (r<0.0) {
      ostrstream msg;
      msg << "radius=" << r << " < 0! [mapped]";
      bomb(string(msg.str()));
    }
    return (r/ASCALE - 1.0)/(r/ASCALE + 1.0);
  } else {
    if (r<0.0)  {
      ostrstream msg;
      msg << "radius=" << r << " < 0!";
      bomb(string(msg.str()));
    }
    return r;
  }
}
    
double EmpCylSL::xi_to_r(double xi)
{
  if (CMAP) {
    if (xi<-1.0) bomb("xi < -1!");
    if (xi>=1.0) bomb("xi >= 1!");

    return (1.0 + xi)/(1.0 - xi) * ASCALE;
  } else {
    return xi;
  }

}

double EmpCylSL::d_xi_to_r(double xi)
{
  if (CMAP) {
    if (xi<-1.0) bomb("xi < -1!");
    if (xi>=1.0) bomb("xi >= 1!");

    return 0.5*(1.0 - xi)*(1.0 - xi) / ASCALE;
  } else {
    return 1.0;
  }
}

void EmpCylSL::bomb(string oops)
{
  cerr << "EmpCylSL: " << oops << endl; 
  MPI_Abort(MPI_COMM_WORLD, -1);
  exit(-1);
}

