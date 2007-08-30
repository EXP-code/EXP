// #define DEBUG 1
// #define DEBUG_NAN 1
#define GHQL 1
// #define STURM 1
// #define EIGEN 1
// #define DEBUG_PCA

#include <iostream>
#include <iomanip>
#include <sstream>

#include <string>
#include <algorithm>

#include <interp.h>

#include "exp_thread.h"

#ifndef STANDALONE
#include "expand.h"
#else  
				// Constants from expand.h
extern int nthrds;
extern double tnow;
extern unsigned multistep;
#endif

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
bool EmpCylSL::enforce_limits=false;
EmpCylSL::EmpModel EmpCylSL::mtype = Exponential;
int EmpCylSL::NUMX=128;
int EmpCylSL::NUMY=64;
int EmpCylSL::NOUT=12;
int EmpCylSL::NUMR=2000;
double EmpCylSL::RMIN=0.001;
double EmpCylSL::RMAX=20.0;
double EmpCylSL::HFAC=0.2;
string EmpCylSL::CACHEFILE = ".eof.cache.file";
string EmpCylSL::TABLEFILE = ".eof.table.file";


EmpCylSL::EmpCylSL(void)
{
  NORDER=0;
  MPIset = false;
  MPIset_eof = false;
  coefs_made = vector<bool>(multistep+1, false);
  eof_made = false;

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

  cylmass = 0.0;
  cylmass_made = false;
  cylmass1 = vector<double>(nthrds);

  mpi_double_buf1 = 0;
  mpi_double_buf2 = 0;
  mpi_double_buf3 = 0;

  hallfile = "";
  hallcount = 0;
  hallfreq = 50;
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

    if (SELECT) {
      for (int nth=0; nth<nthrds; nth++) {
	delete [] accum_cos2[nth];
	delete [] accum_sin2[nth];
      }
    }
      
    delete [] accum_cos;
    delete [] accum_sin;

    delete [] accum_cos2;
    delete [] accum_sin2;
  }
  
  for (int M=0; M<accum_cosL.size(); M++) {
    
    for (int nth=0; nth<nthrds; nth++) {
      delete [] accum_cos0[M][nth];
      delete [] accum_sin0[M][nth];

      delete [] accum_cosL[M][nth];
      delete [] accum_sinL[M][nth];

      delete [] accum_cosN[M][nth];
      delete [] accum_sinN[M][nth];
    }

    delete [] accum_cos0[M];
    delete [] accum_sin0[M];

    delete [] accum_cosL[M];
    delete [] accum_sinL[M];

    delete [] accum_cosN[M];
    delete [] accum_sinN[M];
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
  coefs_made = vector<bool>(multistep+1, false);
  eof_made = false;

  accum_cos = 0;
  accum_sin = 0;

  accum_cos2 = 0;
  accum_sin2 = 0;

  cylmass = 0.0;
  cylmass1 = vector<double>(nthrds);
  cylmass_made = false;

  mpi_double_buf1 = 0;
  mpi_double_buf2 = 0;
  mpi_double_buf3 = 0;

  hallfile = "";
  hallcount = 0;
  hallfreq = 50;
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
  coefs_made = vector<bool>(multistep+1, false);
  eof_made = false;

  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  accum_cos = 0;
  accum_sin = 0;

  accum_cos2 = 0;
  accum_sin2 = 0;

  cylmass = 0.0;
  cylmass1 = vector<double>(nthrds);
  cylmass_made = false;

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
  ostringstream outf;
  outf << "test_adddisk_sl." << myid;
  ofstream out(outf.str().c_str());
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

				// Send to slaves
  if (myid==0) {
    
    for (int M=0; M<=multistep; M++) {
				// Coefficients for current step
      for (int m=0; m<=MMAX; m++) {
	for (int v=0; v<rank3; v++) mpi_double_buf3[v] = accum_cosN[M][0][m][v];
	MPI_Bcast(mpi_double_buf3, rank3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    }
				// Grids in X--Y
    for (int m=0; m<=MMAX; m++) {
      
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

    for (int M=0; M<=multistep; M++) {
				// Coefficients for current step
      for (int m=1; m<=MMAX; m++) {
	for (int v=0; v<rank3; v++) mpi_double_buf3[v] = accum_sinN[M][0][m][v];
	MPI_Bcast(mpi_double_buf3, rank3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    }
    
    for (int m=1; m<=MMAX; m++) {
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

  } else {

    for (int M=0; M<=multistep; M++) {
				// Coefficients for current step
      for (int m=0; m<=MMAX; m++) {
	MPI_Bcast(mpi_double_buf3, rank3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for (int v=0; v<rank3; v++) accum_cosN[M][0][m][v] = mpi_double_buf3[v];
      }
    }

				// Get tables from Master
    for (int m=0; m<=MMAX; m++) {

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

    for (int M=0; M<=multistep; M++) {

				// Coefficients for current step
      for (int m=1; m<=MMAX; m++) {
	
	MPI_Bcast(mpi_double_buf3, rank3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for (int v=0; v<rank3; v++) accum_sinN[M][0][m][v] = mpi_double_buf3[v];
      }
    }
    
    
    for (int m=1; m<=MMAX; m++) {
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
  coefs_made = vector<bool>(multistep+1, false);

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

    out.write((const char *)&MMAX, sizeof(int));
    out.write((const char *)&NUMX, sizeof(int));
    out.write((const char *)&NUMY, sizeof(int));
    out.write((const char *)&NMAX, sizeof(int));
    out.write((const char *)&NORDER, sizeof(int));
    if (DENS) out.write((const char *)&one, sizeof(int));
    else      out.write((const char *)&zero, sizeof(int));
    if (CMAP) out.write((const char *)&one, sizeof(int));
    else      out.write((const char *)&zero, sizeof(int));
    out.write((const char *)&RMIN, sizeof(double));
    out.write((const char *)&RMAX, sizeof(double));
    out.write((const char *)&ASCALE, sizeof(double));
    out.write((const char *)&HSCALE, sizeof(double));
    out.write((const char *)&cylmass, sizeof(double));
    out.write((const char *)&tnow, sizeof(double));

				// Write table

    for (int m=0; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&potC[m][v][ix][iy], sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&rforceC[m][v][ix][iy], sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&zforceC[m][v][ix][iy], sizeof(double));
	  
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      out.write((const char *)&densC[m][v][ix][iy], sizeof(double));

	}
	
      }

    }

    for (int m=1; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&potS[m][v][ix][iy], sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&rforceS[m][v][ix][iy], sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&zforceS[m][v][ix][iy], sizeof(double));
	
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      out.write((const char *)&densS[m][v][ix][iy], sizeof(double));
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

    in.read((char *)&mmax, sizeof(int));
    in.read((char *)&numx, sizeof(int));
    in.read((char *)&numy, sizeof(int));
    in.read((char *)&nmax, sizeof(int));
    in.read((char *)&norder, sizeof(int));
    in.read((char *)&tmp, sizeof(int)); if (tmp) dens = true;
    in.read((char *)&tmp, sizeof(int)); if (tmp) cmap = true;
    in.read((char *)&rmin, sizeof(double));
    in.read((char *)&rmax, sizeof(double));
    in.read((char *)&ascl, sizeof(double));
    in.read((char *)&hscl, sizeof(double));

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
    in.read((char *)&cylmass, sizeof(double));
    in.read((char *)&time, sizeof(double));

				// Read table

    for (int m=0; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&potC[m][v][ix][iy], sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&rforceC[m][v][ix][iy], sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&zforceC[m][v][ix][iy], sizeof(double));
	  
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      in.read((char *)&densC[m][v][ix][iy], sizeof(double));

	}
	
      }

    }

    for (int m=1; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&potS[m][v][ix][iy], sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&rforceS[m][v][ix][iy], sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&zforceS[m][v][ix][iy], sizeof(double));
	
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      in.read((char *)&densS[m][v][ix][iy], sizeof(double));
	}
	
      }

    }
    
    Rtable = M_SQRT1_2 * RMAX;
    XMIN = r_to_xi(RMIN*ASCALE);
    XMAX = r_to_xi(Rtable*ASCALE);
    dX = (XMAX - XMIN)/NUMX;
    
    YMIN = z_to_y(-Rtable*ASCALE);
    YMAX = z_to_y( Rtable*ASCALE);
    dY = (YMAX - YMIN)/NUMY;
    
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

  MPI_Recv(mpi_double_buf1, NORDER*(multistep+1), MPI_DOUBLE, current_source,
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

  for (int M=0; M<=multistep; M++) {
    for (int n=0; n<NORDER; n++) {
      if (type) {
	accum_cosN[M][0][mm][n] = mpi_double_buf1[M*NORDER+n];
      }
      else {
	accum_sinN[M][0][mm][n] = mpi_double_buf1[M*NORDER+n];
      }
    }
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
      
#ifdef DEBUG
  cerr << "Slave " << myid << ": sending type=" << request_id << "   M=" 
       << m << "\n";
#endif
				// First send back part of accum

  for (int M=0; M<=multistep; M++) {
    for (int n=0; n<NORDER; n++) {
      mpi_double_buf1[M*NORDER+n] = 0.0;
      for (int i=1; i<=NMAX*(LMAX-m+1); i++) {
	if (request_id)
	  mpi_double_buf1[M*NORDER+n] += -ef[n+1][i]*accum_cos0[M][0][m][i];
	else
	  mpi_double_buf1[M*NORDER+n] += -ef[n+1][i]*accum_sin0[M][0][m][i];
      }
    }
  }
  
  MPI_Send(mpi_double_buf1, NORDER*(multistep+1), MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
    
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
  if (!accum_cos) {

    accum_cos = new Vector [MMAX+1];
    accum_sin = new Vector [MMAX+1];

    if (SELECT) {
      accum_cos2 = new Vector* [nthrds];
      accum_sin2 = new Vector* [nthrds];

      for (int nth=0; nth<nthrds;nth++) {
	accum_cos2[nth] = new Vector [MMAX+1];
	accum_sin2[nth] = new Vector [MMAX+1];
      }
    }    

  }

  if (accum_cosL.size() == 0) {

    for (int M=0; M<=multistep; M++) {
      accum_cosL.push_back(new Vector* [nthrds]);
      accum_cosN.push_back(new Vector* [nthrds]);
      accum_sinL.push_back(new Vector* [nthrds]);
      accum_sinN.push_back(new Vector* [nthrds]);
      
      for (int nth=0; nth<nthrds; nth++) {
	accum_cosL[M][nth] = new Vector [MMAX+1];
	accum_cosN[M][nth] = new Vector [MMAX+1];
	accum_sinL[M][nth] = new Vector [MMAX+1];
	accum_sinN[M][nth] = new Vector [MMAX+1];
      }

      howmany1.push_back(vector<unsigned>(nthrds, 0));
      howmany.push_back(0);
    }
    
    differC = vector<Matrix>(multistep+1);
    differS = vector<Matrix>(multistep+1);

#ifdef DEBUG_NAN
    cerr << "Slave " << myid << ": tables allocated, MMAX=" << MMAX << "\n";
#endif // DEBUG_NAN
  }

  cylused = cylused1 = 0;
  cylmass = 0.0;
  cylmass_made = false;

  for (int M=0; M<=multistep; M++) {

    for (int nth=0; nth<nthrds; nth++) {

      for (int m=0; m<=MMAX; m++) {
	
	accum_cosN[M][nth][m].setsize(0, NORDER-1);
	accum_cosN[M][nth][m].zero();

	if (accum_cosL[M][nth][m].getlow()  !=0 ||
	    accum_cosL[M][nth][m].gethigh() != NORDER-1)
	  {
	    accum_cosL[M][nth][m].setsize(0, NORDER-1);
	    accum_cosL[M][nth][m].zero();
	  }

	if (m>0) {
	  accum_sinN[M][nth][m].setsize(0, NORDER-1);
	  accum_sinN[M][nth][m].zero();

	  if (accum_sinL[M][nth][m].getlow()  !=0 ||
	      accum_sinL[M][nth][m].gethigh() != NORDER-1)
	    {
	      accum_sinL[M][nth][m].setsize(0, NORDER-1);
	      accum_sinL[M][nth][m].zero();
	    }
	}
      }
    }
  }

  for (int m=0; m<=MMAX; m++) {
    accum_cos[m].setsize(0, NORDER-1);
    accum_cos[m].zero();
    if (m>0) {
      accum_sin[m].setsize(0, NORDER-1);
      accum_sin[m].zero();
    }
  }
  
  if (SELECT) {
    for (int nth=0; nth<nthrds; nth++) {
      for (int m=0; m<=MMAX; m++) {
	accum_cos2[nth][m].setsize(0, NORDER-1);
	accum_cos2[nth][m].zero();
	if (m>0) {
	  accum_sin2[nth][m].setsize(0, NORDER-1);
	  accum_sin2[nth][m].zero();
	}
      }
    }
  }
  
  coefs_made = vector<bool>(multistep+1, false);
}

void EmpCylSL::setup_accumulation(int M)
{
  howmany[M] = 0;

  for (int nth=0; nth<nthrds; nth++) {

    howmany1[M][nth] = 0;

    for (int m=0; m<=MMAX; m++) {
      
      accum_cosN[M][nth][m].setsize(0, NORDER-1);
      accum_cosN[M][nth][m].zero();

      if (m>0) {
	accum_sinN[M][nth][m].setsize(0, NORDER-1);
	accum_sinN[M][nth][m].zero();
	
      }
    }
  }
  
  coefs_made[M] = false;
}

void EmpCylSL::setup_eof()
{
  if (!SC) {

    rank2 = NMAX*(LMAX+1);
    rank3 = NORDER;
    
    ev.setsize(1, rank2);
    ef.setsize(1, rank2, 1, rank2);

    Rtable = M_SQRT1_2 * RMAX;
    XMIN = r_to_xi(RMIN*ASCALE);
    XMAX = r_to_xi(Rtable*ASCALE);
    dX = (XMAX - XMIN)/NUMX;

    YMIN = z_to_y(-Rtable*ASCALE);
    YMAX = z_to_y( Rtable*ASCALE);
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

    SC = new double*** [nthrds];
    SS = new double*** [nthrds];

    for (int nth=0; nth<nthrds; nth++) {
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

    mpi_double_buf1 = new double [NORDER*(multistep+1)];
    mpi_double_buf2 = new double [MPIbufsz*NORDER*MPItable];
    mpi_double_buf3 = new double [rank3];
  }

  if (accum_cos0.size()==0) {
    for (int M=0; M<=multistep; M++) {

      accum_cos0.push_back(new Vector* [nthrds]);
      accum_sin0.push_back(new Vector* [nthrds]);

      for (int nth=0; nth<nthrds; nth++) {

	accum_cos0[M][nth] = new Vector [MMAX+1];
	accum_sin0[M][nth] = new Vector [MMAX+1];

	for (int m=0; m<=MMAX; m++) {
	  accum_cos0[M][nth][m].setsize(1, NMAX*(LMAX-m+1));

	  if (m) accum_sin0[M][nth][m].setsize(1, NMAX*(LMAX-m+1));
	}
      }
    }
  }

  for (int nth=0; nth<nthrds; nth++) {
    for (int M=0; M<=multistep; M++)  {
      for (int m=0; m<=MMAX; m++)  {
	accum_cos0[M][nth][m].zero();
	if (m>0) accum_sin0[M][nth][m].zero();
      }
    }
  }

  for (int nth=0; nth<nthrds; nth++) {
      for (int m=0; m<=MMAX; m++)  {
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
			      int id, int mlevel)
{
  if (eof_made) {
#ifdef DEBUG
    cerr << "accumulate_eof: Process " << myid << ", Thread " 
	 << id << " calling setup_eof()\n";
#endif      
    setup_eof();
  }

  double rr = sqrt(r*r + z*z);

  if (rr/ASCALE>Rtable) return;

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

	// Only the l dependence is important here . . .

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
	  
	  accum_cos0[mlevel][id][m][nn1] += facC[id][ir1][l1-m]*mass*fac0;

	  for (int ir2=1; ir2<=NMAX; ir2++) {
	    for (int l2=m; l2<=LMAX; l2++) {
	      nn2 = ir2 + NMAX*(l2-m);

	      SC[id][m][nn1][nn2] += 
		facC[id][ir1][l1-m]*facC[id][ir2][l2-m] * mass;
	    }
	  }

	} else {

	  accum_cos0[mlevel][id][m][nn1] += facC[id][ir1][l1-m]*mass*fac0;
	  accum_sin0[mlevel][id][m][nn1] += facS[id][ir1][l1-m]*mass*fac0;
	  
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

    for (int M=0; M<=multistep; M++) {

      for (int mm=0; mm<=MMAX; mm++) {

				// Mean
	for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	  accum_cos0[M][0][mm][i] += accum_cos0[M][nth][mm][i];
      }
    }

				// Covariance
    for (int mm=0; mm<=MMAX; mm++) {

      for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	for (int j=i; j<=NMAX*(LMAX-mm+1); j++)
	  SC[0][mm][i][j] += SC[nth][mm][i][j];
  
    }

    for (int M=0; M<=multistep; M++) {
      for (int mm=1; mm<=MMAX; mm++) {
	
				// Mean
	for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	  accum_sin0[M][0][mm][i] += accum_sin0[M][nth][mm][i];

      }
    }

				// Covariance
    for (int mm=1; mm<=MMAX; mm++) {

      for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	for (int j=i; j<=NMAX*(LMAX-mm+1); j++)
	  SS[0][mm][i][j] += SS[nth][mm][i][j];
  
    }

  }

  //
  //  Distribute mean and covariance to all processes
  //

  for (int M=0; M<=multistep; M++) {

    for (int mm=0; mm<=MMAX; mm++) {
				// Mean
      icnt=0;
      for (int i=1; i<=NMAX*(LMAX-mm+1); i++) 
	MPIin_eof[icnt++] = accum_cos0[M][0][mm][i];

      MPI_Allreduce ( MPIin_eof, MPIout_eof, NMAX*(LMAX-mm+1),
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      
      icnt=0;
      for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	accum_cos0[M][0][mm][i] = MPIout_eof[icnt++];
    }

  }

  for (int mm=0; mm<=MMAX; mm++) {
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
  
  for (int M=0; M<=multistep; M++) {

    for (int mm=1; mm<=MMAX; mm++) {

				// Mean
      icnt=0;
      for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	MPIin_eof[icnt++] = accum_sin0[M][0][mm][i];

      MPI_Allreduce ( MPIin_eof, MPIout_eof, NMAX*(LMAX-mm+1),
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=1; i<=NMAX*(LMAX-mm+1); i++)
	accum_sin0[M][0][mm][i] = MPIout_eof[icnt++];
    }
  }

  for (int mm=1; mm<=MMAX; mm++) {
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
    
	try {
#if defined(STURM)
	  ev = Symmetric_Eigenvalues_MSRCH(var[M], ef, NORDER);
	  
#elif defined(GHQL)
	  ev = var[M].Symmetric_Eigenvalues_GHQL(ef);
	  ef = ef.Transpose();
#else
	  ev = var[M].Symmetric_Eigenvalues(ef);
	  ef = ef.Transpose();
#endif
	}
	catch (char const *msg) {
	  cerr << "Process " << myid << ": in eigenvalues problem, M=" << M
	       << ", request=" << request_id
	       << ", error string=" << msg << endl;

	  ostringstream sout;
	  sout << "eigenvalue.error." << myid;
	  ofstream out(sout.str().c_str());

	  out << "# M=" << M << " request=" << request_id << endl;
	  for (int i=1; i<=NMAX*(LMAX-M+1); i++) {
	    for (int j=i; j<=NMAX*(LMAX-M+1); j++)
	      out << setw(18) << var[M][i][j];
	    out << endl;
	  }

	  MPI_Abort(MPI_COMM_WORLD, -1);
	}

      } else {

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
    
	try {
#if defined(STURM)
	  ev = Symmetric_Eigenvalues_MSRCH(var[M], ef, NORDER);

#elif defined(GHQL)
	  ev = var[M].Symmetric_Eigenvalues_GHQL(ef);
	  ef = ef.Transpose();
#else
	  ev = var[M].Symmetric_Eigenvalues(ef);
	  ef = ef.Transpose();
#endif
	}
	catch (char const *msg) {
	  cerr << "Process " << myid << ": in eigenvalues problem, M=" << M
	       << ", request=" << request_id
	       << ", error string=" << msg << endl;

	  ostringstream sout;
	  sout << "eigenvalue.error." << myid;
	  ofstream out(sout.str().c_str());

	  out << "# M=" << M << " request_id=" << request_id << endl;
	  for (int i=1; i<=NMAX*(LMAX-M+1); i++) {
	    for (int j=i; j<=NMAX*(LMAX-M+1); j++)
	      out << setw(18) << var[M][i][j];
	    out << endl;
	  }

	  MPI_Abort(MPI_COMM_WORLD, -1);
	}

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
  coefs_made = vector<bool>(multistep+1, true);
}


void EmpCylSL::accumulate_eof(vector<Particle>& part)
{

  double r, phi, z, mass;

  setup_eof();

  vector<Particle>::iterator p;
  for (p=part.begin(); p!=part.end(); p++) {

    mass = p->mass;
    r = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
    phi = atan2(p->pos[1], p->pos[0]);
    z = p->pos[2];
    
    accumulate_eof(r, z, phi, mass, 0, p->level);
  }

}
  

void EmpCylSL::accumulate(vector<Particle>& part, int mlevel)
{
  double r, phi, z, mass;

  setup_accumulation();

  vector<Particle>::iterator p;
  for (p=part.begin(); p!=part.end(); p++) {

    mass = p->mass;
    r = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
    phi = atan2(p->pos[1], p->pos[0]);
    z = p->pos[2];
    
    accumulate(r, z, phi, mass, 0, mlevel);
  }

}
  

void EmpCylSL::accumulate(double r, double z, double phi, double mass, 
			  int id, int mlevel)
{

  if (coefs_made[mlevel]) {
    cerr << "EmpCylSL::accumulate: Process " << myid << ", Thread " << id 
	 << ": calling setup_accumulation from accumulate\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  double rr = sqrt(r*r+z*z);
  if (rr/ASCALE>Rtable) return;

  howmany1[mlevel][id]++;

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
      accum_cosN[mlevel][id][mm][nn] += hold[id][nn];

    if (SELECT) {
      for (int nn=0; nn<rank3; nn++) 
	accum_cos2[id][mm][nn] += hold[id][nn]*hold[id][nn]/mass;
    }
    if (mm>0) {
      for (int nn=0; nn<rank3; nn++) 
	hold[id][nn] = norm * mass * msin * vs[id][mm][nn];
      for (int nn=0; nn<rank3; nn++) 
	accum_sinN[mlevel][id][mm][nn] += hold[id][nn];
      if (SELECT) {
	for (int nn=0; nn<rank3; nn++) 
	  accum_sin2[id][mm][nn] += hold[id][nn]*hold[id][nn]/mass;
      }
    }
  }

  if (multistep==0 || mstep==Mstep)
    cylmass1[id] += mass;

}


void EmpCylSL::make_coefficients(int M)
{
  if (coefs_made[M]) return;

  int mm, nn;

  if (!MPIset) {
    MPIin  = new double [rank3*(MMAX+1)];
    MPIout = new double [rank3*(MMAX+1)];
    MPIset = true;
  }
  

				// Sum up over threads
  for (int nth=1; nth<nthrds; nth++) {

    howmany1[M][0] += howmany1[M][nth];

    for (mm=0; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++) {
	accum_cosN[M][0][mm][nn] += accum_cosN[M][nth][mm][nn];
      }


    for (mm=1; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++) {
	accum_sinN[M][0][mm][nn] += accum_sinN[M][nth][mm][nn];
      }
  }

				// Begin distribution loop


  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      MPIin[mm*rank3 + nn] = accum_cosN[M][0][mm][nn];
  
  MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      if (multistep)
	accum_cosN[M][0][mm][nn] = MPIout[mm*rank3 + nn];
      else
	accum_cos[mm][nn] = MPIout[mm*rank3 + nn];

  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      MPIin[mm*rank3 + nn] = accum_sinN[M][0][mm][nn];
  
  MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      if (multistep)
	accum_sinN[M][0][mm][nn] = MPIout[mm*rank3 + nn];
      else
	accum_sin[mm][nn] = MPIout[mm*rank3 + nn];
  
  coefs_made[M] = true;
}

void EmpCylSL::reset_mass(void)
{ 
  cylmass=0.0; 
  cylmass_made=false; 
  for (int n=0; n<nthrds; n++) cylmass1[n] = 0.0;
}

void EmpCylSL::make_coefficients(void)
{
  if (!cylmass_made) {
    MPI_Allreduce(&cylmass1[0], &cylmass, 1, MPI_DOUBLE, MPI_SUM, 
		  MPI_COMM_WORLD);
    cylmass_made = true;
  }
  
  if (coefs_made_all()) return;

  int mm, nn;

  if (!MPIset) {
    MPIin  = new double [rank3*(MMAX+1)];
    MPIout = new double [rank3*(MMAX+1)];
    MPIset = true;
  }
  

				// Sum up over threads
  for (int M=0; M<=multistep; M++) {

    if (coefs_made[M]) continue;

    for (int nth=1; nth<nthrds; nth++) {

      howmany1[M][0] += howmany1[M][nth];

      for (mm=0; mm<=MMAX; mm++)
	for (nn=0; nn<rank3; nn++) {
	  accum_cosN[M][0][mm][nn] += accum_cosN[M][nth][mm][nn];
	  if (SELECT && M==0)
	    accum_cos2[0][mm][nn] += accum_cos2[nth][mm][nn];
	}


      for (mm=1; mm<=MMAX; mm++)
	for (nn=0; nn<rank3; nn++) {
	  accum_sinN[M][0][mm][nn] += accum_sinN[M][nth][mm][nn];
	  if (SELECT && M==0)
	    accum_sin2[0][mm][nn] += accum_sin2[nth][mm][nn];
	}
      
    }
  }

				// Begin distribution loop

  for (int M=0; M<=multistep; M++) {

    if (coefs_made[M]) continue;

				// "howmany" is only used for debugging
    MPI_Allreduce ( &howmany1[M][0], &howmany[M], 1, MPI_UNSIGNED,
		    MPI_SUM, MPI_COMM_WORLD);

    for (mm=0; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++)
	MPIin[mm*rank3 + nn] = accum_cosN[M][0][mm][nn];
  
    MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (mm=0; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++)
	if (multistep)
	  accum_cosN[M][0][mm][nn] = MPIout[mm*rank3 + nn];
	else
	  accum_cos[mm][nn] = MPIout[mm*rank3 + nn];
  }
  

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


  for (int M=0; M<=multistep; M++) {
    
    if (coefs_made[M]) continue;

    for (mm=1; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++)
	MPIin[mm*rank3 + nn] = accum_sinN[M][0][mm][nn];
  
    MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

    for (mm=1; mm<=MMAX; mm++)
      for (nn=0; nn<rank3; nn++)
	if (multistep)
	  accum_sinN[M][0][mm][nn] = MPIout[mm*rank3 + nn];
	else
	  accum_sin[mm][nn] = MPIout[mm*rank3 + nn];
  }
  
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
  
  if (SELECT) pca_hall();

  coefs_made = vector<bool>(multistep+1, true);
}


void EmpCylSL::pca_hall(void)
{
  double sqr, var, fac, tot;
  int mm, nn;

#ifdef DEBUG_PCA
  cerr << "Process " << myid << ": made it to pca_hall\n";
#endif

  
				// Need number of particles to compute variance
  MPI_Allreduce ( &cylused1, &cylused, 1,
		  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#ifdef DEBUG_PCA
  cerr << "Process " << myid << ": using " << cylused << " particles\n";
#endif

  ofstream *hout = NULL;
  if (hallcount++%hallfreq==0 && myid==0 && hallfile.length()>0) {
    hout = new ofstream(hallfile.c_str(), ios::out | ios::app);
    if (!hout) {
      cerr << "Could not open <" << hallfile << "> for appending output\n";
    }
    *hout << "# Time = " << tnow << "  Step=" << hallcount << "\n";
    *hout << "#\n";
  }

  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++) {

      tot = 0.0;
      for (int M=0; M<=multistep; M++) tot += accum_cosN[M][0][mm][nn];

      sqr = tot*tot;
      // fac = (accum_cos2[0][mm][nn] - cylmass*sqr + 1.0e-10)/(sqr + 1.0e-18);
      var = accum_cos2[0][mm][nn] - cylmass*sqr;
      fac = sqr/(var/(cylused+1) + sqr + 1.0e-10);
      if (hout) *hout << mm << ", " << nn << ", C:   "
		      << setw(18) << accum_cos2[0][mm][nn] << "  " 
		      << setw(18) << sqr << "  " 
		      << setw(18) << var << "  " 
		      << setw(18) << fac << '\n';

      for (int M=0; M<=multistep; M++) {
	// accum_cosN[M][0][mm][nn] *= 1.0/(1.0 + fac);
	accum_cosN[M][0][mm][nn] *= fac;
      }
    }
  

  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++) {

      tot = 0.0;
      for (int M=0; M<=multistep; M++) tot += accum_sinN[M][0][mm][nn];

      sqr = tot*tot;
      // fac = (accum_sin2[0][mm][nn] - sqr + 1.0e-10)/(sqr + 1.0e-18);
      var = accum_sin2[0][mm][nn] - cylmass*sqr;
      fac = sqr/(var/(cylused+1) + sqr + 1.0e-10);
      if (hout) *hout << mm << ", " << nn << ", S:   "
		      << setw(18) << accum_sin2[0][mm][nn] << "  " 
		      << setw(18) << sqr << "  " 
		      << setw(18) << var << "  " 
		      << setw(18) << fac << '\n';
      for (int M=0; M<=multistep; M++) {
	// accum_sinN[M][0][mm][nn] *= 1.0/(1.0 + fac);
	accum_sinN[M][0][mm][nn] *= fac;
      }
    }
  
#ifdef DEBUG_PCA
  cerr << "Process " << myid << ": exiting to pca_hall\n";
#endif

  if (hout) {
    hout->close();
    delete hout;
  }

}


void EmpCylSL::accumulated_eval(double r, double z, double phi, 
				double &p0, double& p, 
				double& fr, double& fz, double &fp)
{
  if (!coefs_made_all()) make_coefficients();

  fr = 0.0;
  fz = 0.0;
  fp = 0.0;
  p = 0.0;

  double rr = sqrt(r*r + z*z);
  if (rr/ASCALE>Rtable) return;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;
  
  if (ix < 0) {
    ix = 0;
    if (enforce_limits) X = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    if (enforce_limits) Y = 0.0;
  }
  
  if (ix >= NUMX) {
    ix = NUMX-1;
    if (enforce_limits) X = NUMX;
  }
  if (iy >= NUMY) {
    iy = NUMY-1;
    if (enforce_limits) Y = NUMY;
  }

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;
  
  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;
  
  double ccos, ssin=0.0, fac;
  
  for (int mm=0; mm<=MMAX; mm++) {
    
    ccos = cos(phi*mm);
    ssin = sin(phi*mm);

    for (int n=0; n<rank3; n++) {
      
      fac = accum_cos[mm][n] * ccos;
      
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
      
      fac = accum_cos[mm][n] * ssin;
      
      fp += fac * mm *
	(
	 potC[mm][n][ix  ][iy  ] * c00 +
	 potC[mm][n][ix+1][iy  ] * c10 +
	 potC[mm][n][ix  ][iy+1] * c01 +
	 potC[mm][n][ix+1][iy+1] * c11 
	 );
      
      
      if (mm) {
	
	fac = accum_sin[mm][n] * ssin;
	
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
	
	fac = -accum_sin[mm][n] * ccos;
	
	fp += fac * mm *
	  (
	   potS[mm][n][ix  ][iy  ] * c00 +
	   potS[mm][n][ix+1][iy  ] * c10 +
	   potS[mm][n][ix  ][iy+1] * c01 +
	   potS[mm][n][ix+1][iy+1] * c11 
	   );
	
      }
      
    }
    
    if (mm==0) p0 = p;

  }
  
}


double EmpCylSL::accumulated_dens_eval(double r, double z, double phi, 
				       double& d0)
{
  if (!DENS) return 0.0;

  if (!coefs_made_all()) make_coefficients();

  double ans = 0.0;

  double rr = sqrt(r*r + z*z);

  if (rr/ASCALE > Rtable) return ans;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;

  if (ix < 0) {
    ix = 0;
    if (enforce_limits) X = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    if (enforce_limits) Y = 0.0;
  }
  
  if (ix >= NUMX) {
    ix = NUMX-1;
    if (enforce_limits) X = NUMX;
  }
  if (iy >= NUMY) {
    iy = NUMY-1;
    if (enforce_limits) Y = NUMY;
  }

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

      fac = accum_cos[mm][n]*ccos;

      ans += fac *
	(
	 densC[mm][n][ix  ][iy  ] * c00 +
	 densC[mm][n][ix+1][iy  ] * c10 +
	 densC[mm][n][ix  ][iy+1] * c01 +
	 densC[mm][n][ix+1][iy+1] * c11 
	 );

      if (mm) {

	fac = accum_sin[mm][n]*ssin;

	ans += fac *
	  (
	   densS[mm][n][ix  ][iy  ] * c00 +
	   densS[mm][n][ix+1][iy  ] * c10 +
	   densS[mm][n][ix  ][iy+1] * c01 +
	   densS[mm][n][ix+1][iy+1] * c11 
	   );
      }

    }

    if (mm==0) d0 = ans;

  }

  return ans;
}


  
void EmpCylSL::get_pot(Matrix& Vc, Matrix& Vs, double r, double z)
{
  Vc.setsize(0, max(1,MMAX), 0, rank3-1);
  Vs.setsize(0, max(1,MMAX), 0, rank3-1);

  if (z/ASCALE > Rtable) z =  Rtable*ASCALE;
  if (z/ASCALE <-Rtable) z = -Rtable*ASCALE;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;

  if (ix < 0) {
    ix = 0;
    if (enforce_limits) X = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    if (enforce_limits) Y = 0.0;
  }
  
  if (ix >= NUMX) {
    ix = NUMX-1;
    if (enforce_limits) X = NUMX;
  }
  if (iy >= NUMY) {
    iy = NUMY-1;
    if (enforce_limits) Y = NUMY;
  }

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
  if (!coefs_made_all()) make_coefficients();

  fr = 0.0;
  fz = 0.0;
  fp = 0.0;
  p = 0.0;
  d = 0.0;

  double rr = sqrt(r*r + z*z);

  if (rr/ASCALE>Rtable) {
    p = -cylmass/(rr+1.0e-16);
    fr = p*r/(rr+1.0e-16)/(rr+1.0e-16);
    fz = p*z/(rr+1.0e-16)/(rr+1.0e-16);

    return;
  }

  if (z/ASCALE > Rtable) z =  Rtable*ASCALE;
  if (z/ASCALE <-Rtable) z = -Rtable*ASCALE;

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

  for (int M=0; M<=multistep; M++) {

    for (int mm=0; mm<=MMAX; mm++) {

      for (int l=mm; l<=LMAX; l++) {

	out << setw(4) << M << setw(4) << mm << setw(4) << l << endl;

	for (int j=1; j<=NMAX; j++)
	  out << " " << setw(15) << accum_cos0[M][0][mm][j + NMAX*(l-mm)]*znorm;
	out << endl;

	if (mm) {

	  out << setw(4) << M << setw(4) << mm << setw(4) << l << endl;

	  for (int j=1; j<=NMAX; j++)
	    out << " " << setw(15) << accum_sin0[M][0][mm][j + NMAX*(l-mm)]*znorm;
	  out << endl;
	  
	}
      }
    }
  }

  out << endl;

  for (int mm=0; mm<=MMAX; mm++) {

    out << setw(4) << mm << setw(5) << "C" << endl;

    for (int n=0; n<rank3; n++) 
      out << " " << setw(15) << accum_cos[mm][n];
    out << endl;

    if (mm) {

      out << setw(4) << mm << setw(5) << "S" << endl;

      for (int n=0; n<rank3; n++) 
	out << " " << setw(15) << accum_sin[mm][n];
      out << endl;
	  
    }
  }

  out << endl;

  for (int M=0; M<=multistep; M++) {

    for (int mm=0; mm<=MMAX; mm++) {

      out << setw(4) << M << setw(4) << mm << setw(5) << "C" << endl;

      for (int n=0; n<rank3; n++)
	out << " " << setw(15) << accum_cosN[M][0][mm][n];
      out << endl;

      if (mm) {
	  
	out << setw(4) << M << setw(4) << mm << setw(5) << "S" << endl;
	
	for (int n=0; n<rank3; n++)
	  out << " " << setw(15) << accum_sinN[M][0][mm][n];
	out << endl;
	  
      }
    }
  }

}


#ifdef STANDALONE
#include <coef.H>
static CoefHeader coefheader;
static CoefHeader2 coefheader2;
#endif

void EmpCylSL::dump_coefs_binary_last(ostream& out, double time)
{
  double p, znorm = 1.0;

  coefheader.time = time;
  coefheader.mmax = MMAX;
  coefheader.nord = NORDER;
  coefheader.nmax = NMAX;

  out.write((const char *)&coefheader, sizeof(CoefHeader));

  for (int mm=0; mm<=MMAX; mm++) {

    for (int l=mm; l<=LMAX; l++) {
      
      for (int j=1; j<=NMAX; j++)
	out.write((const char *)&(p=accum_cos0[M][0][mm][j+NMAX*(l-mm)]*znorm), sizeof(double));
      
      if (mm) {

	for (int j=1; j<=NMAX; j++)
	  out.write((const char *)&(p=accum_sin0[M][0][mm][j+NMAX*(l-mm)]*znorm), sizeof(double));

      }
      
    }
  }
}

void EmpCylSL::dump_coefs_binary_curr(ostream& out, double time)
{
  coefheader2.time = time;
  coefheader2.mmax = MMAX;
  coefheader2.nmax = rank3;

  out.write((const char *)&coefheader2, sizeof(CoefHeader2));

  for (int mm=0; mm<=MMAX; mm++) {

    for (int j=0; j<rank3; j++)
      out.write((const char *)&accum_cos[mm][j], sizeof(double));
    
    if (mm) {

      for (int j=0; j<rank3; j++)
	out.write((const char *)&accum_sin[mm][j], sizeof(double));
      
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

  double fac=1;
  int n, mm;
  ofstream** outC = new ofstream* [MPItable];
  ofstream** outS = new ofstream* [MPItable];
  
  for (mm=0; mm<=MMAX; mm++) {

    for (n=0; n<=min<int>(NOUT, rank3-1); n++) {

				// Make output streams
      for (int i=0; i<MPItable; i++) {

	ostringstream ins;
	ins << name << ".C." << labels[i] << mm << "." << n
	    << "." << step;
	
	outC[i] = new ofstream(ins.str().c_str());
	outC[i]->write((const char *)&numx, sizeof(int));
	outC[i]->write((const char *)&numy, sizeof(int));
	outC[i]->write((const char *)&(zz=  0.0), sizeof(float));
	outC[i]->write((const char *)&(zz= rmax), sizeof(float));
	outC[i]->write((const char *)&(zz=-rmax), sizeof(float));
	outC[i]->write((const char *)&(zz= rmax), sizeof(float));
      }

      if (mm) {

	for (int i=0; i<MPItable; i++) {

	  ostringstream ins;
	  ins << name << ".S." << labels[i] << mm << "." << n 
	      << "." << step;
	
	  outS[i] = new ofstream(ins.str().c_str());
	  outS[i]->write((const char *)&numx,       sizeof(int));
	  outS[i]->write((const char *)&numy,       sizeof(int));
	  outS[i]->write((const char *)&(zz=  0.0), sizeof(float));
	  outS[i]->write((const char *)&(zz= rmax), sizeof(float));
	  outS[i]->write((const char *)&(zz=-rmax), sizeof(float));
	  outS[i]->write((const char *)&(zz= rmax), sizeof(float));
	}

      }


				// Ok, write data

      for (int k=0; k<numy; k++) {

	z = -rmax + dz*k;

	for (int j=0; j<numx; j++) {
	  
	  r = dr*j;

	  double X = (r_to_xi(r) - XMIN)/dX;
	  double Y = ( z_to_y(z) - YMIN)/dY;

	  int ix = (int)X;
	  int iy = (int)Y;

	  if (ix < 0) {
	    ix = 0;
	    if (enforce_limits) X = 0.0;
	  }
	  if (iy < 0) {
	    iy = 0;
	    if (enforce_limits) Y = 0.0;
	  }
	  
	  if (ix >= NUMX) {
	    ix = NUMX-1;
	    if (enforce_limits) X = NUMX;
	  }
	  if (iy >= NUMY) {
	    iy = NUMY-1;
	    if (enforce_limits) Y = NUMY;
	  }
	  
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

	  outC[0]->write((const char *)&zz, sizeof(float));
	    
	  zz = fac*(
	    rforceC[mm][n][ix  ][iy  ] * c00 +
	    rforceC[mm][n][ix+1][iy  ] * c10 +
	    rforceC[mm][n][ix  ][iy+1] * c01 +
	    rforceC[mm][n][ix+1][iy+1] * c11 );

	  outC[1]->write((const char *)&zz, sizeof(float));

	  zz = fac*(
	    zforceC[mm][n][ix  ][iy  ] * c00 +
	    zforceC[mm][n][ix+1][iy  ] * c10 +
	    zforceC[mm][n][ix  ][iy+1] * c01 +
	    zforceC[mm][n][ix+1][iy+1] * c11 );

	  outC[2]->write((const char *)&zz, sizeof(float));

	  if (DENS) {
	    zz = fac*(
	      densC[mm][n][ix  ][iy  ] * c00 +
	      densC[mm][n][ix+1][iy  ] * c10 +
	      densC[mm][n][ix  ][iy+1] * c01 +
	      densC[mm][n][ix+1][iy+1] * c11 );

	    outC[3]->write((const char *)&zz, sizeof(float));
	  }

	  if (mm) {
      
	    zz = fac*(
	      potS[mm][n][ix  ][iy  ] * c00 +
	      potS[mm][n][ix+1][iy  ] * c10 +
	      potS[mm][n][ix  ][iy+1] * c01 +
	      potS[mm][n][ix+1][iy+1] * c11 );

	    outS[0]->write((const char *)&zz, sizeof(float));
	    
	    zz = fac*(
	      rforceS[mm][n][ix  ][iy  ] * c00 +
	      rforceS[mm][n][ix+1][iy  ] * c10 +
	      rforceS[mm][n][ix  ][iy+1] * c01 +
	      rforceS[mm][n][ix+1][iy+1] * c11 );

	    outS[1]->write((const char *)&zz, sizeof(float));

	    zz = fac*(
	      zforceS[mm][n][ix  ][iy  ] * c00 +
	      zforceS[mm][n][ix+1][iy  ] * c10 +
	      zforceS[mm][n][ix  ][iy+1] * c01 +
	      zforceS[mm][n][ix+1][iy+1] * c11 );
	    
	    outS[2]->write((const char *)&zz, sizeof(float));
	    
	    if (DENS) {
	      zz = fac*(
		densS[mm][n][ix  ][iy  ] * c00 +
		densS[mm][n][ix+1][iy  ] * c10 +
		densS[mm][n][ix  ][iy+1] * c01 +
		densS[mm][n][ix+1][iy+1] * c11 );
	      
	      outS[3]->write((const char *)&zz, sizeof(float));
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

void EmpCylSL::dump_images(const string& OUTFILE,
			   double XYOUT, double ZOUT, int OUTR, int OUTZ,
			   bool logscale)
{
  if (myid!=0) return;
  
  double p, d, rf, zf, pf;
  double dr, dz = 2.0*ZOUT/(OUTZ-1);
  double rmin = RMIN*ASCALE;
  
  if (logscale) 
    dr = (log(XYOUT) - log(rmin))/(OUTR-1);
  else
    dr = (XYOUT - rmin)/(OUTR-1);
  
  string Name;
  int Number  = 14;
  int Number2 = 8;
  string Types[] = {".pot", ".dens", ".fr", ".fz", ".fp",
				// Finite difference
		    ".empfr", ".empfz", ".empfp",
				// Compare with recursion
		    ".diffr", ".diffz", ".diffp",
				// Relative difference
		    ".relfr", ".relfz", ".relfp"};
  
  //============
  // Open files
  //============
  ofstream *out = new ofstream [Number];
  for (int j=0; j<Number; j++) {
    Name = OUTFILE + Types[j] + ".eof_recon";
    out[j].open(Name.c_str());
    if (!out[j]) {
      cerr << "Couldn't open <" << Name << ">\n";
      break;
    }
  }
  
  double del, tmp, rP, rN, zP, zN, pP, pN, p0, d0;
  double potpr, potnr, potpz, potnz, potpp, potnp;

  
  double hr = HFAC*dr;
  double hz = HFAC*dz;
  double hp = 2.0*M_PI/64.0;
  
  if (logscale) hr = HFAC*0.5*(exp(dr) - exp(-dr));
  
  double r=0.0, z=0.0;
  
  for (int iz=0; iz<OUTZ; iz++) {
    
    z = -ZOUT + dz*iz;

    for (int ir=0; ir<OUTR; ir++) {
      
      if (logscale)
	r = rmin*exp(dr*ir);
      else
	r = rmin + dr*ir;
	  
      accumulated_eval(r, z, 0.0, p0, p, rf, zf, pf);
      d = accumulated_dens_eval(r, z, 0.0, d0);
      
      out[0]  << setw(15) << r << setw(15) << z << setw(15) << p;
      out[1]  << setw(15) << r << setw(15) << z << setw(15) << d;
      out[2]  << setw(15) << r << setw(15) << z << setw(15) << rf;
      out[3]  << setw(15) << r << setw(15) << z << setw(15) << zf;
      out[4]  << setw(15) << r << setw(15) << z << setw(15) << pf;

      
      //===================
      // Finite difference
      //===================
      
      if (logscale) {
	rP = r*(1.0 + 0.5*hr);
	rN = max<double>(1.0e-5, r*(1.0-0.5*hr));
	del = rP - rN;
      }
      else {
	rP = dr*ir+0.5*hr;
	rN = max<double>(1.0e-5, dr*ir-0.5*hr);
	del = hr;
      }
      zP = -ZOUT + dz*iz + 0.5*hz;
      zN = -ZOUT + dz*iz - 0.5*hz;
      
      pP =  0.5*hp;
      pN = -0.5*hp;

      accumulated_eval(rP, z, 0.0, tmp, potpr, tmp, tmp, tmp);
      accumulated_eval(rN, z, 0.0, tmp, potnr, tmp, tmp, tmp);
      accumulated_eval(r, zP, 0.0, tmp, potpz, tmp, tmp, tmp);
      accumulated_eval(r, zN, 0.0, tmp, potnz, tmp, tmp, tmp);
      accumulated_eval(r,  z,  pP, tmp, potpp, tmp, tmp, tmp);
      accumulated_eval(r,  z,  pN, tmp, potnp, tmp, tmp, tmp);
      
      //==================================================

      out[5] << setw(15) << r << setw(15) << z 
	     << setw(15) << -(potpr - potnr)/del;
      
      out[6] << setw(15) << r << setw(15) << z << setw(15)
	     << -(potpz - potnz)/hz;
      
      out[7] << setw(15) << r << setw(15) << z << setw(15)
	     << -(potpp - potnp)/hp;
      
      //==================================================
      

      out[8] << setw(15) << r << setw(15) << z 
	     << setw(15) << (potnr - potpr)/del - rf;
      
      out[9] << setw(15) << r << setw(15) << z << setw(15)
	     << (potnz - potpz)/hz - zf;
      
      out[10] << setw(15) << r << setw(15) << z << setw(15)
	      << (potpp - potnp)/hp - pf;
      
      //==================================================
	  
	  
      out[11] << setw(15) << r << setw(15) << z 
	      << setw(15) << ((potnr - potpr)/del - rf)/(fabs(rf)+1.0e-18);
      
      out[12] << setw(15) << r << setw(15) << z << setw(15)
	      << ((potnz - potpz)/hz - zf)/(fabs(zf)+1.0e-18);
      
      out[13] << setw(15) << r << setw(15) << z << setw(15)
	      << ((potpp - potnp)/hp - pf)/(fabs(pf)+1.0e-18);
      
      //==================================================

      for (int j=0; j<Number; j++) out[j] << endl;

    }
    
    for (int j=0; j<Number; j++) out[j] << endl;
    
  }
  
  //
  // Close current files
  //
  for (int j=0; j<Number; j++) out[j].close();

  //
  // Open files (face on)
  //
  for (int j=0; j<Number2; j++) {
    Name = OUTFILE + Types[j] + ".eof_recon_face";
    out[j].open(Name.c_str());
    if (!out[j]) {
      cerr << "Couldn't open <" << Name << ">\n";
      break;
    }
  }
  
  double rr, pp, xx, yy, dxy = 2.0*XYOUT/(OUTR-1);
  
  for (int iy=0; iy<OUTR; iy++) {
    yy = -XYOUT + dxy*iy;
    
    for (int ix=0; ix<OUTR; ix++) {
      xx = -XYOUT + dxy*ix;
      
      rr = sqrt(xx*xx + yy*yy);
      pp = atan2(yy, xx);
      
      accumulated_eval(rr, 0.0, pp, p0, p, rf, zf, pf);
      d = accumulated_dens_eval(rr, 0.0, pp, d0);
      
      out[0]  << setw(15) << xx << setw(15) << yy << setw(15) << p;
      out[1]  << setw(15) << xx << setw(15) << yy << setw(15) << d;
      out[2]  << setw(15) << xx << setw(15) << yy << setw(15) << rf;
      out[3]  << setw(15) << xx << setw(15) << yy << setw(15) << zf;
      out[4]  << setw(15) << xx << setw(15) << yy << setw(15) << pf;

      //===================
      // Finite difference
      //===================
      
      if (logscale) {
	rP = rr*(1.0 + 0.5*hr);
	rN = max<double>(1.0e-5, rr*(1.0-0.5*hr));
	del = rP - rN;
      }
      else {
	rP = rr+0.5*hr;
	rN = max<double>(1.0e-5, rr-0.5*hr);
	del = hr;
      }
      zP =  0.5*hz;
      zN = -0.5*hz;
      
      pP = pp + 0.5*hp;
      pN = pp - 0.5*hp;

      accumulated_eval(rP, 0.0, pp, tmp, potpr, tmp, tmp, tmp);
      accumulated_eval(rN, 0.0, pp, tmp, potnr, tmp, tmp, tmp);
      accumulated_eval(rr,  zP, pp, tmp, potpz, tmp, tmp, tmp);
      accumulated_eval(rr,  zN, pp, tmp, potnz, tmp, tmp, tmp);
      accumulated_eval(rr, 0.0, pP, tmp, potpp, tmp, tmp, tmp);
      accumulated_eval(rr, 0.0, pN, tmp, potnp, tmp, tmp, tmp);
      
      //==================================================

      out[5] << setw(15) << xx << setw(15) << yy 
	     << setw(15) << -(potpr - potnr)/del;
      
      out[6] << setw(15) << xx << setw(15) << yy << setw(15)
	     << -(potpz - potnz)/hz;
      
      out[7] << setw(15) << xx << setw(15) << yy << setw(15)
	     << -(potpp - potnp)/hp;
      
      //==================================================

      for (int j=0; j<Number2; j++) out[j] << endl;

    }
    
    for (int j=0; j<Number2; j++) out[j] << endl;
    
  }
  
  //
  // Close current files and delete
  //
  for (int j=0; j<Number2; j++) out[j].close();
  delete [] out;
}


void EmpCylSL::dump_images_basis(const string& OUTFILE,
				 double XYOUT, double ZOUT, 
				 int OUTR, int OUTZ, bool logscale,
				 int M1, int M2, int N1, int N2)
{
  if (myid!=0) return;
  
  double p, d, rf, zf, pf;
  double dr, dz = 2.0*ZOUT/(OUTZ-1);
  double rmin = RMIN*ASCALE;
  
  if (logscale) 
    dr = (log(XYOUT) - log(rmin))/(OUTR-1);
  else
    dr = (XYOUT - rmin)/(OUTR-1);
  
  int Number  = 10;
  string Types[] = {".pot", ".dens", ".fr", ".fz", ".empfr", ".empfz", 
		    ".diffr", ".diffz", ".relfr", ".relfz"};
  
  double tmp, rP, rN, zP, zN, r, z, del;
  double potpr, potnr, potpz, potnz;
  
  ofstream *out = new ofstream [Number];
  
  double hr = HFAC*dr;
  double hz = HFAC*dz;
  
  if (logscale) hr = HFAC*0.5*(exp(dr) - exp(-dr));
  
  for (int m=M1; m<=M2; m++) {
    
    for (int n=N1; n<=N2; n++) {
      
      
      //============
      // Open files
      //============
      
      for (int j=0; j<Number; j++) {
	ostringstream fname;
	fname << OUTFILE << Types[j] << '.' << m << '.' << n << ".eof_basis";
	out[j].open(fname.str().c_str());
	if (out[j].bad()) {
	  cout << "EmpCylSL::dump_images_basis: failed to open " 
	       << fname.str() << endl;
	  return;
	}
      }
      
      for (int iz=0; iz<OUTZ; iz++) {
	for (int ir=0; ir<OUTR; ir++) {
	  
	  z = -ZOUT + dz*iz;
	  if (logscale)
	    r = rmin*exp(dr*ir);
	  else
	    r = rmin + dr*ir;
	  
	  get_all(m, n, r, z, 0.0, p, d, rf, zf, pf);
	  
	  out[0] << setw(15) << r << setw(15) << z << setw(15) << p;
	  out[1] << setw(15) << r << setw(15) << z << setw(15) << d;
	  out[2] << setw(15) << r << setw(15) << z << setw(15) << rf;
	  out[3] << setw(15) << r << setw(15) << z << setw(15) << zf;
	  
	  //===================
	  // Finite difference
	  //===================
	  
	  if (logscale) {
	    rP = r*(1.0 + 0.5*hr);
	    rN = max<double>(1.0e-5, r*(1.0-0.5*hr));
	    del = rP - rN;
	  }
	  else {
	    rP = dr*ir+0.5*hr;
	    rN = max<double>(1.0e-5, dr*ir-0.5*hr);
	    del = hr;
	  }
	  zP = -ZOUT + dz*iz + 0.5*hz;
	  zN = -ZOUT + dz*iz - 0.5*hz;
	  
	  get_all(m, n, rP, z, 0.0, potpr, tmp, tmp, tmp, tmp);
	  get_all(m, n, rN, z, 0.0, potnr, tmp, tmp, tmp, tmp);
	  get_all(m, n, r, zP, 0.0, potpz, tmp, tmp, tmp, tmp);
	  get_all(m, n, r, zN, 0.0, potnz, tmp, tmp, tmp, tmp);
	  
	  out[4] << setw(15) << r << setw(15) << z 
		 << setw(15) << -(potpr - potnr)/del ;
	  
	  out[5] << setw(15) << r << setw(15) << z << setw(15)
		 << -(potpz - potnz)/hz  << setw(15) << hz;
	  
	  out[6] << setw(15) << r << setw(15) << z << setw(15)
		 << -(potpr - potnr)/del - rf;
	  
	  out[7] << setw(15) << r << setw(15) << z << setw(15)
		 << -(potpz - potnz)/hz  - zf;
	  
	  
	  out[8] << setw(15) << r << setw(15) << z << setw(15)
		 << (-(potpr - potnr)/del - rf)/(fabs(rf)+1.0e-18);
	  
	  out[9] << setw(15) << r << setw(15) << z << setw(15)
		 << (-(potpz - potnz)/hz  - zf)/(fabs(zf)+1.0e-18);
	  
	  for (int j=0; j<Number; j++) out[j] << endl;
	  
	}
	
	for (int j=0; j<Number; j++) out[j] << endl;
	
      }
      
      for (int j=0; j<Number; j++) out[j].close();
    }
    
  }
  
  delete [] out;
}


double EmpCylSL::r_to_xi(double r)
{
  if (CMAP) {
    if (r<0.0) {
      ostringstream msg;
      msg << "radius=" << r << " < 0! [mapped]";
      bomb(msg.str());
    }
    return (r/ASCALE - 1.0)/(r/ASCALE + 1.0);
  } else {
    if (r<0.0)  {
      ostringstream msg;
      msg << "radius=" << r << " < 0!";
      bomb(msg.str());
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

#define MINEPS 1.0e-10

void EmpCylSL::legendre_R(int lmax, double x, Matrix& p)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;

  p[0][0] = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[m][m] = pll;
      if (isnan(p[m][m]))
	cerr << "legendre_R: p[" << m << "][" << m << "]: pll=" << pll << "\n";
      fact += 2.0;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p[m][m];
    p[m+1][m] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[l][m] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      if (isnan(p[l][m]))
	cerr << "legendre_R: p[" << l << "][" << m << "]: pll=" << pll << "\n";

      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (isnan(x))
    cerr << "legendre_R: x\n";
  for(l=0; l<=lmax; l++)
    for (m=0; m<=l; m++)
      if (isnan(p[l][m]))
	cerr << "legendre_R: p[" << l << "][" << m << "] lmax=" 
	     << lmax << "\n";

}

void EmpCylSL::dlegendre_R(int lmax, double x, Matrix &p, Matrix &dp)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;

  p[0][0] = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[m][m] = pll;
      fact += 2.0;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p[m][m];
    p[m+1][m] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[l][m] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (1.0-fabs(x) < MINEPS) {
    if (x>0) x =   1.0 - MINEPS;
    else     x = -(1.0 - MINEPS);
  }

  somx2 = 1.0/(x*x - 1.0);
  dp[0][0] = 0.0;
  for (l=1; l<=lmax; l++) {
    for (m=0; m<l; m++)
      dp[l][m] = somx2*(x*l*p[l][m] - (l+m)*p[l-1][m]);
    dp[l][l] = somx2*x*l*p[l][l];
  }
}

void EmpCylSL::sinecosine_R(int mmax, double phi, Vector& c, Vector& s)
{
  int m;

  c[0] = 1.0;
  s[0] = 0.0;

  c[1] = cos(phi);
  s[1] = sin(phi);

  for (m=2; m<=mmax; m++) {
    c[m] = 2.0*c[1]*c[m-1] - c[m-2];
    s[m] = 2.0*c[1]*s[m-1] - s[m-2];
  }
}


void EmpCylSL::multistep_update_begin()
{
				// Clear the update matricies
  for (int M=0; M<=multistep; M++) {
    differC[M].setsize(0, MMAX, 0, rank3-1);
    differS[M].setsize(0, MMAX, 0, rank3-1);

    for (int mm=0; mm<=MMAX; mm++) {
      for (int nn=0; nn<rank3; nn++) {
	differC[M][mm][nn] = differS[M][mm][nn] = 0.0;
      }
    }
  }
}

void EmpCylSL::multistep_update_finish()
{
  vector<double> work(rank3);
				// Combine the update matricies
  for (int M=0; M<=multistep; M++) {

    for (int mm=0; mm<=MMAX; mm++) {
      
      MPI_Allreduce (&differC[M][mm][0], &work[0],
		     rank3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for (int nn=0; nn<rank3; nn++) accum_cos[mm][nn] += work[nn];

      if (mm) {
	MPI_Allreduce (&differS[M][mm][1], &work[0],
		       rank3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for (int nn=0; nn<rank3; nn++) accum_sin[mm][nn] += work[nn];
      }
    }

  }

}

void EmpCylSL::multistep_update(int from, int to, double r, double z, double phi, double mass)
{
  double rr = sqrt(r*r+z*z);
  if (rr/ASCALE>Rtable) return;

  double msin, mcos;

  double norm = -4.0*M_PI;
  
  get_pot(vc[0], vs[0], r, z);

  for (int mm=0; mm<=MMAX; mm++) {

    mcos = cos(phi*mm);
    msin = sin(phi*mm);

    for (int nn=0; nn<rank3; nn++) 
      hold[0][nn] = norm * mass * mcos * vc[0][mm][nn];

    for (int nn=0; nn<rank3; nn++) {
      differC[from][mm][nn] -= hold[0][nn];
      differC[to  ][mm][nn] += hold[0][nn];
    }

    if (mm>0) {
      for (int nn=0; nn<rank3; nn++) 
	hold[0][nn] = norm * mass * msin * vs[0][mm][nn];
      for (int nn=0; nn<rank3; nn++) {
	differS[from][mm][nn] -= hold[0][nn];
	differS[to  ][mm][nn] += hold[0][nn];
      }
    }
  }

}



void EmpCylSL::compute_multistep_coefficients(unsigned mlevel)
{
				// Clean coefficient matrix
				// 
  for (int mm=0; mm<=MMAX; mm++) {
    for (int nn=0; nn<rank3; nn++) {
      accum_cos[mm][nn] = 0.0;
      if (mm) accum_sin[mm][nn] = 0.0;
    }
  }
  
				// Interpolate to get coefficients above
  double a, b;			// 
  for (int M=0; M<mlevel; M++) {
    b = (double)(mstep - stepL[M])/(double)(stepN[M] - stepL[M]);
    a = 1.0 - b;
    for (int mm=0; mm<=MMAX; mm++) {
      for (int nn=0; nn<rank3; nn++) {
	accum_cos[mm][nn] += a*accum_cosL[M][0][mm][nn] + b*accum_cosN[M][0][mm][nn];
	if (mm)
	  accum_sin[mm][nn] += a*accum_sinL[M][0][mm][nn] + b*accum_sinN[M][0][mm][nn];
      }
    }
    // Sanity debug check
    if (a<0.0 && a>1.0) {
      cout << "Process " << myid << ": interpolation error in multistep [a]" << endl;
    }
    if (b<0.0 && b>1.0) {
      cout << "Process " << myid << ": interpolation error in multistep [b]" << endl;
    }
  }
				// Add coefficients at or below this level
				// 
  for (int M=mlevel; M<=multistep; M++) {
    for (int mm=0; mm<=MMAX; mm++) {
      for (int nn=0; nn<rank3; nn++) {
	accum_cos[mm][nn] += accum_cosN[M][0][mm][nn];
	if (mm)
	  accum_sin[mm][nn] += accum_sinN[M][0][mm][nn];
      }
    }
  }

}

//
// Swap pointers rather than copy
//
void EmpCylSL::multistep_swap(unsigned M)
{
  Vector **p;

  p = accum_cosL[M];
  accum_cosL[M] = accum_cosN[M];
  accum_cosN[M] = p;

  p = accum_sinL[M];
  accum_sinL[M] = accum_sinN[M];
  accum_sinN[M] = p;

  setup_accumulation(M);
}

void EmpCylSL::multistep_debug()
{
  for (int n=0; n<numprocs; n++) {
    if (n==myid) {
      cout << "Process " << myid << endl
	   << "   accum_cos[0][0]="
	   << accum_cos[0][0] << endl
	   << "   accum_sin[1][0]="
	   << accum_sin[1][0] << endl;

      cout.setf(ios::scientific);
      int c = cout.precision(2);
      for (int M=0; M<=multistep; M++) {
	cout << "   M=" << M << ": #part=" << howmany[M] << endl;

	cout << "   M=" << M << ", m=0, C_N: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << accum_cosN[M][0][0][j];
	cout << endl;

	cout << "   M=" << M << ", m=0, C_L: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << accum_cosL[M][0][0][j];
	cout << endl;

	cout << "   M=" << M << ", m=1, C_N: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << accum_cosN[M][0][1][j];
	cout << endl;

	cout << "   M=" << M << ", m=1, C_L: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << accum_cosL[M][0][1][j];
	cout << endl;

	cout << "   M=" << M << ", m=1, S_N: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << accum_sinN[M][0][1][j];
	cout << endl;

	cout << "   M=" << M << ", m=1, S_L: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << accum_sinL[M][0][1][j];
	cout << endl;

	cout.precision(c);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

