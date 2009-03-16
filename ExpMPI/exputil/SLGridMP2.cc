using namespace std;

// #define DEBUG 1
// #define DEBUG_SLEDGE 1
// #define DEBUG_NAN 1
#define NOTIFY 1
#define SLEDGE_VERBOSE 1
#define USE_TABLE 1

#define CYLFAC 1.01
// #define EXPONCYL
#define PLUMMERCYL
// #define JAFFECYL
// #define MESTELCYL

#define XOFFSET (1.0e-8)

#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>

#include <SLGridMP2.h>
#include <massmodel.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

// For fortran call
// (This should work both 32-bit and 64-bit . . . )
//
typedef int	logical;
typedef double	doublereal;
typedef int	integer;

MPI_Status status;

int SLGridCyl::mpi = 0;		// initially off
int SLGridCyl::cache = 1;	// initially yes
double SLGridCyl::A = 1.0;

extern "C" {
  int sledge_(logical* job, doublereal* cons, logical* endfin, 
	      integer* invec, doublereal* tol, logical* type, 
	      doublereal* ev, integer* numx, doublereal* xef, doublereal* ef, 
	      doublereal* pdef, doublereal* t, doublereal* rho, 
	      integer* iflag, doublereal* store);
  double iv(double, double);
  double kn(int, double);
  double i0(double);
  double i1(double);
  double k0(double);
  double k1(double);
}


// Unit density exponential disk with scale length A


double cylpot(double r)
{
				// Jaffe-like cylinder
#ifdef JAFFECYL
  double r2 = r/SLGridCyl::A;
  return log(1.0 + r2) - log(1.0 + CYLFAC*cylrfac);
#endif
				// Mestel-like cylinder
#ifdef MESTELCYL
  double a2 = SLGridCyl::A * SLGridCyl::A;
  double r2 = r*r/a2;
  return log(1.0 + sqrt(1.0 + r2)) - log(1.0 + sqrt(1.0 + CYLFAC*cylrfac));
#endif

				// Thin expontial
#ifdef EXPONCYL
  double y = 0.5 * r / SLGridCyl::A;
  return -2.0*M_PI*SLGridCyl::A*y*(i0(y)*k1(y) - i1(y)*k0(y));
#endif
				// Plummer-like cylinder
#ifdef PLUMMERCYL
  double a2 = SLGridCyl::A*SLGridCyl::A;
  double r2 = r*r/a2;

  return -1.0/sqrt(1.0 + r2);
#endif
}

double cyldpot(double r)
{
				// Jaffe-like cylinder
#ifdef JAFFECYL
  return 1.0/(SLGridCyl::A + r);
#endif
				// Mestel-like cylinder
#ifdef MESTELCYL
  double a2 = SLGridCyl::A * SLGridCyl::A;
  double fac = sqrt(1.0 + r*r/a2);
  return r/a2/(fac*(1.0 + fac));
#endif

				// Thin expontial
#ifdef EXPONCYL
  double y = 0.5 * r / SLGridCyl::A;
  return 4.0*M_PI*SLGridCyl::A*y*y*(i0(y)*k0(y) - i1(y)*k1(y));
#endif
				// Plummer-like cylinder
#ifdef PLUMMERCYL
  double a2 = SLGridCyl::A*SLGridCyl::A;
  double fac = sqrt(1.0 + r*r/a2);

  return r/a2/(fac*fac*fac);
#endif
}

double cyldens(double r)
{
				// Jaffe-like cylinder
#ifdef JAFFECYL
  double fac = r + SLGridCyl::A;
  return SLGridCyl::A/( fac*fac*r );
#endif

				// Mestel-like cylinder
#ifdef MESTELCYL
  double a2 = SLGridCyl::A * SLGridCyl::A;
  double r2 = r*r/a2;
  double fac1 = sqrt(1.0 + r2);
  double fac2 = 1.0 + fac1;

  return 1.0/(fac1*fac2*a2) * (2.0 - r2/(fac1*fac2) - r2/(fac1*fac1));
#endif

  // This 4pi from Poisson's eqn
  //        |
  //        |       /-- This begins the true projected density profile
  //        |       |
  //        v       v
#ifdef EXPONCYL
  return 4.0*M_PI * exp(-r/SLGridCyl::A);
#endif

#ifdef PLUMMERCYL
  /*
  double a2 = SLGridCyl::A*SLGridCyl::A;
  double r2 = r*r/a2;

  return (2.0 - r2)/pow(1.0 + r2, 2.5)/a2;
  */
  return 4.0*M_PI * exp(-r/SLGridCyl::A);
#endif
}


double cylpotsl(double r)
{
				// Jaffe cylinder
#ifdef JAFFECYL
  return cyldens(r);
#endif
				// Mestel cylinder
#ifdef MESTELCYL
  return cyldens(r);
#endif
				// Exponential
#ifdef EXPONCYL
  double y = 0.5 * r / SLGridCyl::A;
  return M_PI*(2.0*SLGridCyl::A*i0(y)*k0(y) - r*i0(y)*k1(y) +
	       r*i1(y)*k0(y))/(SLGridCyl::A*SLGridCyl::A);
#endif
  
				// Plummer-like
#ifdef PLUMMERCYL
  double a2 = SLGridCyl::A*SLGridCyl::A;
  double r2 = r*r/a2;

  return (2.0 - r2)/pow(1.0 + r2, 2.5)/a2;
#endif
}

void SLGridCyl::bomb(string oops)
{
  cerr << "SLGridCyl: " << oops << endl; 
  exit(-1);
}

				// Constructors

SLGridCyl::SLGridCyl(int MMAX, int NMAX, int NUMR, int NUMK, 
	       double RMIN, double RMAX, double L, int CMAP, double SCALE)
{
  int m, k;

  mmax = MMAX;
  nmax = NMAX;
  numr = NUMR;
  numk = NUMK;

  rmin = RMIN;
  rmax = RMAX;
  l = L;

  cmap = CMAP;
  scale = SCALE;

#ifdef JAFFECYL
  cylrfac = RMAX/A;
#endif
#ifdef MESTELCYL
  cylrfac = RMAX*RMAX/A*A;
#endif

  kv.setsize(0, NUMK);
  // dk = M_PI/L;
  dk = 0.5*M_PI/L;
  for (k=0; k<=NUMK; k++) kv[k] = dk*k;

  init_table();


#ifdef DEBUG
  if (mpi)
    cout << "Process " << myid << ": MPI is on!\n";
  else
    cout << "Process " << myid << ": MPI is off!\n";
#endif

  if (mpi) {

    table =  new TableCyl* [mmax+1];
    for (m=0; m<=mmax; m++) table[m] = new TableCyl [numk+1];

    mpi_setup();

    if (mpi_myid) {

      compute_table_slave();

      //
      // <Receive completed table from master>
      //

      for (m=0; m<=mmax; m++) {
	for (k=0; k<=numk; k++) {
	  /*
	  MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, 0,
		   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		   */
	  MPI_Bcast(mpi_buf, mpi_bufsz, MPI_PACKED, 0, MPI_COMM_WORLD);
    
	  mpi_unpack_table();      
	}
      }
      
    }
    else {			// BEGIN Master

      int slave = 0;
      int request_id = 1;

      if (!read_cached_table()) {

	double K;
	m=0; k=0;

	while (m<=mmax) {

	  if (slave<mpi_numprocs-1 && m<=mmax) { // Send request to slave
	    slave++;
      
#ifdef DEBUG    
	    cout << "Master sending orders to Slave " << slave << ": (m,k)=("
		 << m << ", " << k << ")\n";
#endif
	    MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	    MPI_Send(&m, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	    MPI_Send(&k, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
      
#ifdef DEBUG    
	    cout << "Master gave orders to Slave " << slave << ": (m,k)=("
		 << m << ", " << ")\n";
#endif

				// Increment counters
	    k++;
	    if (k>numk) {
	      k=0;
	    m++;
	    }
	  }
	  
	  if (slave == mpi_numprocs-1 && m<=mmax) {
	    
	    //
	    // <Wait and receive>
	    //
	
	    MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 
		   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
	    int retid = status.MPI_SOURCE;

	    mpi_unpack_table();      

	    //
	    // <Send new request>
	    //

	    K = dk*k;

	    MPI_Send(&request_id, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	    MPI_Send(&m, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	    MPI_Send(&k, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
      
#ifdef DEBUG    
	    cout << "Master gave orders to Slave " << retid << ": (m,k)=("
		 << m << ", " << k << ")\n";
#endif
				// Increment counters
	    k++;
	    if (k>numk) {
	      k=0;
	      m++;
	    }
	  }
	}
      
	//
	// <Wait for all slaves to return>
	//
  
	while (slave) {
	
	
	  MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 
		   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
	  mpi_unpack_table();      

	  slave--;
	}

	if (cache) write_cached_table();

      }


      //
      // <Tell slaves to continue>
      //

      request_id = -1;
      for (slave=1; slave < mpi_numprocs; slave++)
	MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);


      //
      // <Send table to slaves>
      //

      for (m=0; m<=mmax; m++) {
	for (k=0; k<=numk; k++) {
	  int position = mpi_pack_table(&table[m][k], m, k);
	  /*
	  for (slave=1; slave < mpi_numprocs; slave++)
	    MPI_Send(mpi_buf, position, MPI_PACKED, slave, 11, MPI_COMM_WORLD);
	    */
	  MPI_Bcast(mpi_buf, position, MPI_PACKED, 0, MPI_COMM_WORLD);
	}
      }


    } // END Master

  }
  else {

    table =  new TableCyl* [mmax+1];
    for (m=0; m<=mmax; m++) table[m] = new TableCyl [numk+1];

    if (!cache || !read_cached_table()) {
      for (m=0; m<=mmax; m++) {
	for (k=0; k<=numk; k++) {
	  cerr << "Begin [" << m << ", " << k << "] . . .\n";
	  compute_table(&(table[m][k]), m, k);
	  cerr << ". . . done\n";
	}
      }
      if (cache) write_cached_table();
    }
  }
}

const string cyl_cache_name = ".slgrid_cyl_cache";

int SLGridCyl::read_cached_table(void)
{
  if (!cache) return 0;

  ifstream in(cyl_cache_name.c_str());
  if (!in) return 0;

  int MMAX, NMAX, NUMR, NUMK, i, j, CMAP;
  double RMIN, RMAX, L, AA, SCL;

  if (myid==0)
    cerr << "SLGridCyl::read_cached_table: trying to read cached table . . ."
	 << endl;

  in.read((char *)&MMAX, sizeof(int));		if(!in || MMAX!=mmax) return 0;
  in.read((char *)&NMAX, sizeof(int));		if(!in || NMAX!=nmax) return 0;
  in.read((char *)&NUMR, sizeof(int));		if(!in || NUMR!=numr) return 0;
  in.read((char *)&NUMK, sizeof(int));		if(!in || NUMK!=numk) return 0;
  in.read((char *)&RMIN, sizeof(double));	if(!in || RMIN!=rmin) return 0;
  in.read((char *)&RMAX, sizeof(double));	if(!in || RMAX!=rmax) return 0;
  in.read((char *)&L, sizeof(double));		if(!in || L!=l) return 0;
  in.read((char *)&AA, sizeof(double));		if(!in || AA!=A) return 0;
  in.read((char *)&CMAP, sizeof(double));	if(!in || CMAP!=cmap) return 0;
  in.read((char *)&SCL, sizeof(double));	if(!in || SCL!=scale) return 0;

  for (int m=0; m<=mmax; m++) {
    for (int k=0; k<=numk; k++) {

      in.read((char *)&table[m][k].m, sizeof(int));
      in.read((char *)&table[m][k].k, sizeof(int));

				// Double check
      if (table[m][k].m != m) {
	cerr << "SLGridCyl: error reading <" << cyl_cache_name << ">\n";
	cerr << "SLGridCyl: m: read value (" << table[m][k].m << ") != internal value (" << m << ")\n";
	return 0;
      }
      if (table[m][k].k != k) {
	cerr << "SLGridCyl: error reading <" << cyl_cache_name << ">\n";
	cerr << "SLGridCyl: k: read value (" << table[m][k].k << ") != internal value (" << k << ")\n";
	return 0;
      }

      table[m][k].ev.setsize(1, nmax);
      table[m][k].ef.setsize(1, nmax, 0, numr-1);

      for (j=1; j<=nmax; j++) in.read((char *)&table[m][k].ev[j], sizeof(double));

      for (j=1; j<=nmax; j++)
	for (i=0; i<numr; i++)
	  in.read((char *)&table[m][k].ef[j][i], sizeof(double));
    }
  }

  if (myid==0)
    cerr << "SLGridCyl::read_cached_table: Success!!" << endl;
  return 1;
}


void SLGridCyl::write_cached_table(void)
{
  ofstream out(cyl_cache_name.c_str());
  if (!out) {
    cerr << "SLGridCyl: error writing <" << cyl_cache_name << ">\n";
    return;
  }

  int i, j;

  out.write((char *)&mmax, sizeof(int));
  out.write((char *)&nmax, sizeof(int));
  out.write((char *)&numr, sizeof(int));
  out.write((char *)&numk, sizeof(int));
  out.write((char *)&rmin, sizeof(double));
  out.write((char *)&rmax, sizeof(double));
  out.write((char *)&l, sizeof(double));
  out.write((char *)&A, sizeof(double));
  out.write((char *)&cmap, sizeof(int));
  out.write((char *)&scale, sizeof(double));

  for (int m=0; m<=mmax; m++) {
    for (int k=0; k<=numk; k++) {

      out.write((char *)&table[m][k].m, sizeof(int));
      out.write((char *)&table[m][k].k, sizeof(int));

      for (j=1; j<=nmax; j++)
	out.write((char *)&table[m][k].ev[j], sizeof(double));

      for (j=1; j<=nmax; j++)
	for (i=0; i<numr; i++)
	  out.write((char *)&table[m][k].ef[j][i], sizeof(double));
    }
  }

  cerr << "SLGridCyl::write_cached_table: done!!\n";
  return ;
}


SLGridCyl::~SLGridCyl()
{
  for (int m=0; m<=mmax; m++) delete [] table[m];
  delete [] table;
}

				// Members

double SLGridCyl::r_to_xi(double r)
{
  if (cmap) {
    if (r<0.0) bomb("radius < 0!");
    return (r/scale-1.0)/(r/scale+1.0);
  } else {
    if (r<0.0) bomb("radius < 0!");
    return r;
  }
}
    
double SLGridCyl::xi_to_r(double xi)
{
  if (cmap) {
    if (xi<-1.0) bomb("xi < -1!");
    if (xi>=1.0) bomb("xi >= 1!");

    return (1.0+xi)/(1.0 - xi) * scale;
  } else {
    return xi;
  }

}

double SLGridCyl::d_xi_to_r(double xi)
{
  if (cmap) {
    if (xi<-1.0) bomb("xi < -1!");
    if (xi>=1.0) bomb("xi >= 1!");

    return 0.5*(1.0-xi)*(1.0-xi)/scale;
  } else {
    return 1.0;
  }
}

double SLGridCyl::get_pot(double x, int m, int n, int k, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

#ifdef USE_TABLE
  return (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
    sqrt(fabs(table[m][k].ev[n])) * (x1*p0[indx] + x2*p0[indx+1]);
#else
  return (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
    sqrt(fabs(table[m][k].ev[n])) * cylpot(xi_to_r(x));
#endif
}


double SLGridCyl::get_dens(double x, int m, int n, int k, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  
#ifdef USE_TABLE
  return (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1]) *
    sqrt(fabs(table[m][k].ev[n])) * (x1*d0[indx] + x2*d0[indx+1])
    * table[m][k].ev[n]/fabs(table[m][k].ev[n]);
#else
  return (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1]) *
    sqrt(fabs(table[m][k].ev[n])) * cyldens(xi_to_r(x))
    * table[m][k].ev[n]/fabs(table[m][k].ev[n]);
#endif

}

double SLGridCyl::get_force(double x, int m, int n, int k, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  
				// Use three point formula

				// Point -1: indx-1
				// Point  0: indx
				// Point  1: indx+1

  return d_xi_to_r(x)/dxi * (
			     (p - 0.5)*table[m][k].ef[n][indx-1]*p0[indx-1]
			     -2.0*p*table[m][k].ef[n][indx]*p0[indx]
			     + (p + 0.5)*table[m][k].ef[n][indx+1]*p0[indx+1]
			     ) / sqrt(fabs(table[m][k].ev[n]));
}


void SLGridCyl::get_pot(Matrix& mat, double x, int m, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

  mat.setsize(0, numk, 1, nmax);

				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int k=0; k<=numk; k++) {
    for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
      mat[k][n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
	sqrt(fabs(table[m][k].ev[n])) * (x1*p0[indx] + x2*p0[indx+1]);
#else
      mat[k][n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
	sqrt(fabs(table[m][k].ev[n])) * cylpot(xi_to_r(x));
#endif
#ifdef DEBUG_NAN
      if (isnan(mat[k][n]) || isinf(mat[k][n]) ) {
	cerr << "SLGridCyl::get_pot: invalid value\n";
      }
#endif
    }
  }

}


void SLGridCyl::get_dens(Matrix& mat, double x, int m, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

  mat.setsize(0, numk, 1, nmax);


				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int k=0; k<=numk; k++) {
    for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
      mat[k][n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])*
	sqrt(fabs(table[m][k].ev[n])) * (x1*d0[indx] + x2*d0[indx+1])
      * table[m][k].ev[n]/fabs(table[m][k].ev[n]);
#else
      mat[k][n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])*
	sqrt(fabs(table[m][k].ev[n])) * cyldens(xi_to_r(x))
	* table[m][k].ev[n]/fabs(table[m][k].ev[n]);
#endif
    }
  }

}


void SLGridCyl::get_force(Matrix& mat, double x, int m, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

  mat.setsize(0, numk, 1, nmax);


				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int k=0; k<=numk; k++) {
    for (int n=1; n<=nmax; n++) {
      mat[k][n] = fac * (
			 (p - 0.5)*table[m][k].ef[n][indx-1]*p0[indx-1]
			 -2.0*p*table[m][k].ef[n][indx]*p0[indx]
			 + (p + 0.5)*table[m][k].ef[n][indx+1]*p0[indx+1]
			 ) / sqrt(fabs(table[m][k].ev[n]));
#ifdef DEBUG_NAN
      if (isnan(mat[k][n]) || isinf(mat[k][n]) ) {
	cerr << "SLGridCyl::get_force: invalid value\n";
      }
#endif
    }
  }

}


void SLGridCyl::get_pot(Vector& vec, double x, int m, int k, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

  vec.setsize(1, nmax);


				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
      sqrt(fabs(table[m][k].ev[n])) * (x1*p0[indx] + x2*p0[indx+1]);
#else
    vec[n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
      sqrt(fabs(table[m][k].ev[n])) * cylpot(xi_to_r(x));
#endif
  }

}


void SLGridCyl::get_dens(Vector& vec, double x, int m, int k, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

  vec.setsize(1, nmax);


				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])*
      sqrt(fabs(table[m][k].ev[n])) * (x1*d0[indx] + x2*d0[indx+1])
      * table[m][k].ev[n]/fabs(table[m][k].ev[n]);
#else
    vec[n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])*
      sqrt(fabs(table[m][k].ev[n])) * cyldens(xi_to_r(x))
      * table[m][k].ev[n]/fabs(table[m][k].ev[n]);
#endif
  }

}


void SLGridCyl::get_force(Vector& vec, double x, int m, int k, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

  vec.setsize(1, nmax);


				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int n=1; n<=nmax; n++) {
    vec[n] = fac * (
		    (p - 0.5)*table[m][k].ef[n][indx-1]*p0[indx-1]
		    -2.0*p*table[m][k].ef[n][indx]*p0[indx]
		    + (p + 0.5)*table[m][k].ef[n][indx+1]*p0[indx+1]
		    ) / sqrt(fabs(table[m][k].ev[n]));
  }

}


void SLGridCyl::get_pot(Matrix* mat, double x, int mMin, int mMax, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

  if (mmax < mMax) mMax = mmax;

  for (int m=mMin; m<=mMax; m++)
    mat[m].setsize(0, numk, 1, nmax);


				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int m=mMin; m<=mMax; m++) {
    for (int k=0; k<=numk; k++) {
      for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
	mat[m][k][n] = (x1*table[m][k].ef[n][indx] + 
			x2*table[m][k].ef[n][indx+1])/
	  sqrt(fabs(table[m][k].ev[n])) * (x1*p0[indx] + x2*p0[indx+1]);
#else
	mat[m][k][n] = (x1*table[m][k].ef[n][indx] + 
			x2*table[m][k].ef[n][indx+1])/
	  sqrt(fabs(table[m][k].ev[n])) * cylpot(xi_to_r(x));
#endif
#ifdef DEBUG_NAN
	if (isnan(mat[m][k][n]) || isinf(mat[m][k][n]) ) {
	  cerr << "SLGridCyl::get_pot: invalid value\n";
	  cerr <<   "  x1=" << x1
	       << "\n  x2=" << x2
	       << "\n  t0=" << table[m][k].ef[n][indx]
	       << "\n  p0=" << p0[indx]
	       << "\n  tp=" << table[m][k].ef[n][indx+1]
	       << "\n  pp=" << p0[indx+1]
	       << "\n  ev=" << fabs(table[m][k].ev[n])
	       << "\n val=" << (x1*table[m][k].ef[n][indx] + 
				x2*table[m][k].ef[n][indx+1])/
	    sqrt(fabs(table[m][k].ev[n])) * (x1*p0[indx] + x2*p0[indx+1])
	       << "\n";
	}
#endif
      }
    }
  }

}


void SLGridCyl::get_dens(Matrix* mat, double x, int mMin, int mMax, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

  if (mmax < mMax) mMax = mmax;

  for (int m=mMin; m<=mMax; m++)
    mat[m].setsize(0, numk, 1, nmax);


				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int m=mMin; m<=mMax; m++) {
    for (int k=0; k<=numk; k++) {
      for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
	mat[m][k][n] = (x1*table[m][k].ef[n][indx] + 
			x2*table[m][k].ef[n][indx+1])*
	  sqrt(fabs(table[m][k].ev[n])) * (x1*d0[indx] + x2*d0[indx+1])
	  * table[m][k].ev[n]/fabs(table[m][k].ev[n]);
#else
	mat[m][k][n] = (x1*table[m][k].ef[n][indx] + 
			x2*table[m][k].ef[n][indx+1])*
	  sqrt(fabs(table[m][k].ev[n])) * cyldens(xi_to_r(x))
	  * table[m][k].ev[n]/fabs(table[m][k].ev[n]);
#endif
      }
    }
  }

}


void SLGridCyl::get_force(Matrix* mat, double x, int mMin, int mMax, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-XOFFSET;
  }

  if (mmax < mMax) mMax = mmax;

  for (int m=mMin; m<=mMax; m++)
    mat[m].setsize(0, numk, 1, nmax);


				// XI grid is same for all k

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int m=mMin; m<=mMax; m++) {
    for (int k=0; k<=numk; k++) {
      for (int n=1; n<=nmax; n++) {
	mat[m][k][n] = fac * (
			      (p - 0.5)*table[m][k].ef[n][indx-1]*p0[indx-1]
			      -2.0*p*table[m][k].ef[n][indx]*p0[indx]
			      + (p + 0.5)*table[m][k].ef[n][indx+1]*p0[indx+1]
			      ) / sqrt(fabs(table[m][k].ev[n]));
#ifdef DEBUG_NAN
	if (isnan(mat[m][k][n]) || isinf(mat[m][k][n]) ) {
	  cerr << "SLGridCyl::get_force: invalid value\n";
	  cerr <<   "   p=" << p
	       << "\n  tm=" << table[m][k].ef[n][indx-1]
	       << "\n  pm=" << p0[indx-1]
	       << "\n  t0=" << table[m][k].ef[n][indx]
	       << "\n  p0=" << p0[indx]
	       << "\n  tp=" << table[m][k].ef[n][indx+1]
	       << "\n  pp=" << p0[indx+1]
	       << "\n  ev=" << fabs(table[m][k].ev[n])
	       << "\n val=" << fac * (
				      (p - 0.5)*table[m][k].ef[n][indx-1]*p0[indx-1]
				      -2.0*p*table[m][k].ef[n][indx]*p0[indx]
				      + (p + 0.5)*table[m][k].ef[n][indx+1]*p0[indx+1]
				      ) / sqrt(fabs(table[m][k].ev[n]))
	       << "\n";
	}
#endif
      }
    }
  }
  
}


static double L2, M2, K2;
static int sl_dim;

void SLGridCyl::compute_table(struct TableCyl* table, int m, int k)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol[6] = {1.0e-4*scale,1.0e-5,  
		   1.0e-4*scale,1.0e-5,  
		   1.0e-4*scale,1.0e-5};
  int i, j, VERBOSE=0;
  integer NUM, N, M;
  logical type[8];
  logical endfin[2] = {1, 1};
  
#ifdef DEBUG_SLEDGE
  if (myid==1) VERBOSE = SLEDGE_VERBOSE;
#endif

  cons[6] = rmin;
  cons[7] = rmax;
  M2 = m*m;
  K2 = kv[k]*kv[k];
  NUM = numr;
  N = nmax;
  M = m;

  // integer iflag[nmax], invec[nmax+3];
  integer *iflag = new integer [nmax];
  integer *invec = new integer [nmax+3];

  // double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
  // pdef[NUM*N];
  double *t=0, *rho=0;
  double *ev = new double [N];
  double *store = new double [26*(NUM+16)];
  double *xef = new double [NUM+16];
  double *ef = new double [NUM*N];
  double *pdef = new double [NUM*N];
  double f;

				// Inner  BC
  f = cylpot(cons[6]);
  if (M==0) {
    cons[0] = cyldpot(cons[6])/f;
    cons[2] = 1.0/(cons[6]*f*f);
  }
  else
    cons[0] = 1.0;

				// Outer BC
  f = cylpot(cons[7]);
  // cons[4] = (1.0+M)/cons[7] + cyldpot(cons[7])/f;
				// TEST
  cons[4] = (1.0+M)/(cons[7]*cons[7]) + cyldpot(cons[7])/f/cons[7];
  cons[5] = 1.0/(cons[7]*f*f);

  
  //
  //     Initialize the vector INVEC(*):
  //       estimates for the eigenvalues/functions specified
  //

  invec[0] = VERBOSE;		// little printing (1), no printing (0)
  invec[1] = 3;			// spectrum is ignored
  invec[2] = N;			// estimates for N eigenvalues/functions

  for (i=0; i<N; i++) invec[3+i] = i;

  //
  //     Set the JOB(*) vector:
  //        estimate both eigenvalues and eigenvectors,
  //        don't estimate the spectral density function,
  //        classify,
  //        let SLEDGE choose the initial mesh
  //
  logical job[5] = {0,1,0,0,0};

  //
  //     Output mesh
  //
  for (i=0; i<NUM; i++) xef[i] = r[i];

  //     
  //     Open file for output.
  //
  sl_dim = 2;

  sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	   t, rho, iflag, store);
  //
  //     Print results:
  //
  
  cout.precision(6);
  cout.setf(ios::scientific);

  for (i=0; i<N; i++) {
    cout << setw(15) << invec[3+i] 
	 << setw(15) << ev[i]
	 << setw(5) << iflag[i]
	 << endl;
  
    if (VERBOSE) {

      if (iflag[i] > -10) {
	cout << setw(14) << "x"
	     << setw(25) << "u(x)"
	     << setw(25) << "(pu`)(x)"
	     << endl;
	k = NUM*i;
	for (j=0; j<NUM; j++) {
	  cout << setw(25) << xef[j]
	       << setw(25) << ef[j+k]
	       << setw(25) << pdef[j+k]
	       << endl;
	}
      }

    }

  }
  
				// Load table

  table->ev.setsize(1, N);
  for (i=0; i<N; i++) table->ev[i+1] = ev[i];

  table->ef.setsize(1, N, 0, numr-1);

  for (i=0; i<numr; i++) {
    for(j=0; j<N; j++) 
      table->ef[j+1][i] = ef[j*NUM+i];
  }

  table->m = m;
  table->k = k;

  delete [] iflag;
  delete [] invec;
  delete [] ev;
  delete [] store;
  delete [] xef;
  delete [] ef;
  delete [] pdef;

}


void SLGridCyl::init_table(void)
{

  int i;

  xi.setsize(0, numr-1);
  r.setsize(0, numr-1);
  p0.setsize(0, numr-1);
  d0.setsize(0, numr-1);

  if (cmap) {
    xmin = (rmin/scale - 1.0)/(rmin/scale + 1.0);
    xmax = (rmax/scale - 1.0)/(rmax/scale + 1.0);
    dxi = (xmax-xmin)/(numr-1);
  } else {
    xmin = rmin;
    xmax = rmax;
    dxi = (xmax-xmin)/(numr-1);
  }

  for (i=0; i<numr; i++) {
    xi[i] = xmin + dxi*i;
    r[i] = xi_to_r(xi[i]);
    p0[i] = cylpot(r[i]);
    d0[i] = cyldens(r[i]);
  }

}


void SLGridCyl::compute_table_slave(void)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol[6] = {1.0e-4*scale,1.0e-5,  
		   1.0e-4*scale,1.0e-5, 
		   1.0e-4*scale,1.0e-5};
  int i, j, VERBOSE=0;
  integer NUM;
  logical type[8];
  logical endfin[2] = {1, 1};
  
  struct TableCyl table;
  int M, N, K;

#ifdef DEBUG_SLEDGE
  if (myid==1) VERBOSE = SLEDGE_VERBOSE;
#endif

#ifdef DEBUG
    cout << "Slave " << mpi_myid << " begins . . .\n";
#endif

  //
  // <Wait for orders>
  //
  
  int request_id;

  while(1) {

    MPI_Recv(&request_id, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (request_id < 0) break;	// Good-bye

    MPI_Recv(&M, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&K, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);


#ifdef DEBUG
    cout << "Slave " << mpi_myid << ": ordered to compute (m, k)=("
	 << M << ", " << K << ")\n";
#endif

    cons[0] = cons[1] = cons[2] = cons[3] = cons[4] = cons[5] = 0.0;
    cons[6] = rmin;
    cons[7] = rmax;
    M2 = M*M;
    K2 = kv[K]*kv[K];
    NUM = numr;
    N = nmax;

    // integer iflag[nmax], invec[nmax+3];
    integer *iflag = new integer [nmax];
    integer *invec = new integer [nmax+3];

    // double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
    // pdef[NUM*N];
    double *t=0, *rho=0;
    double *ev = new double [N];
    double *store = new double [26*(NUM+16)];
    double *xef = new double [NUM+16];
    double *ef = new double [NUM*N];
    double *pdef = new double [NUM*N];
    double f;

    f = cylpot(cons[6]);
    cons[2] = -1.0/(cons[6]*f);
    cons[4] = M/cons[7];
    f = cylpot(cons[7]);
    cons[5] = 1.0/(cons[7]*f*f);

    //
    //     Initialize the vector INVEC(*):
    //       estimates for the eigenvalues/functions specified
    //

    invec[0] = VERBOSE;		// little printing (1), no printing (0)
    invec[1] = 3;		// spectrum is ignored
    invec[2] = N;		// estimates for N eigenvalues/functions

    for (i=0; i<N; i++) invec[3+i] = i;

    //
    //     Set the JOB(*) vector:
    //        estimate both eigenvalues and eigenvectors,
    //        don't estimate the spectral density function,
    //        classify,
    //        let SLEDGE choose the initial mesh
    //
    logical job[5] = {0,1,0,0,0};

    //
    //     Output mesh
    //
    for (i=0; i<NUM; i++) xef[i] = r[i];

    //     
    //     Open file for output.
    //
    sl_dim = 2;

    sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	    t, rho, iflag, store);
    //
    //     Print results:
    //
    
#ifdef DEBUG
    
    if (type[0]) cout << "Slave " << myid 
		      << ": Inner endpoint is regular\n";
    if (type[1]) cout << "Slave " << myid 
		      << ": Inner endpoint is limit circle\n";
    if (type[2]) cout << "Slave " << myid 
		      << ": Inner endpoint is nonoscillatory for all EV\n";
    if (type[3]) cout << "Slave " << myid 
		      << ": Inner endpoint is oscillatory for all EV\n";
    if (type[4]) cout << "Slave " << myid 
		      << ": Outer endpoint is regular\n";
    if (type[5]) cout << "Slave " << myid 
		      << ": Outer endpoint is limit circle\n";
    if (type[6]) cout << "Slave " << myid 
		      << ": Outer endpoint is nonoscillatory for all EV\n";
    if (type[7]) cout << "Slave " << myid 
		      << ": Outer endpoint is oscillatory for all EV\n";
		      
    cout << "Slave " << mpi_myid << ": computed (m, k)=("
		 << M << ", " << K << ")\n";

    cout.precision(6);
    cout.setf(ios::scientific);

    for (i=0; i<N; i++) {
      cout << setw(15) << invec[3+i] 
	   << setw(15) << ev[i]
	   << setw(5) << iflag[i]
	   << endl;
  
      if (VERBOSE) {

	if (iflag[i] > -10) {
	  cout << setw(14) << "x"
	       << setw(25) << "u(x)"
	       << setw(25) << "(pu`)(x)"
	       << endl;
	  int k = NUM*i;
	  for (j=0; j<NUM; j++) {
	    cout << setw(25) << xef[j]
		 << setw(25) << ef[j+k]
		 << setw(25) << pdef[j+k]
		 << endl;
	  }
	}
	
      }

    }
  
#endif

#ifdef NOTIFY
    bool ok = true;
    for (i=0; i<N; i++) {
      if (iflag[i] < 0) {
	ok = false;
	cerr << "***** Slave " << setw(3) << myid 
	     << "  Level " << setw(3) << i << ": " << iflag[i] 
	     << "  M=" << M << "  kv[" << K << "]=" << kv[K] << endl;
      }
    }
    if (!ok) cerr << "***** Slave " << setw(3) << myid 
		  << ": if error=-2, consider increasing zmax" << endl;
#endif

				// Load table

    table.ev.setsize(1, N);
    for (i=0; i<N; i++) table.ev[i+1] = ev[i];

    table.ef.setsize(1, N, 0, numr-1);
    for (i=0; i<numr; i++) {
      for(j=0; j<N; j++) 
	table.ef[j+1][i] = ef[j*NUM+i];
    }

    table.m = M;
    table.k = K;
  
    int position = mpi_pack_table(&table, M, K);
    MPI_Send(mpi_buf, position, MPI_PACKED, 0, 11, MPI_COMM_WORLD);

#ifdef DEBUG
    cout << "Slave " << mpi_myid << ": sent to master (m, k)=("
		 << M << ", " << K << ")\n";
#endif

    delete [] iflag;
    delete [] invec;
    delete [] ev;
    delete [] store;
    delete [] xef;
    delete [] ef;
    delete [] pdef;

  }

}


void SLGridCyl::mpi_setup(void)
{
				// Get MPI id

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myid);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);


  // Setup for pack/unpack

  int buf1, buf2;
  MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &buf1);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &buf2);

  mpi_bufsz = 2*buf1 +		// m, k
    nmax*buf2 +			// ev
    nmax*numr*buf2 ;		// ef

  mpi_buf = new char [mpi_bufsz];
}


int SLGridCyl::mpi_pack_table(struct TableCyl* table, int m, int k)
{
  int i, j, position = 0;

  MPI_Pack( &m, 1, MPI_INT, mpi_buf, mpi_bufsz, 
	    &position, MPI_COMM_WORLD);
  MPI_Pack( &k, 1, MPI_INT, mpi_buf, mpi_bufsz, 
	    &position, MPI_COMM_WORLD);

  for (j=1; j<=nmax; j++)
    MPI_Pack( &table->ev[j], 1, MPI_DOUBLE, mpi_buf, mpi_bufsz, 
	      &position, MPI_COMM_WORLD);

  for (j=1; j<=nmax; j++)
    for (i=0; i<numr; i++)
      MPI_Pack( &table->ef[j][i], 1, MPI_DOUBLE, mpi_buf, mpi_bufsz, 
		&position, MPI_COMM_WORLD);

  return position;
}


void SLGridCyl::mpi_unpack_table(void)
{
  int i, j, length, position = 0;
  int m, k;

#ifdef DEBUG
  int retid = status.MPI_SOURCE;
#endif

  /*
  MPI_Get_count( &status, MPI_PACKED, &length);
  */
  length = mpi_bufsz;


  MPI_Unpack( mpi_buf, length, &position, &m, 1, MPI_INT,
	      MPI_COMM_WORLD);
  MPI_Unpack( mpi_buf, length, &position, &k, 1, MPI_INT,
	      MPI_COMM_WORLD);

#ifdef DEBUG    
  cout << "Process " << mpi_myid << ": unpacking table entry from Process " 
       << retid << ": (m, k)=(" << m << ", " << k << ")\n";
#endif


  table[m][k].m = m;
  table[m][k].k = k;
  table[m][k].ev.setsize(1, nmax);
  table[m][k].ef.setsize(1, nmax, 0, numr-1);


  for (j=1; j<=nmax; j++)
    MPI_Unpack( mpi_buf, length, &position, &table[m][k].ev[j], 1, MPI_DOUBLE,
		MPI_COMM_WORLD);

  for (j=1; j<=nmax; j++)
    for (i=0; i<numr; i++)
      MPI_Unpack( mpi_buf, length, &position, &table[m][k].ef[j][i], 1, 
		  MPI_DOUBLE, MPI_COMM_WORLD);
}

//======================================================================
//======================================================================
//======================================================================


int SLGridSph::mpi = 0;		// initially off
int SLGridSph::cache = 1;	// initially yes

string SLGridSph::sph_cache_name = ".slgrid_sph_cache";
string SLGridSph::model_file_name = "SLGridSph.model";

extern "C" {
  int sledge_(logical* job, doublereal* cons, logical* endfin, 
	      integer* invec, doublereal* tol, logical* type, 
	      doublereal* ev, integer* numx, doublereal* xef, doublereal* ef, 
	      doublereal* pdef, doublereal* t, doublereal* rho, 
	      integer* iflag, doublereal* store);
  double iv(double, double);
  double kn(int, double);
}


static AxiSymModel *model;

double sphpot(double r)
{
  return model->get_pot(r);
}

double sphdpot(double r)
{
  return model->get_dpot(r);
}

double sphdens(double r)
{
  // This 4pi from Poisson's eqn
  //        |
  //        |       /-- This begins the true density profile
  //        |       |
  //        v       v
  return 4.0*M_PI * model->get_density(r);
}


void SLGridSph::bomb(string oops)
{
  cerr << "SLGridSph [#=" << myid << "]: " << oops << endl; 
  exit(-1);
}

				// Constructors

SLGridSph::SLGridSph(int LMAX, int NMAX, int NUMR,
		     double RMIN, double RMAX, int CMAP, double SCALE,
		     int DIVERGE, double DFAC)
{
  model = new SphericalModelTable(model_file_name, DIVERGE, DFAC);

  initialize(LMAX, NMAX, NUMR, RMIN, RMAX, CMAP, SCALE);
}

SLGridSph::SLGridSph(int LMAX, int NMAX, int NUMR,
		     double RMIN, double RMAX, SphericalModelTable *mod,
		     int CMAP, double SCALE)
{
  model = mod;

  initialize(LMAX, NMAX, NUMR, RMIN, RMAX, CMAP, SCALE);
}


void SLGridSph::initialize(int LMAX, int NMAX, int NUMR,
			   double RMIN, double RMAX, int CMAP, double SCALE)
{
  int l;

  lmax = LMAX;
  nmax = NMAX;
  numr = NUMR;

  rmin = RMIN;
  rmax = RMAX;

  cmap = CMAP;
  scale = SCALE;

  init_table();


#ifdef DEBUG
  if (mpi)
    cout << "Process " << myid << ": MPI is on!\n";
  else
    cout << "Process " << myid << ": MPI is off!\n";
#endif

  if (mpi) {

    table =  new TableSph [lmax+1];

    mpi_setup();

    if (mpi_myid) {
      compute_table_slave();

      //
      // <Receive completed table from master>
      //

      for (l=0; l<=lmax; l++) {

	/*
	MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, 0,
		 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	*/

	MPI_Bcast(mpi_buf, mpi_bufsz, MPI_PACKED, 0, MPI_COMM_WORLD);
    
	mpi_unpack_table();      
      }
      
    }
    else {			// BEGIN Master

      int slave = 0;
      int request_id = 1;

      if (!read_cached_table()) {


	l=0;

	while (l<=lmax) {

	  if (slave<mpi_numprocs-1) { // Send request to slave
	    slave++;
      
#ifdef DEBUG    

	    cout << "Master sending orders to Slave " << slave 
		 << ": l=" << l << '\n'; 
#endif
	    MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	    MPI_Send(&l, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
      
#ifdef DEBUG    
	    cout << "Master gave orders to Slave " << slave 
		 << ": l=" << l << '\n'; 
#endif

				// Increment counters
	    l++;
	  }

	  if (slave == mpi_numprocs-1 && l<=lmax) {
	  
	    //
	    // <Wait and receive>
	    //
	
	    MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 
		     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    
	    int retid = status.MPI_SOURCE;

	    mpi_unpack_table();      

	    //
	    // <Send new request>
	    //

	    MPI_Send(&request_id, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	    MPI_Send(&l, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
      
#ifdef DEBUG    
	    cout << "Master g orders to Slave " << retid
		 << ": l=" << l << '\n'; 
#endif
				// Increment counters
	    l++;
	  }
	}
      
	//
	// <Wait for all slaves to return>
	//
  
	while (slave) {
	
	
	  MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 
		   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
	  mpi_unpack_table();      

	  slave--;
	}

	if (cache) write_cached_table();

      }


      //
      // <Tell slaves to continue>
      //

      request_id = -1;
      for (slave=1; slave < mpi_numprocs; slave++)
	MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);


      //
      // <Send table to slaves>
      //

      for (l=0; l<=lmax; l++) {
	int position = mpi_pack_table(&table[l], l);

	/*
	for (slave=1; slave < mpi_numprocs; slave++)
	  MPI_Send(mpi_buf, position, MPI_PACKED, slave, 11, MPI_COMM_WORLD);
	*/

	MPI_Bcast(mpi_buf, position, MPI_PACKED, 0, MPI_COMM_WORLD);
      }


    } // END Master

  }
  else {

    table =  new TableSph [lmax+1];

    if (!cache || !read_cached_table()) {

      for (l=0; l<=lmax; l++) {
	cerr << "Begin [" << l << "] . . .\n";
	compute_table(&(table[l]), l);
	cerr << ". . . done\n";
      }
      if (cache) write_cached_table();
    }
  }

#ifdef DEBUG
  cerr << "Process " << myid << ": exiting constructor\n";
#endif

}

void check_vector_values_SL(const Vector& v)
{
  for (int i=v.getlow(); i<=v.gethigh(); i++)
    if (isinf(v[i]) || isnan(v[i]))
      {
	cerr << "check_vector: Illegal value\n";
      }
}

int SLGridSph::read_cached_table(void)
{
  if (!cache) return 0;

  ifstream in(sph_cache_name.c_str());
  if (!in) return 0;

  int LMAX, NMAX, NUMR, i, j, CMAP;
  double RMIN, RMAX, SCL;

  if (myid==0) 
    cerr << "SLGridSph::read_cached_table: trying to read cached table . . ."
	 << endl;

  in.read((char *)&LMAX, sizeof(int));		if(!in || LMAX!=lmax) return 0;
  in.read((char *)&NMAX, sizeof(int));		if(!in || NMAX!=nmax) return 0;
  in.read((char *)&NUMR, sizeof(int));		if(!in || NUMR!=numr) return 0;
  in.read((char *)&CMAP, sizeof(int));		if(!in || CMAP!=cmap) return 0;
  in.read((char *)&RMIN, sizeof(double));	if(!in || RMIN!=rmin) return 0;
  in.read((char *)&RMAX, sizeof(double));	if(!in || RMAX!=rmax) return 0;
  in.read((char *)&SCL, sizeof(double));	if(!in || SCL!=scale) return 0;

  for (int l=0; l<=lmax; l++) {

    in.read((char *)&table[l].l, sizeof(int));

				// Double check
    if (table[l].l != l) {
      if (myid==0)
	cerr << "SLGridSph: error reading <" << sph_cache_name << ">" << endl
	     << "SLGridSph: l: read value (" << table[l].l 
	     << ") != internal value (" << l << ")" << endl;
	return 0;
    }

    table[l].ev.setsize(1, nmax);
    table[l].ef.setsize(1, nmax, 0, numr-1);

    for (j=1; j<=nmax; j++) in.read((char *)&table[l].ev[j], sizeof(double));

#ifdef DEBUG_NAN
    check_vector_values_SL(table[l].ev);
#endif

    for (j=1; j<=nmax; j++) {
      for (i=0; i<numr; i++)
	in.read((char *)&table[l].ef[j][i], sizeof(double));
#ifdef DEBUG_NAN
      check_vector_values_SL(table[l].ef[j]);
#endif
    }
  }

  if (myid==0)
    cerr << "SLGridSph::read_cached_table: Success!!\n";

  return 1;
}


void SLGridSph::write_cached_table(void)
{
  ofstream out(sph_cache_name.c_str());
  if (!out) {
    cerr << "SLGridSph: error writing <" << sph_cache_name << ">\n";
    return;
  }

  int i, j;

  out.write((char *)&lmax, sizeof(int));
  out.write((char *)&nmax, sizeof(int));
  out.write((char *)&numr, sizeof(int));
  out.write((char *)&cmap, sizeof(int));
  out.write((char *)&rmin, sizeof(double));
  out.write((char *)&rmax, sizeof(double));
  out.write((char *)&scale, sizeof(double));

  for (int l=0; l<=lmax; l++) {

    out.write((char *)&table[l].l, sizeof(int));

    for (j=1; j<=nmax; j++)
      out.write((char *)&table[l].ev[j], sizeof(double));

    for (j=1; j<=nmax; j++)
      for (i=0; i<numr; i++)
	out.write((char *)&table[l].ef[j][i], sizeof(double));
  }

  cerr << "SLGridSph::write_cached_table: done!!\n";
  return ;
}


SLGridSph::~SLGridSph()
{
  delete [] table;
}

				// Members

double SLGridSph::r_to_xi(double r)
{
  double ret;

  if (cmap==1) {
    if (r<0.0) bomb("radius < 0!");
    ret =  (r/scale-1.0)/(r/scale+1.0);
  } else if (cmap==2) {
    if (r<=0.0) bomb("radius <= 0!");
    ret = log(r);
  } else {
    ret = r;
  }    

  return ret;
}
    
double SLGridSph::xi_to_r(double xi)
{
  double ret;

  if (cmap==1) {
    if (xi<-1.0) bomb("xi < -1!");
    if (xi>=1.0) bomb("xi >= 1!");

    ret =(1.0+xi)/(1.0 - xi) * scale;
  } else if (cmap==2) {
    ret = exp(xi);
  } else {
    if (xi<0.0) bomb("xi < 0!");

    ret = xi;
  }

  return ret;
}

double SLGridSph::d_xi_to_r(double xi)
{
  double ret;

  if (cmap==1) {
    if (xi<-1.0) bomb("xi < -1!");
    if (xi>=1.0) bomb("xi >= 1!");

    ret = 0.5*(1.0-xi)*(1.0-xi)/scale;
  } else if (cmap==2) {
    ret = exp(-xi);
  } else {
    if (xi<0.0) bomb("xi < 0!");
    ret = 1.0;
  }

  return ret;
}

double SLGridSph::get_pot(double x, int l, int n, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }


  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

#ifdef USE_TABLE
  return (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
    sqrt(table[l].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
  return (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
    sqrt(table[l].ev[n]) * sphpot(xi_to_r(x));
#endif
}


double SLGridSph::get_dens(double x, int l, int n, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  
#ifdef USE_TABLE
  return (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1]) *
    sqrt(table[l].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
  return (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1]) *
    sqrt(table[l].ev[n]) * sphdens(xi_to_r(x));
#endif

}

double SLGridSph::get_force(double x, int l, int n, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }


  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  
				// Use three point formula

				// Point -1: indx-1
				// Point  0: indx
				// Point  1: indx+1

  return d_xi_to_r(x)/dxi * (
			     (p - 0.5)*table[l].ef[n][indx-1]*p0[indx-1]
			     -2.0*p*table[l].ef[n][indx]*p0[indx]
			     + (p + 0.5)*table[l].ef[n][indx+1]*p0[indx+1]
			     ) / sqrt(table[l].ev[n]);
}


void SLGridSph::get_pot(Matrix& mat, double x, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  mat.setsize(0, lmax, 1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int l=0; l<=lmax; l++) {
    for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
      mat[l][n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
	sqrt(table[l].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
      mat[l][n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
	sqrt(table[l].ev[n]) * sphpot(xi_to_r(x));
#endif
    }
  }

}


void SLGridSph::get_dens(Matrix& mat, double x, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  mat.setsize(0, lmax, 1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int l=0; l<=lmax; l++) {
    for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
      mat[l][n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])*
	sqrt(table[l].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
      mat[l][n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])*
	sqrt(table[l].ev[n]) * sphdens(xi_to_r(x));
#endif
    }
  }

}


void SLGridSph::get_force(Matrix& mat, double x, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  mat.setsize(0, lmax, 1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int l=0; l<=lmax; l++) {
    for (int n=1; n<=nmax; n++) {
      mat[l][n] = fac * (
			 (p - 0.5)*table[l].ef[n][indx-1]*p0[indx-1]
			 -2.0*p*table[l].ef[n][indx]*p0[indx]
			 + (p + 0.5)*table[l].ef[n][indx+1]*p0[indx+1]
			 ) / sqrt(table[l].ev[n]);
    }
  }
  
}


void SLGridSph::get_pot(Vector& vec, double x, int l, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  vec.setsize(1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
      sqrt(table[l].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
    vec[n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
      sqrt(table[l].ev[n]) * sphpot(xi_to_r(x));
#endif
  }

}


void SLGridSph::get_dens(Vector& vec, double x, int l, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  vec.setsize(1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numr-2) indx = numr - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])*
      sqrt(table[l].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
    vec[n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])*
      sqrt(table[l].ev[n]) * sphdens(xi_to_r(x));
#endif
  }

}


void SLGridSph::get_force(Vector& vec, double x, int l, int which)
{
  if (which || !cmap)
    x = r_to_xi(x);
  else {
    if (cmap==1) {
      if (x<-1.0) x=-1.0;
      if (x>=1.0) x=1.0-XOFFSET;
    }
    if (cmap==2) {
      if (x<xmin) x=xmin;
      if (x>xmax) x=xmax;
    }
  }

  vec.setsize(1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numr-2) indx = numr - 2;


  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int n=1; n<=nmax; n++) {
    vec[n] = fac * (
		    (p - 0.5)*table[l].ef[n][indx-1]*p0[indx-1]
		    -2.0*p*table[l].ef[n][indx]*p0[indx]
		    + (p + 0.5)*table[l].ef[n][indx+1]*p0[indx+1]
		    ) / sqrt(table[l].ev[n]);
  }

}

void SLGridSph::compute_table(struct TableSph* table, int l)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol[6] = {1.0e-4*scale,1.0e-6,  
		   1.0e-4*scale,1.0e-6,  
		   1.0e-4*scale,1.0e-6};
  int i, j, k, VERBOSE=0;
  integer NUM, N;
  logical type[8] = {0, 0, 1, 0, 0, 0, 1, 0};
  logical endfin[2] = {1, 1};
  
#ifdef DEBUG_SLEDGE
  if (myid==1) VERBOSE = SLEDGE_VERBOSE;
#endif

  cons[6] = rmin;
  cons[7] = rmax;
  L2 = l*(l+1);
  NUM = numr;
  N = nmax;
  
  // integer iflag[nmax], invec[nmax+3];
  integer *iflag = new integer [nmax];
  integer *invec = new integer [nmax+3];

  // double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
  // pdef[NUM*N];
  double *t=0, *rho=0;
  double *ev = new double [N];
  double *store = new double [26*(NUM+16)];
  double *xef = new double [NUM+16];
  double *ef = new double [NUM*N];
  double *pdef = new double [NUM*N];
  double f;

				// Inner BC
  f = sphpot(cons[6]);
  if (l==0) {
    cons[0] = sphdpot(cons[6])/f;
    cons[2] = 1.0/(cons[6]*cons[6]*f*f);
  }
  else
    cons[0] = 1.0;

				// Outer BC
  f = sphpot(cons[7]);
  cons[4] = (1.0 + l)/cons[7] + sphdpot(cons[7])/f;
  cons[5] = 1.0/(cons[7]*cons[7]*f*f);

  //
  //     Initialize the vector INVEC(*):
  //       estimates for the eigenvalues/functions specified
  //

  invec[0] = VERBOSE;		// little printing (1), no printing (0)
  invec[1] = 3;			// spectrum is ignored
  invec[2] = N;			// estimates for N eigenvalues/functions

  for (i=0; i<N; i++) invec[3+i] = i;

  //
  //     Set the JOB(*) vector:
  //        estimate both eigenvalues and eigenvectors,
  //        don't estimate the spectral density function,
  //        classify,
  //        let SLEDGE choose the initial mesh
  //
  logical job[5] = {0,1,0,1,0};

  //
  //     Output mesh
  //
  for (i=0; i<NUM; i++) xef[i] = r[i];

  //     
  //     Open file for output.
  //
  sl_dim = 3;

  sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	   t, rho, iflag, store);
  //
  //     Print results:
  //
  
  cout.precision(6);
  cout.setf(ios::scientific);

  for (i=0; i<N; i++) {
    cout << setw(15) << invec[3+i] 
	 << setw(15) << ev[i]
	 << setw(5) << iflag[i]
	 << endl;
  
    if (VERBOSE) {

      if (iflag[i] > -10) {
	cout << setw(14) << "x"
	     << setw(25) << "u(x)"
	     << setw(25) << "(pu`)(x)"
	     << endl;
	k = NUM*i;
	for (j=0; j<NUM; j++) {
	  cout << setw(25) << xef[j]
	       << setw(25) << ef[j+k]
	       << setw(25) << pdef[j+k]
	       << endl;
	}
      }

    }

  }
  
				// Load table

  table->ev.setsize(1, N);
  for (i=0; i<N; i++) table->ev[i+1] = ev[i];

  table->ef.setsize(1, N, 0, numr-1);

  for (i=0; i<numr; i++) {
    for(j=0; j<N; j++) 
      table->ef[j+1][i] = ef[j*NUM+i];
  }

  table->l = l;

  delete [] iflag;
  delete [] invec;
  delete [] ev;
  delete [] store;
  delete [] xef;
  delete [] ef;
  delete [] pdef;

}


void SLGridSph::init_table(void)
{

  int i;

  xi.setsize(0, numr-1);
  r.setsize(0, numr-1);
  p0.setsize(0, numr-1);
  d0.setsize(0, numr-1);

  if (cmap==1) {
    xmin = (rmin/scale - 1.0)/(rmin/scale + 1.0);
    xmax = (rmax/scale - 1.0)/(rmax/scale + 1.0);
  }
  else if (cmap==2) {
    xmin = log(rmin);
    xmax = log(rmax);
  } else {
    xmin = rmin;
    xmax = rmax;
  }
  dxi = (xmax-xmin)/(numr-1);
    
  for (i=0; i<numr; i++) {
    xi[i] = xmin + dxi*i;
    r[i] = xi_to_r(xi[i]);
    p0[i] = sphpot(r[i]);
    d0[i] = sphdens(r[i]);
  }

}


void SLGridSph::compute_table_slave(void)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol[6] = {1.0e-1*scale,1.0e-6,  
		   1.0e-1*scale,1.0e-6,  
		   1.0e-1*scale,1.0e-6};

  int i, j, VERBOSE=0;
  integer NUM;
  logical type[8] = {0, 0, 1, 0, 0, 0, 1, 0};
  logical endfin[2] = {1, 1};
  
  struct TableSph table;
  int L, N;

#ifdef DEBUG_SLEDGE
  if (myid==1) VERBOSE = SLEDGE_VERBOSE;
#endif

#ifdef DEBUG
  cout << "Slave " << mpi_myid << " begins . . .\n";
#endif

  //
  // <Wait for orders>
  //
  
  int request_id;

  while(1) {

    MPI_Recv(&request_id, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (request_id < 0) break;	// Good-bye

    MPI_Recv(&L, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);


#ifdef DEBUG
    cout << "Slave " <<  mpi_myid << ": ordered to compute l = " << L << "\n";
#endif

    cons[0] = cons[1] = cons[2] = cons[3] = cons[4] = cons[5] = 0.0;
    cons[6] = rmin;
    cons[7] = rmax;
    L2 = L*(L+1);
    NUM = numr;
    N = nmax;

    // integer iflag[nmax], invec[nmax+3];
    integer *iflag = new integer [nmax];
    integer *invec = new integer [nmax+3];

    // double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
    // pdef[NUM*N];
    double *t=0, *rho=0;
    double *ev = new double [N];
    double *store = new double [26*(NUM+16)];
    double *xef = new double [NUM+16];
    double *ef = new double [NUM*N];
    double *pdef = new double [NUM*N];
    double f;

    f = sphpot(cons[6]);
    cons[2] = -1.0/(cons[6]*cons[6]*f*f);
    cons[4] = L/cons[7];
    f = sphpot(cons[7]);
    cons[5] = 1.0/(cons[7]*cons[7]*f*f);

    //
    //     Initialize the vector INVEC(*):
    //       estimates for the eigenvalues/functions specified
    //

    invec[0] = VERBOSE;		// little printing (1), no printing (0)
    invec[1] = 3;		// spectrum is ignored
    invec[2] = N;		// estimates for N eigenvalues/functions

    for (i=0; i<N; i++) invec[3+i] = i;

    //
    //     Set the JOB(*) vector:
    //        estimate both eigenvalues and eigenvectors,
    //        don't estimate the spectral density function,
    //        classify,
    //        let SLEDGE choose the initial mesh
    //
    logical job[5] = {0,1,0,1,0};

    //
    //     Output mesh
    //
    for (i=0; i<NUM; i++) xef[i] = r[i];

    //     
    //     Open file for output.
    //
    sl_dim = 3;

    sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	    t, rho, iflag, store);
    //
    //     Print results:
    //
    
#ifdef DEBUG
    cout << "Slave " <<  mpi_myid << ": computed l = " << L << "\n";

    cout.precision(6);
    cout.setf(ios::scientific);

    for (i=0; i<N; i++) {
      cout << setw(15) << invec[3+i] 
	   << setw(15) << ev[i]
	   << setw(5) << iflag[i]
	   << endl;
  
      if (VERBOSE) {

	if (iflag[i] > -10) {
	  cout << setw(14) << "x"
	       << setw(25) << "u(x)"
	       << setw(25) << "(pu`)(x)"
	       << endl;
	  int k = NUM*i;
	  for (j=0; j<NUM; j++) {
	    cout << setw(25) << xef[j]
		 << setw(25) << ef[j+k]
		 << setw(25) << pdef[j+k]
		 << endl;
	  }
	}
	
      }

    }
  
#endif
				// Load table

    table.ev.setsize(1, N);
    for (i=0; i<N; i++) table.ev[i+1] = ev[i];

    table.ef.setsize(1, N, 0, numr-1);
    for (i=0; i<numr; i++) {
      for(j=0; j<N; j++) 
	table.ef[j+1][i] = ef[j*NUM+i];
    }

    table.l = L;
  
    int position = mpi_pack_table(&table, L);
    MPI_Send(mpi_buf, position, MPI_PACKED, 0, 11, MPI_COMM_WORLD);

#ifdef DEBUG
    cout << "Slave " <<  mpi_myid << ": send to master l = " << L << "\n";
#endif

    delete [] iflag;
    delete [] invec;
    delete [] ev;
    delete [] store;
    delete [] xef;
    delete [] ef;
    delete [] pdef;

  }

}


void SLGridSph::mpi_setup(void)
{
				// Get MPI id

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myid);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);


  // Setup for pack/unpack

  int buf1, buf2;
  MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &buf1);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &buf2);

  mpi_bufsz = buf1 +		// l
    nmax*buf2 +			// ev
    nmax*numr*buf2 ;		// ef

  mpi_buf = new char [mpi_bufsz];
}


int SLGridSph::mpi_pack_table(struct TableSph* table, int l)
{
  int i, j, position = 0;

  MPI_Pack( &l, 1, MPI_INT, mpi_buf, mpi_bufsz, 
	    &position, MPI_COMM_WORLD);

  for (j=1; j<=nmax; j++)
    MPI_Pack( &table->ev[j], 1, MPI_DOUBLE, mpi_buf, mpi_bufsz, 
	      &position, MPI_COMM_WORLD);

  for (j=1; j<=nmax; j++)
    for (i=0; i<numr; i++)
      MPI_Pack( &table->ef[j][i], 1, MPI_DOUBLE, mpi_buf, mpi_bufsz, 
		&position, MPI_COMM_WORLD);

  return position;
}


void SLGridSph::mpi_unpack_table(void)
{
  int i, j, length, position = 0;
  int l;

  /*
  MPI_Get_count( &status, MPI_PACKED, &length);
  */
  length = mpi_bufsz;

#ifdef DEBUG
  int retid = status.MPI_SOURCE;
#endif

  MPI_Unpack( mpi_buf, length, &position, &l, 1, MPI_INT,
	      MPI_COMM_WORLD);

#ifdef DEBUG    
  cout << "Process " <<  mpi_myid << ": unpacking table entry from Process " 
       << retid << "  l = " << l << "\n";
#endif


  table[l].l = l;
  table[l].ev.setsize(1, nmax);
  table[l].ef.setsize(1, nmax, 0, numr-1);

  for (j=1; j<=nmax; j++)
    MPI_Unpack( mpi_buf, length, &position, &table[l].ev[j], 1, MPI_DOUBLE,
		MPI_COMM_WORLD);

  for (j=1; j<=nmax; j++)
    for (i=0; i<numr; i++)
      MPI_Unpack( mpi_buf, length, &position, &table[l].ef[j][i], 1, 
		  MPI_DOUBLE, MPI_COMM_WORLD);
}


//======================================================================
//======================================================================
//======================================================================


int SLGridSlab::mpi = 0;	// initially off
int SLGridSlab::cache = 1;	// initially yes
double SLGridSlab::H = 0.1;	// Scale height
double SLGridSlab::L = 1.0;	// Periodic box size
double SLGridSlab::ZBEG = 0.0;	// Offset on from origin
double SLGridSlab::ZEND = 0.0;	// Offset on potential zero

static double KKZ;

				// Isothermal slab with
				// G = M = 1
static double poffset=0.0;

double slabpot(double z)
{
  return 2.0*M_PI*SLGridSlab::H*log(cosh(z/SLGridSlab::H)) - poffset;
}

double slabdpot(double z)
{
  return 2.0*M_PI*tanh(z/SLGridSlab::H);
}

double slabdens(double z)
{
  double tmp = 1.0/cosh(z/SLGridSlab::H);
  return 4.0*M_PI * 0.5/SLGridSlab::H * tmp*tmp;
}


void SLGridSlab::bomb(string oops)
{
  cerr << "SLGridSlab: " << oops << endl; 
  exit(-1);
}

				// Constructors

SLGridSlab::SLGridSlab(int NUMK, int NMAX, int NUMZ, double ZMAX)
{
  int kx, ky;

  numk = NUMK;
  nmax = NMAX;
  numz = NUMZ;

  zmax = ZMAX;

  poffset = 0.0;
  poffset = slabpot((1.0+ZEND)*zmax);

  init_table();


#ifdef DEBUG
  if (mpi)
    cout << "Process " << myid << ": MPI is on!\n";
  else
    cout << "Process " << myid << ": MPI is off!\n";
#endif

  table =  new TableSlab* [numk+1];
  for (kx=0; kx<=numk; kx++)
    table[kx] =  new TableSlab [kx+1];

  if (mpi) {

    mpi_setup();

    if (mpi_myid) {
      compute_table_slave();

      //
      // <Receive completed table from master>
      //

      for (kx=0; kx<=numk; kx++) {
	for (ky=0; ky<=kx; ky++) {
	  MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, 0,
		   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
	  mpi_unpack_table();      
	}
      }
      
    }
    else {			// BEGIN Master

      int slave = 0;
      int request_id = 1;

      if (!read_cached_table()) {

	kx = 0;
	ky = 0;

	while (kx<=numk) {

	  if (slave<mpi_numprocs-1) { // Send request to slave
	    slave++;
      
#ifdef DEBUG    
	    cout << "Master sending orders to Slave " << slave 
		 << ": Kx=" << kx << ", Ky=" << ky << '\n';
#endif
	    MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	    MPI_Send(&kx, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	    MPI_Send(&ky, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
      
#ifdef DEBUG    
	    cout << "Master gave orders to Slave " << slave 
		 << ": Kx=" << kx << ", Ky=" << ky << '\n';
#endif

				// Increment counters
	    ky++;
	    if (ky>kx) {
	      kx++;
	      ky = 0;
	    }
	    
	  }

	  if (slave == mpi_numprocs-1 && kx<=numk) {
	  
	    //
	    // <Wait and receive>
	    //
	
	    MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 
		     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    
	    int retid = status.MPI_SOURCE;

	    mpi_unpack_table();      

	    //
	    // <Send new request>
	    //

	    MPI_Send(&request_id, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	    MPI_Send(&kx, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	    MPI_Send(&ky, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
      
#ifdef DEBUG    
	    cout << "Master gave orders to Slave " << retid 
		 << ": Kx=" << kx << ", Ky=" << ky << '\n';
#endif
				// Increment counters
	    ky++;
	    if (ky>kx) {
	      kx++;
	      ky = 0;
	    }

	  }
	}
      
	//
	// <Wait for all slaves to return>
	//
  
	while (slave) {
	
	
	  MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, 
		   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
	  mpi_unpack_table();      

	  slave--;
	}

	if (cache) write_cached_table();

      }


      //
      // <Tell slaves to continue>
      //

      request_id = -1;
      for (slave=1; slave < mpi_numprocs; slave++)
	MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);


      //
      // <Send table to slaves>
      //

      for (kx=0; kx<=numk; kx++) {
	for (ky=0; ky<=kx; ky++) {
	  int position = mpi_pack_table(&table[kx][ky], kx, ky);
	  for (slave=1; slave < mpi_numprocs; slave++)
	    MPI_Send(mpi_buf, position, MPI_PACKED, slave, 11, MPI_COMM_WORLD);
	}
      }


    } // END Master

  }
  else {

    if (!read_cached_table()) {

      for (kx=0; kx<=numk; kx++) {
	for (ky=0; ky<=kx; ky++) {
	  cerr << "Begin [" << kx << ", " << ky << "] . . .\n";
	  compute_table(&(table[kx][ky]), kx, ky);
	  cerr << ". . . done\n";
	}
      }

      if (cache) write_cached_table();
    }
  }

#ifdef DEBUG
  cerr << "Process " << myid << ": exiting constructor\n";
#endif
}


const string Slab_cache_name = ".slgrid_slab_cache";

int SLGridSlab::read_cached_table(void)
{
  if (!cache) return 0;

  ifstream in(Slab_cache_name.c_str());
  if (!in) return 0;

  int NUMK, NMAX, NUMZ, i, j;
  double ZMAX, HH, LL, zbeg, zend;

  if (myid==0)
    cerr << "SLGridSlab::read_cached_table: trying to read cached table . . ."
	 << endl;

  in.read((char *)&NUMK, sizeof(int));		if(!in || NUMK!=numk) return 0;
  in.read((char *)&NMAX, sizeof(int));		if(!in || NMAX!=nmax) return 0;
  in.read((char *)&NUMZ, sizeof(int));		if(!in || NUMZ!=numz) return 0;
  in.read((char *)&ZMAX, sizeof(double));	if(!in || ZMAX!=zmax) return 0;
  in.read((char *)&HH, sizeof(double));		if(!in || HH!=H)      return 0;
  in.read((char *)&LL, sizeof(double));		if(!in || LL!=L)      return 0;
  in.read((char *)&zbeg, sizeof(double));	if(!in || zbeg!=ZBEG) return 0;
  in.read((char *)&zend, sizeof(double));	if(!in || zend!=ZEND) return 0;

  for (int kx=0; kx<=numk; kx++) {
    for (int ky=0; ky<=kx; ky++) {

      in.read((char *)&table[kx][ky].kx, sizeof(int));
      in.read((char *)&table[kx][ky].ky, sizeof(int));

				// Double check
      if (table[kx][ky].kx != kx) {
	if (myid==0)
	  cerr << "SLGridSlab: error reading <" << Slab_cache_name << ">"
	       << endl
	       << "SLGridSlab: kx: read value (" << table[kx][ky].kx 
	       << ") != internal value (" << kx << ")" << endl;
	return 0;
      }
      if (table[kx][ky].ky != ky) {
	if (myid==0) 
	  cerr << "SLGridSlab: error reading <" << Slab_cache_name << ">"
	       << endl
	       << "SLGridSlab: ky: read value (" << table[kx][ky].ky 
	       << ") != internal value (" << ky << ")" << endl;
	return 0;
      }

      table[kx][ky].ev.setsize(1, nmax);
      table[kx][ky].ef.setsize(1, nmax, 0, numz-1);

      for (j=1; j<=nmax; j++) in.read((char *)&table[kx][ky].ev[j], sizeof(double));
      
#ifdef DEBUG_NAN
      check_vector_values_SL(table[kx][ky].ev);
#endif

      for (j=1; j<=nmax; j++) {
	for (i=0; i<numz; i++)
	  in.read((char *)&table[kx][ky].ef[j][i], sizeof(double));
#ifdef DEBUG_NAN
	check_vector_values_SL(table[kx][ky].ef[j]);
#endif
      }
    }
  }

  if (myid==0)
    cerr << "SLGridSlab::read_cached_table: Success!!\n";

  return 1;
}


void SLGridSlab::write_cached_table(void)
{
  ofstream out(Slab_cache_name.c_str());
  if (!out) {
    cerr << "SLGridSlab: error writing <" << Slab_cache_name << ">\n";
    return;
  }

  int i, j;

  out.write((char *)&numk, sizeof(int));
  out.write((char *)&nmax, sizeof(int));
  out.write((char *)&numz, sizeof(int));
  out.write((char *)&zmax, sizeof(double));
  out.write((char *)&H,    sizeof(double));
  out.write((char *)&L,    sizeof(double));
  out.write((char *)&ZBEG, sizeof(double));
  out.write((char *)&ZEND, sizeof(double));

  for (int kx=0; kx<=numk; kx++) {
    for (int ky=0; ky<=kx; ky++) {

      out.write((char *)&table[kx][ky].kx, sizeof(int));
      out.write((char *)&table[kx][ky].ky, sizeof(int));

      for (j=1; j<=nmax; j++)
	out.write((char *)&table[kx][ky].ev[j], sizeof(double));

      for (j=1; j<=nmax; j++)
	for (i=0; i<numz; i++)
	  out.write((char *)&table[kx][ky].ef[j][i], sizeof(double));
    }
  }

  cerr << "SLGridSlab::write_cached_table: done!!\n";
  return ;
}


SLGridSlab::~SLGridSlab()
{
  for (int kx=0; kx<=numk; kx++) delete [] table[kx];
  delete [] table;
}

				// Members

/*
double SLGridSlab::z_to_xi(double z)
{
  return tanh(z/H);
}
    
double SLGridSlab::xi_to_z(double xi)
{
  return H*atanh(xi);
}

double SLGridSlab::d_xi_to_z(double xi)
{
  return H/(1.0 - xi*xi);
}
*/


/*
double SLGridSlab::z_to_xi(double z)
{
  return z/sqrt(z*z + H*H);
}
    
double SLGridSlab::xi_to_z(double xi)
{
  return xi*H/sqrt(1.0 - xi*xi);
}

double SLGridSlab::d_xi_to_z(double xi)
{
  return pow(1.0 - xi*xi, 1.5)/H;
}
*/

				// Simple cartesian coordinates seem
				// to work best here; this transformation
				// is the identity . . . 

double SLGridSlab::z_to_xi(double z)
{
  return z;
}

double SLGridSlab::xi_to_z(double xi)
{
  return xi;
}

double SLGridSlab::d_xi_to_z(double xi)
{
  return 1.0;
}



double SLGridSlab::get_pot(double x, int kx, int ky, int n, int which)
{
  int hold;

				// Flip sign for antisymmetric basis functions
  int sign=1;
  if (x<0 && 2*(n/2)==n) sign=-1;
  x = fabs(x);

  if (which)
    x = z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

#ifdef USE_TABLE
  return (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])/
    sqrt(table[kx][ky].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]) * sign;
#else
  return (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])/
    sqrt(table[kx][ky].ev[n]) * slabpot(xi_to_z(x)) * sign;
#endif
}


double SLGridSlab::get_dens(double x, int kx, int ky, int n, int which)
{
  int hold;

  int sign=1;
  if (x<0 && 2*(n/2)==n) sign=-1;
  x = fabs(x);
  
  if (which)
    x = z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  
#ifdef USE_TABLE
  return (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1]) *
    sqrt(table[kx][ky].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]) * sign;
#else
  return (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1]) *
    sqrt(table[kx][ky].ev[n]) * slabdens(xi_to_z(x)) * sign;
#endif

}

double SLGridSlab::get_force(double x, int kx, int ky, int n, int which)
{
  int hold;

  int sign=1;
  if (x<0 && 2*(n/2)!=n) sign = -1;
  x = fabs(x);

  if (which)
    x = z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numz-2) indx = numz - 2;


  double p = (x - xi[indx])/dxi;
  
				// Use three point formula

				// Point -1: indx-1
				// Point  0: indx
				// Point  1: indx+1

  return d_xi_to_z(x)/dxi * (
			     (p - 0.5)*table[kx][ky].ef[n][indx-1]*p0[indx-1]
			     -2.0*p*table[kx][ky].ef[n][indx]*p0[indx]
			     + (p + 0.5)*table[kx][ky].ef[n][indx+1]*p0[indx+1]
			     ) / sqrt(table[kx][ky].ev[n]) * sign;
}


void SLGridSlab::get_pot(Matrix& mat, double x, int which)
{
  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)
    x = z_to_xi(x);

  int ktot = (numk+1)*(numk+2)/2;
  mat.setsize(0, ktot, 1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;


  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  int l=0;
  for (int kx=0; kx<=numk; kx++) {
    for (int ky=0; ky<=kx; ky++) {
      sign2 = 1;
      for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
	mat[l][n] = (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])/
	  sqrt(table[kx][ky].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]) * sign2;
#else
	mat[l][n] = (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])/
	  sqrt(table[kx][ky].ev[n]) * slabpot(xi_to_z(x)) * sign2;
#endif
	sign2 *= sign;
      }
      l++;
    }    
  }

}


void SLGridSlab::get_dens(Matrix& mat, double x, int which)
{
  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)
    x = z_to_xi(x);

  int ktot = (numk+1)*(numk+2)/2;
  mat.setsize(0, ktot, 1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  
  int l=0;
  for (int kx=0; kx<=numk; kx++) {
    for (int ky=0; ky<=kx; ky++) {
      sign2 = 1;
      for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
	mat[l][n] = (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])*
	  sqrt(table[kx][ky].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]) * sign2;
#else
	mat[l][n] = (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])*
	  sqrt(table[kx][ky].ev[n]) * slabdens(xi_to_z(x)) * sign2;
#endif
	sign2 *= sign;
      }
      l++;
    }
  }

}


void SLGridSlab::get_force(Matrix& mat, double x, int which)
{
  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)
    x = z_to_xi(x);

  int ktot = (numk+1)*(numk+2)/2;
  mat.setsize(0, ktot, 1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numz-2) indx = numz - 2;


  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_z(x)/dxi;

  int l=0;
  for (int kx=0; kx<=numk; kx++) {
    for (int ky=0; ky<=kx; ky++) {
      sign2 = sign;
      for (int n=1; n<=nmax; n++) {
	mat[l][n] = fac * (
			   (p - 0.5)*table[kx][ky].ef[n][indx-1]*p0[indx-1]
			   -2.0*p*table[kx][ky].ef[n][indx]*p0[indx]
			   + (p + 0.5)*table[kx][ky].ef[n][indx+1]*p0[indx+1]
			   ) / sqrt(table[kx][ky].ev[n]) * sign2;
	sign2 *= -sign;
      }
      l++;
    }
  }
  
}


void SLGridSlab::get_pot(Vector& vec, double x, int kx, int ky, int which)
{
  int hold;

  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)
    x = z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  vec.setsize(1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  sign2 = 1;
  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])/
      sqrt(table[kx][ky].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]) * sign2;
#else
    vec[n] = (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])/
      sqrt(table[kx][ky].ev[n]) * slabpot(xi_to_z(x)) * sign2;
#endif
    sign2 *= sign;
  }

}


void SLGridSlab::get_dens(Vector& vec, double x, int kx, int ky, int which)
{
  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)
    x = z_to_xi(x);

  vec.setsize(1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<0) indx = 0;
  if (indx>numz-2) indx = numz - 2;

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  sign2 = 1;
  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])*
      sqrt(table[kx][ky].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]) * sign2;
#else
    vec[n] = (x1*table[kx][ky].ef[n][indx] + x2*table[kx][ky].ef[n][indx+1])*
      sqrt(table[kx][ky].ev[n]) * slabdens(xi_to_z(x)) * sign2;
#endif
    sign2 *= sign;
  }

}


void SLGridSlab::get_force(Vector& vec, double x, int kx, int ky, int which)
{
  int hold;

  int sign=1, sign2;
  if (x<0) sign = -1;
  x = fabs(x);
  
  if (which)
    x = z_to_xi(x);

  if (ky > kx) {
    hold = ky;
    ky = kx;
    kx = hold;
  }

  vec.setsize(1, nmax);

  int indx = (int)( (x-xmin)/dxi );
  if (indx<1) indx = 1;
  if (indx>numz-2) indx = numz - 2;


  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_z(x)/dxi;

  sign2 = sign;
  for (int n=1; n<=nmax; n++) {
    vec[n] = fac * (
		    (p - 0.5)*table[kx][ky].ef[n][indx-1]*p0[indx-1]
		    -2.0*p*table[kx][ky].ef[n][indx]*p0[indx]
		    + (p + 0.5)*table[kx][ky].ef[n][indx+1]*p0[indx+1]
		    ) / sqrt(table[kx][ky].ev[n]) * sign2;
    sign2 *= sign;
  }

}

void SLGridSlab::compute_table(struct TableSlab* table, int KX, int KY)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol[6] = {1.0e-4,1.0e-5,  1.0e-4,1.0e-5,  1.0e-4,1.0e-5};
  int i, j, k, VERBOSE=0;
  integer NUM, N;
  logical type[8] = {1, 0, 0, 0, 1, 0, 0, 0};
  logical endfin[2] = {1, 1};
  
#ifdef DEBUG_SLEDGE
  if (myid==1) VERBOSE = SLEDGE_VERBOSE;
#endif

  cons[6] =  ZBEG;
  cons[7] =  zmax;
  NUM = numz;
				// Divide total functions into symmetric
				// and antisymmetric (keeping equal number
				// of each or one more symmetric
  N = (int)( 0.5*nmax + 0.501);
  
  integer *iflag = new integer [nmax];
  integer *invec = new integer [nmax+3];

  // double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
  // pdef[NUM*N];
  double *t=0, *rho=0;
  double *ev = new double [N];
  double *store = new double [26*(NUM+16)];
  double *xef = new double [NUM+16];
  double *ef = new double [NUM*N];
  double *pdef = new double [NUM*N];
  double f, df;

  KKZ = 2.0*M_PI/L * sqrt((double)(KX*KX + KY*KY));

				// Even BC, inner has zero gradient
  f = slabpot(cons[6]);
  cons[2] = -1.0/(f*f);

				// Outer
  if (KKZ>1.0e-4) {
    f = slabpot(cons[7]);
    df = slabdpot(cons[7]);
    cons[4] = (df + KKZ*f)*f;
  }
  cons[5] = 1.0;
  //  cons[5] = 1.0/(f*f);

  //
  //     Initialize the vector INVEC(*):
  //       estimates for the eigenvalues/functions specified
  //

  invec[0] = VERBOSE;		// little printing (1), no printing (0)
  invec[1] = 3;			// spectrum is ignored
  invec[2] = N;			// estimates for N eigenvalues/functions

  for (i=0; i<N; i++) invec[3+i] = i;

  //
  //     Set the JOB(*) vector:
  //        estimate both eigenvalues and eigenvectors,
  //        don't estimate the spectral density function,
  //        classify,
  //        let SLEDGE choose the initial mesh
  //
  logical job[5] = {0,1,0,1,0};

  //
  //     Output mesh
  //
  for (i=0; i<NUM; i++) xef[i] = z[i];

  //     
  //     Open file for output.
  //
  sl_dim = 1;

  sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	   t, rho, iflag, store);
  //
  //     Print results:
  //
  
  cout.precision(6);
  cout.setf(ios::scientific);

  cout << "Even:\n";
  for (i=0; i<N; i++) {
    cout << setw(15) << invec[3+i] 
	 << setw(15) << ev[i]
	 << setw(5) << iflag[i]
	 << endl;
  
    if (VERBOSE) {

      if (iflag[i] > -10) {
	cout << setw(14) << "x"
	     << setw(25) << "u(x)"
	     << setw(25) << "(pu`)(x)"
	     << endl;
	k = NUM*i;
	for (j=0; j<NUM; j++) {
	  cout << setw(25) << xef[j]
	       << setw(25) << ef[j+k]
	       << setw(25) << pdef[j+k]
	       << endl;
	}
      }

    }

  }
  
				// Load table

  table->ev.setsize(1, nmax);
  for (i=0; i<N; i++) table->ev[i*2+1] = ev[i];

  table->ef.setsize(1, nmax, 0, numz-1);
  for (i=0; i<numz; i++) {
    for(j=0; j<N; j++) 
      table->ef[j*2+1][i] = ef[j*NUM+i];
  }


				// Odd BC, Inner zero value
  cons[0] = 1.0;
  cons[2] = 0.0;

				// Redo to get antisymmetric functions

  sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	  t, rho, iflag, store);


  //
  //     Print results:
  //
  
  cout.precision(6);
  cout.setf(ios::scientific);

  cout << "Odd:\n";
  for (i=0; i<N; i++) {
    cout << setw(15) << invec[3+i] 
	 << setw(15) << ev[i]
	 << setw(5) << iflag[i]
	 << endl;
  
    if (VERBOSE) {

      if (iflag[i] > -10) {
	cout << setw(14) << "x"
	     << setw(25) << "u(x)"
	     << setw(25) << "(pu`)(x)"
	     << endl;
	k = NUM*i;
	for (j=0; j<NUM; j++) {
	  cout << setw(25) << xef[j]
	       << setw(25) << ef[j+k]
	       << setw(25) << pdef[j+k]
	       << endl;
	}
      }

    }

  }

				// Load table

  N = nmax - N;

  for (i=0; i<N; i++) table->ev[i*2+2] = ev[i];

  for (i=0; i<numz; i++) {
    for(j=0; j<N; j++) 
      table->ef[j*2+2][i] = ef[j*NUM+i];
  }

				// Correct for symmetrizing
  table->ef *= 7.071067811865475e-01;

  table->kx = KX;
  table->ky = KY;

  delete [] iflag;
  delete [] invec;
  delete [] ev;
  delete [] store;
  delete [] xef;
  delete [] ef;
  delete [] pdef;
}


void SLGridSlab::init_table(void)
{

  int i;

  xi.setsize(0, numz-1);
  z.setsize(0, numz-1);
  p0.setsize(0, numz-1);
  d0.setsize(0, numz-1);

  xmin = z_to_xi( ZBEG);
  xmax = z_to_xi( zmax);
  dxi = (xmax-xmin)/(numz-1);


  for (i=0; i<numz; i++) {
    xi[i] = xmin + dxi*i;
    z[i] = xi_to_z(xi[i]);
    p0[i] = slabpot(z[i]);
    d0[i] = slabdens(z[i]);
  }

}


void SLGridSlab::compute_table_slave(void)
{

  //  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  //  double tol[6] = {1.0e-4,1.0e-5,  1.0e-4,1.0e-5,  1.0e-4,1.0e-5};
  double cons[8];
  double tol[6] = {1.0e-6,1.0e-7,  1.0e-6,1.0e-7,  1.0e-6,1.0e-7};
  int i, j, VERBOSE=0;
  integer NUM;
  logical type[8] = {1, 0, 0, 0, 1, 0, 0, 0};
  logical endfin[2] = {1, 1};
  
  struct TableSlab table;
  int KX, KY, N;

#ifdef DEBUG_SLEDGE
  if (myid==1) VERBOSE = SLEDGE_VERBOSE;
#endif

#ifdef DEBUG
  cout << "Slave " << mpi_myid << " begins . . .\n";
#endif

  //
  // <Wait for orders>
  //
  
  int request_id;

  while(1) {

    MPI_Recv(&request_id, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (request_id < 0) break;	// Good-bye

    MPI_Recv(&KX, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    MPI_Recv(&KY, 1, 
	     MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);


#ifdef DEBUG
    cout << "Slave " << mpi_myid << ": ordered to compute Kx, Ky = "
	 << KX << ", " << KY << "\n";
#endif

    cons[0] = cons[1] = cons[2] = cons[3] = cons[4] = cons[5] = 0.0;
    cons[6] =  ZBEG;
    cons[7] =  zmax;
    NUM = numz;
				// Divide total functions into symmetric
				// and antisymmetric (keeping equal number
				// of each or one more symmetric
    N = (int)( 0.5*nmax + 0.501);

    // integer iflag[nmax], invec[nmax+3];
    integer *iflag = new integer [nmax];
    integer *invec = new integer [nmax+3];

    // double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
    // pdef[NUM*N];
    double *t=0, *rho=0;
    double *ev = new double [N];
    double *store = new double [26*(NUM+16)];
    double *xef = new double [NUM+16];
    double *ef = new double [NUM*N];
    double *pdef = new double [NUM*N];
    double f, df;

    KKZ = 2.0*M_PI/L * sqrt((double)(KX*KX + KY*KY));

				// Even BC, inner has zero gradient
    f = slabpot(cons[6]);
    cons[2] = -1.0/(f*f);

#ifdef DEBUG
    cout << "Slave " << mpi_myid << ": Kx, Ky = " 
	 << KX << ", " << KY << " [Even inputs]\n";
    for (int n=0; n<8; n++)
      cout << setw(5) << n << setw(15) << cons[n] << endl;
#endif

				// Outer
    if (KKZ>1.0e-4) {
      f = slabpot(cons[7]);
      df = slabdpot(cons[7]);
      cons[4] = (df + KKZ*f)*f;
    }
    cons[5] = 1.0;

    //
    //     Initialize the vector INVEC(*):
    //       estimates for the eigenvalues/functions specified
    //

    invec[0] = VERBOSE;		// little printing (1), no printing (0)
    invec[1] = 3;		// spectrum is ignored
    invec[2] = N;		// estimates for N eigenvalues/functions

    for (i=0; i<N; i++) invec[3+i] = i;

    //
    //     Set the JOB(*) vector:
    //        estimate both eigenvalues and eigenvectors,
    //        don't estimate the spectral density function,
    //        don't classify,
    //        we choose the initial mesh
    //
    logical job[5] = {0,1,0,1,0};

    //
    //     Output mesh
    //
    for (i=0; i<NUM; i++) xef[i] = z[i];

    //     
    //     Open file for output.
    //
    sl_dim = 1;

    sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	    t, rho, iflag, store);
    //
    //     Print results:
    //
    
#ifdef DEBUG
    cout << "Slave " << mpi_myid << ": computed Kx, Ky = " 
	 << KX << ", " << KY << " [Even]\n";

    cout.precision(6);
    cout.setf(ios::scientific);

    for (i=0; i<N; i++) {
      cout << setw(15) << invec[3+i] 
	   << setw(15) << ev[i]
	   << setw(5) << iflag[i]
	   << endl;
  
      if (VERBOSE) {

	if (iflag[i] > -10) {
	  cout << setw(14) << "x"
	       << setw(25) << "u(x)"
	       << setw(25) << "(pu`)(x)"
	       << endl;
	  int k = NUM*i;
	  for (j=0; j<NUM; j++) {
	    cout << setw(25) << xef[j]
		 << setw(25) << ef[j+k]
		 << setw(25) << pdef[j+k]
		 << endl;
	  }
	}
	
      }

    }
  
#endif
				// Load table

    table.ev.setsize(1, nmax);
    for (i=0; i<N; i++) table.ev[i*2+1] = ev[i];

    table.ef.setsize(1, nmax, 0, numz-1);
    for (i=0; i<numz; i++) {
      for(j=0; j<N; j++) 
	table.ef[j*2+1][i] = ef[j*NUM+i];
    }


				// Odd BC, inner zero value
    cons[0] = 1.0;
    cons[2] = 0.0;

				// Redo to get antisymmetric functions

    sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	    t, rho, iflag, store);
    //
    //     Print results:
    //
    
#ifdef DEBUG
    cout << "Slave " << mpi_myid << ": computed Kx, Ky = " 
	 << KX << ", " << KY << " [Odd]\n";

    cout.precision(6);
    cout.setf(ios::scientific);

    for (i=0; i<N; i++) {
      cout << setw(15) << invec[3+i] 
	   << setw(15) << ev[i]
	   << setw(5) << iflag[i]
	   << endl;
  
      if (VERBOSE) {

	if (iflag[i] > -10) {
	  cout << setw(14) << "x"
	       << setw(25) << "u(x)"
	       << setw(25) << "(pu`)(x)"
	       << endl;
	  int k = NUM*i;
	  for (j=0; j<NUM; j++) {
	    cout << setw(25) << xef[j]
		 << setw(25) << ef[j+k]
		 << setw(25) << pdef[j+k]
		 << endl;
	  }
	}
	
      }

    }
  
#endif
				// Load table

    N = nmax - N;

    for (i=0; i<N; i++) table.ev[i*2+2] = ev[i];

    for (i=0; i<numz; i++) {
      for(j=0; j<N; j++) 
	table.ef[j*2+2][i] = ef[j*NUM+i];
    }

				// Correct for symmetrizing
    table.ef *= 7.071067811865475e-01;
    table.kx = KX;
    table.ky = KY;
  
    int position = mpi_pack_table(&table, KX, KY);
    MPI_Send(mpi_buf, position, MPI_PACKED, 0, 11, MPI_COMM_WORLD);

#ifdef DEBUG
    cout << "Slave " << mpi_myid << ": sent to master Kx, Ky = "
	 << KX << ", " <<  KY << '\n';
#endif

    delete [] iflag;
    delete [] invec;
    delete [] ev;
    delete [] store;
    delete [] xef;
    delete [] ef;
    delete [] pdef;

  }

}


void SLGridSlab::mpi_setup(void)
{
				// Get MPI id

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myid);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);


  // Setup for pack/unpack

  int buf1, buf2;
  MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &buf1);
  MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &buf2);

  mpi_bufsz = 2*buf1 +		// kx, ky
    nmax*buf2 +			// ev
    nmax*numz*buf2 ;		// ef

  mpi_buf = new char [mpi_bufsz];
}


int SLGridSlab::mpi_pack_table(struct TableSlab* table, int kx, int ky)
{
  int i, j, position = 0;

  MPI_Pack( &kx, 1, MPI_INT, mpi_buf, mpi_bufsz, 
	    &position, MPI_COMM_WORLD);

  MPI_Pack( &ky, 1, MPI_INT, mpi_buf, mpi_bufsz, 
	    &position, MPI_COMM_WORLD);

  for (j=1; j<=nmax; j++)
    MPI_Pack( &table->ev[j], 1, MPI_DOUBLE, mpi_buf, mpi_bufsz, 
	      &position, MPI_COMM_WORLD);

  for (j=1; j<=nmax; j++)
    for (i=0; i<numz; i++) {
      MPI_Pack( &table->ef[j][i], 1, MPI_DOUBLE, mpi_buf, mpi_bufsz, 
		&position, MPI_COMM_WORLD);
    }

  return position;
}


void SLGridSlab::mpi_unpack_table(void)
{
  int i, j, length, position = 0;
  int kx, ky;

#ifdef DEBUG
  int retid = status.MPI_SOURCE;
#endif

  MPI_Get_count( &status, MPI_PACKED, &length);


  MPI_Unpack( mpi_buf, length, &position, &kx, 1, MPI_INT,
	      MPI_COMM_WORLD);

  MPI_Unpack( mpi_buf, length, &position, &ky, 1, MPI_INT,
	      MPI_COMM_WORLD);

#ifdef DEBUG    
  cout << "Process " << mpi_myid << ": unpacking table entry from Process "
       << retid << ": kx=" << kx << ", " << ky << "\n";
#endif


  table[kx][ky].kx = kx;
  table[kx][ky].ky = ky;
  table[kx][ky].ev.setsize(1, nmax);
  table[kx][ky].ef.setsize(1, nmax, 0, numz-1);

  for (j=1; j<=nmax; j++)
    MPI_Unpack( mpi_buf, length, &position, &table[kx][ky].ev[j], 1, MPI_DOUBLE,
		MPI_COMM_WORLD);

  for (j=1; j<=nmax; j++)
    for (i=0; i<numz; i++)
      MPI_Unpack( mpi_buf, length, &position, &table[kx][ky].ef[j][i], 1, 
		  MPI_DOUBLE, MPI_COMM_WORLD);
}






//======================================================================
//======================================================================
//======================================================================


extern "C" int coeff_(doublereal* x, doublereal* px, doublereal* qx, 
	doublereal* rx)
{
  double f,rho;

  if (sl_dim==1) {		// 1-d slab

    f = slabpot(*x);
    rho = slabdens(*x);

    *px = f*f;
    *qx = (KKZ*KKZ*f - rho)*f;
    *rx = -rho*f;
  }
  else if (sl_dim==2) {		// Cylindrical

    f = cylpot(*x);
    rho = cyldens(*x);

    *px = (*x)*f*f;
    *qx = (M2*f/(*x) + K2*f*(*x) - cylpotsl(*x)*(*x))*f;
    *rx = -rho*(*x)*f;
  }
  else {			// Spherical

    f = sphpot(*x);
    rho = sphdens(*x);

    *px = (*x)*(*x)*f*f;
    *qx = (L2*f - rho*(*x)*(*x))*f;
    *rx = -rho*(*x)*(*x)*f;

    if (*px<=0) 
      cerr << "Process " << myid << ": "
	   << "px<=0: x=" << *x << " f=" << f << "\n";
    if (*rx<=0)
      cerr << "Process " << myid << ": "
	   << "rx<=0: x=" << *x << " f=" << f << " rho=" << rho <<  "\n";

  }

  return 0;
}

