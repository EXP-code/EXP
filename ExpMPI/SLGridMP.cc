#define DEBUG 1
#define USE_TABLE 1

#include <stdlib.h>
#include <iomanip.h>
#include <f2c.h>

#include <SLGridMP.h>

MPI_Status status;

int SLGrid::mpi = 0;		// initially off

extern "C" {
  int sledge_(logical* job, doublereal* cons, logical* endfin, 
	      integer* invec, doublereal* tol, logical* type, 
	      doublereal* ev, integer* numx, doublereal* xef, doublereal* ef, 
	      doublereal* pdef, doublereal* t, doublereal* rho, 
	      integer* iflag, doublereal* store);
  double iv(double, double);
  double kn(int, double);
}


double exppot(double r)
{
  double y = 0.5 * r;

  return -y*(iv(0, y)*kn(1, y) - iv(1, y)*kn(0, y));
}

double expdens(double r)
{
  return 4.0*M_PI * 2.0*exp(-r);
}


void SLGrid::bomb(String oops)
{
  cerr << "SLGrid: " << oops << endl; 
  exit(-1);
}

				// Constructors

SLGrid::SLGrid(int MMAX, int NMAX, int NUMR, int NUMK, 
	       double RMIN, double RMAX, double L)
{
  int m, k;

  mmax = MMAX;
  nmax = NMAX;
  numr = NUMR;
  numk = NUMK;

  rmin = RMIN;
  rmax = RMAX;
  l = L;

  kv.setsize(1, NUMK);
  dk = 0.5*M_PI/L;
  for (k=1; k<=NUMK; k++) kv[k] = dk*k;

  init_table();


#ifdef DEBUG
  if (mpi)
    cout.form("Process %d: MPI is on!\n", myid);
  else
    cout.form("Process %d: MPI is off!\n", myid);

  cout.form("Process %d: MMAX=%d  NMAX=%d  NUMR=%d  NUMK=%d  RMIN=%f  RMAX=%f  L=%f\n", mpi_myid, MMAX, NMAX, NUMR, NUMK, RMIN, RMAX, L);

#endif

  if (mpi) {

    table =  new Table* [mmax+1];
    for (m=0; m<=mmax; m++) table[m] = new Table [numk] - 1;

    mpi_setup();

    if (mpi_myid) {
      compute_table_slave();

      //
      // <Receive completed table from master>
      //

      for (m=0; m<=mmax; m++) {
	for (k=1; k<=numk; k++) {
	  MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, 0,
		   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
	  mpi_unpack_table();      
	}
      }
      
    }
    else {			// BEGIN Master

      int slave = 0;
      int request_id = 1;
      double K;
      m=0; k=1;

      while (m<=mmax) {

	if (slave<mpi_numprocs-1 && m<=mmax) {	// Send request to slave
	  slave++;
      
#ifdef DEBUG    
	  cout.form("Master sending orders to Slave %d: (m,k)=(%d, %d)\n", 
		    slave, m, k);
#endif
	  MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	  MPI_Send(&m, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	  MPI_Send(&k, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
      
#ifdef DEBUG    
	  cout.form("Master gave orders to Slave %d: (m,k)=(%d, %d)\n", 
		    slave, m, k);
#endif

				// Increment counters
	  k++;
	  if (k>numk) {
	    k=1;
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
	  cout.form("Master gave orders to Slave %d: (m,k)=(%d, %d)\n", 
		    retid, m, k);
#endif
				// Increment counters
	  k++;
	  if (k>numk) {
	    k=1;
	    m++;
	  }
	}
      }
      
      //
      // <Wait for all slaves to return>
      //
  
      while (slave) {
	
	
	MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, 
		 MPI_COMM_WORLD, &status);
    
	mpi_unpack_table();      

	slave--;
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
	for (k=1; k<=numk; k++) {
	  int position = mpi_pack_table(&table[m][k], m, k);
	  for (slave=1; slave < mpi_numprocs; slave++)
	    MPI_Send(mpi_buf, position, MPI_PACKED, slave, 11, MPI_COMM_WORLD);
	}
      }


    } // END Master

  }
  else {

    table =  new Table* [mmax+1];
    for (m=0; m<=mmax; m++) {
      table[m] = new Table [numk] - 1;
      for (k=1; k<=numk; k++) {
	cerr.form("Begin [%d, %d] . . .\n", m, k);
	compute_table(&(table[m][k]), m, k);
	cerr.form(". . . done\n");
      }
    }
  }
}

SLGrid::~SLGrid()
{
  for (int m=0; m<=mmax; m++) delete [] (table[m]+1);
  delete [] table;
}

				// Members

double SLGrid::r_to_xi(double r)
{
  if (r<0.0) bomb("radius < 0!");
  return (r-1.0)/(r+1.0);
}
    
double SLGrid::xi_to_r(double xi)
{
  if (xi<-1.0) bomb("xi < -1!");
  if (xi>=1.0) bomb("xi >= 1!");

  return (1.0+xi)/(1.0 - xi);
}

double SLGrid::d_xi_to_r(double xi)
{
  if (xi<-1.0) bomb("xi < -1!");
  if (xi>=1.0) bomb("xi >= 1!");

  return 0.5*(1.0-xi)*(1.0-xi);
}

double SLGrid::get_pot(double x, int m, int n, int k, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  int indx;

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

#ifdef USE_TABLE
  return (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
    sqrt(table[m][k].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
  return (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
    sqrt(table[m][k].ev[n]) * exppot(xi_to_r(x));
#endif
}


double SLGrid::get_dens(double x, int m, int n, int k, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  int indx;

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  
#ifdef USE_TABLE
  return (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1]) *
    sqrt(table[m][k].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
  return (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1]) *
    sqrt(table[m][k].ev[n]) * expdens(xi_to_r(x));
#endif

}

double SLGrid::get_force(double x, int m, int n, int k, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  int indx;

  if (x<=xmin+dxi)
    indx = 1;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double p = (x - xi[indx])/dxi;
  
				// Use three point formula

				// Point -1: indx-1
				// Point  0: indx
				// Point  1: indx+1

  return d_xi_to_r(x)/dxi * (
			     (p - 0.5)*table[m][k].ef[n][indx-1]*p0[indx-1]
			     -2.0*p*table[m][k].ef[n][indx]*p0[indx]
			     + (p + 0.5)*table[m][k].ef[n][indx+1]*p0[indx+1]
			     ) / sqrt(table[m][k].ev[n]);
}


void SLGrid::get_pot(Matrix& mat, double x, int m, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  mat.setsize(1, numk, 1, nmax);

  int indx;

				// XI grid is same for all k

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int k=1; k<=numk; k++) {
    for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
      mat[k][n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
	sqrt(table[m][k].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
      mat[k][n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
	sqrt(table[m][k].ev[n]) * exppot(xi_to_r(x));
#endif
    }
  }

}


void SLGrid::get_dens(Matrix& mat, double x, int m, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  mat.setsize(1, numk, 1, nmax);

  int indx;

				// XI grid is same for all k

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int k=1; k<=numk; k++) {
    for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
      mat[k][n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])*
	sqrt(table[m][k].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
      mat[k][n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])*
	sqrt(table[m][k].ev[n]) * expdens(xi_to_r(x));
#endif
    }
  }

}


void SLGrid::get_force(Matrix& mat, double x, int m, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  mat.setsize(1, numk, 1, nmax);

  int indx;

				// XI grid is same for all k

  if (x<=xmin+dxi)
    indx = 1;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int k=1; k<=numk; k++) {
    for (int n=1; n<=nmax; n++) {
      mat[k][n] = fac * (
			 (p - 0.5)*table[m][k].ef[n][indx-1]*p0[indx-1]
			 -2.0*p*table[m][k].ef[n][indx]*p0[indx]
			 + (p + 0.5)*table[m][k].ef[n][indx+1]*p0[indx+1]
			 ) / sqrt(table[m][k].ev[n]);
    }
  }

}


void SLGrid::get_pot(Vector& vec, double x, int m, int k, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  vec.setsize(1, nmax);

  int indx;

				// XI grid is same for all k

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
      sqrt(table[m][k].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
    vec[n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])/
      sqrt(table[m][k].ev[n]) * exppot(xi_to_r(x));
#endif
  }

}


void SLGrid::get_dens(Vector& vec, double x, int m, int k, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  vec.setsize(1, nmax);

  int indx;

				// XI grid is same for all k

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    vec[n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])*
      sqrt(table[m][k].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
    vec[n] = (x1*table[m][k].ef[n][indx] + x2*table[m][k].ef[n][indx+1])*
      sqrt(table[m][k].ev[n]) * expdens(xi_to_r(x));
#endif
  }

}


void SLGrid::get_force(Vector& vec, double x, int m, int k, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  vec.setsize(1, nmax);

  int indx;

				// XI grid is same for all k

  if (x<=xmin+dxi)
    indx = 1;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int n=1; n<=nmax; n++) {
    vec[n] = fac * (
		    (p - 0.5)*table[m][k].ef[n][indx-1]*p0[indx-1]
		    -2.0*p*table[m][k].ef[n][indx]*p0[indx]
		    + (p + 0.5)*table[m][k].ef[n][indx+1]*p0[indx+1]
		    ) / sqrt(table[m][k].ev[n]);
  }

}


void SLGrid::get_pot(Matrix* mat, double x, int mMin, int mMax, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  if (mmax < mMax) mMax = mmax;

  for (int m=mMin; m<=mMax; m++)
    mat[m].setsize(1, numk, 1, nmax);

  int indx;

				// XI grid is same for all k

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int m=mMin; m<=mMax; m++) {
    for (int k=1; k<=numk; k++) {
      for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
	mat[m][k][n] = (x1*table[m][k].ef[n][indx] + 
			x2*table[m][k].ef[n][indx+1])/
	  sqrt(table[m][k].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
	mat[m][k][n] = (x1*table[m][k].ef[n][indx] + 
			x2*table[m][k].ef[n][indx+1])/
	  sqrt(table[m][k].ev[n]) * exppot(xi_to_r(x));
#endif
      }
    }
  }

}


void SLGrid::get_dens(Matrix* mat, double x, int mMin, int mMax, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  if (mmax < mMax) mMax = mmax;

  for (int m=mMin; m<=mMax; m++)
    mat[m].setsize(1, numk, 1, nmax);

  int indx;

				// XI grid is same for all k

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int m=mMin; m<=mMax; m++) {
    for (int k=1; k<=numk; k++) {
      for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
	mat[m][k][n] = (x1*table[m][k].ef[n][indx] + 
			x2*table[m][k].ef[n][indx+1])*
	  sqrt(table[m][k].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
	mat[m][k][n] = (x1*table[m][k].ef[n][indx] + 
			x2*table[m][k].ef[n][indx+1])*
	  sqrt(table[m][k].ev[n]) * expdens(xi_to_r(x));
#endif
      }
    }
  }

}


void SLGrid::get_force(Matrix* mat, double x, int mMin, int mMax, int which)
{
  if (which)
    x = r_to_xi(x);
  else {
    if (x<-1.0) x=-1.0;
    if (x>=1.0) x=1.0-1.0e-3;
  }

  if (mmax < mMax) mMax = mmax;

  for (int m=mMin; m<=mMax; m++)
    mat[m].setsize(1, numk, 1, nmax);

  int indx;

				// XI grid is same for all k

  if (x<=xmin+dxi)
    indx = 1;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double p = (x - xi[indx])/dxi;
  double fac = d_xi_to_r(x)/dxi;

  for (int m=mMin; m<=mMax; m++) {
    for (int k=1; k<=numk; k++) {
      for (int n=1; n<=nmax; n++) {
	mat[m][k][n] = fac * (
			      (p - 0.5)*table[m][k].ef[n][indx-1]*p0[indx-1]
			      -2.0*p*table[m][k].ef[n][indx]*p0[indx]
			      + (p + 0.5)*table[m][k].ef[n][indx+1]*p0[indx+1]
			      ) / sqrt(table[m][k].ev[n]);
      }
    }
  }
  
}


static double M2, K2;

void SLGrid::compute_table(struct Table* table, int m, int k)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol[6] = {1.0e-4,1.0e-5,  1.0e-4,1.0e-5,  1.0e-4,1.0e-5};
  int i, j, VERBOSE=0;
  integer NUM, N, M;
  logical type[8];
  logical endfin[2] = {1, 1};
  
#ifdef DEBUG
  //  VERBOSE=1;
#endif

  cons[6] = rmin;
  cons[7] = rmax;
  M2 = m*m;
  K2 = kv[k]*kv[k];
  NUM = numr;
  N = nmax;
  M = m;

  integer iflag[nmax], invec[nmax+3];
  double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
    pdef[NUM*N];
  double f;

  f = exppot(cons[6]);
  cons[2] = -1.0/(cons[6]*f);
  cons[4] = M/cons[7];
  f = exppot(cons[7]);
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
}


void SLGrid::init_table(void)
{

  int i;

  xi.setsize(0, numr-1);
  r.setsize(0, numr-1);
  p0.setsize(0, numr-1);
  d0.setsize(0, numr-1);

  xmin = (rmin - 1.0)/(rmin + 1.0);
  xmax = (rmax - 1.0)/(rmax + 1.0);
  dxi = (xmax-xmin)/(numr-1);

  for (i=0; i<numr; i++) {
    xi[i] = xmin + dxi*i;
    r[i] = xi_to_r(xi[i]);
    p0[i] = exppot(r[i]);
    d0[i] = expdens(r[i]);
  }

}


void SLGrid::compute_table_slave(void)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol[6] = {1.0e-4,1.0e-5,  1.0e-4,1.0e-5,  1.0e-4,1.0e-5};
  int i, j, k, VERBOSE=0;
  integer NUM;
  logical type[8];
  logical endfin[2] = {1, 1};
  
  struct Table table;
  int M, N, K;

#ifdef DEBUG
  //  VERBOSE=1;
#endif

#ifdef DEBUG
    cout.form("Slave %d begins . . .\n", mpi_myid);
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
    cout.form("Slave %d: ordered to compute (m, k) = (%d, %d)\n", mpi_myid,
	      M, K);
#endif

    cons[6] = rmin;
    cons[7] = rmax;
    M2 = M*M;
    K2 = kv[K]*kv[K];
    NUM = numr;
    N = nmax;

    integer iflag[N], invec[N+3];
    double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
      pdef[NUM*N];
    double f;

    f = exppot(cons[6]);
    cons[2] = -1.0/(cons[6]*f);
    cons[4] = M/cons[7];
    f = exppot(cons[7]);
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

#ifdef DEBUG
    cout.form("Slave %d: xef[0]=%f  xef[%d]=%f,  Rmin=%f  Rmax=%f\n", 
	      mpi_myid, xef[0], NUM-1, xef[NUM-1], cons[6], cons[7]);
#endif

    //     
    //     Open file for output.
    //
    sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	    t, rho, iflag, store);
    //
    //     Print results:
    //
    
#ifdef DEBUG
    cout.form("Slave %d: computed (m, k) = (%d, %d)\n", mpi_myid,
	      M, K);

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
    cout.form("Slave %d: sent to master (m, k) = (%d, %d)\n", mpi_myid, M, K);
#endif

  }

}


void SLGrid::mpi_setup(void)
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


int SLGrid::mpi_pack_table(struct Table* table, int m, int k)
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


void SLGrid::mpi_unpack_table(void)
{
  int i, j, length, position = 0;
  int m, k;

  int retid = status.MPI_SOURCE;

  MPI_Get_count( &status, MPI_PACKED, &length);


  MPI_Unpack( mpi_buf, length, &position, &m, 1, MPI_INT,
	      MPI_COMM_WORLD);
  MPI_Unpack( mpi_buf, length, &position, &k, 1, MPI_INT,
	      MPI_COMM_WORLD);

#ifdef DEBUG    
  cout.form("Process %d unpacking table entry from Process %d: (m, k)=(%d, %d)\n", 
	    mpi_myid, retid, m, k);
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



extern "C" int coeff_(doublereal* x, doublereal* px, doublereal* qx, 
	doublereal* rx)
{
  double f,rho;

  f = exppot(*x);
  rho = expdens(*x);

  *px = (*x)*f*f;
  *qx = (M2*f/(*x) + K2*f*(*x) - rho*(*x))*f;
  *rx = -rho*(*x)*f;

  return 0;
}

