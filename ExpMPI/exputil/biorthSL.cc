#pragma implementation

// #define DEBUG 1
// #define DEBUG_SLEDGE 1
#define SLEDGE_VERBOSE 1
#define USE_TABLE 1

#include <stdlib.h>
#include <iomanip>
#include <f2c.h>

#include <biorthSL.h>

MPI_Status status;


int SphereSL::mpi = 0;		// initially off
int SphereSL::cache = 1;	// initially yes
double SphereSL::rscale = 1;

extern "C" {
  int sledge_(logical* job, doublereal* cons, logical* endfin, 
	      integer* invec, doublereal* tol, logical* type, 
	      doublereal* ev, integer* numx, doublereal* xef, doublereal* ef, 
	      doublereal* pdef, doublereal* t, doublereal* rho, 
	      integer* iflag, doublereal* store);
  double iv(double, double);
  double kn(int, double);
}


static double L2, M2, K2;
static int sl_dim;

func_1d sphPot, sphDens, cylPot, cylPotSL, cylDens;

/*
void SphereSL::bomb(String oops)
{
  cerr << "SphereSL: " << oops << endl; 
  exit(-1);
}
*/
				// Constructors

SphereSL::SphereSL(int LMAX, int NMAX, int NUMR,
		   double RMIN, double RMAX,
		   func_1d pot, func_1d dpot, func_1d dens)
{
  BiorthID = "SphereSL";
  dof = 3;

  int l;

  lmax = LMAX;
  nmax = NMAX;
  numr = NUMR;

  rmin = RMIN;
  rmax = RMAX;

  sphpot = pot;
  sphdpot = dpot;
  sphdens = dens;

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
	MPI_Recv(mpi_buf, mpi_bufsz, MPI_PACKED, 0,
		 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
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
		 << ": l=" << l << endl;

#endif
	    MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	    MPI_Send(&l, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
      
#ifdef DEBUG    
	    cout << "Master gave orders to Slave " << slave 
		 << ": l=" << l << "\n";
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
	    cout << "Master gave orders to Slave " << retid 
		 << ": l=" << l << endl;
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
	for (slave=1; slave < mpi_numprocs; slave++)
	  MPI_Send(mpi_buf, position, MPI_PACKED, slave, 11, MPI_COMM_WORLD);
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
}

const string sph_cache_name = ".slgrid_sph_cache";

int SphereSL::read_cached_table(void)
{
  ifstream in(sph_cache_name.c_str());
  if (!in) return 0;

  int LMAX, NMAX, NUMR, i, j;
  double RMIN, RMAX;

  cerr << "SphereSL::read_cached_table: trying to read cached table . . . \n";

  in.read((char *)&LMAX, sizeof(int));		if(!in || LMAX!=lmax) return 0;
  in.read((char *)&NMAX, sizeof(int));		if(!in || NMAX!=nmax) return 0;
  in.read((char *)&NUMR, sizeof(int));		if(!in || NUMR!=numr) return 0;
  in.read((char *)&RMIN, sizeof(double));	if(!in || RMIN!=rmin) return 0;
  in.read((char *)&RMAX, sizeof(double));	if(!in || RMAX!=rmax) return 0;

  for (int l=0; l<=lmax; l++) {

    in.read((char *)&table[l].l, sizeof(int));

				// Double check
    if (table[l].l != l) {
      cerr << "SphereSL: error reading <" << sph_cache_name << ">\n";
      cerr << "SphereSL: l: read value (" << table[l].l << ") != internal value (" << l << ")\n";
	return 0;
    }

    table[l].ev.setsize(1, nmax);
    table[l].ef.setsize(1, nmax, 0, numr-1);

    for (j=1; j<=nmax; j++) in.read((char *)&table[l].ev[j], sizeof(double));

    for (j=1; j<=nmax; j++)
      for (i=0; i<numr; i++)
	in.read((char *)&table[l].ef[j][i], sizeof(double));
  }

  cerr << "SphereSL::read_cached_table: Success!!\n";
  return 1;
}


void SphereSL::write_cached_table(void)
{
  ofstream out(sph_cache_name.c_str());
  if (!out) {
    cerr << "SphereSL: error writing <" << sph_cache_name << ">\n";
    return;
  }

  int i, j;

  out.write((char *)&lmax, sizeof(int));
  out.write((char *)&nmax, sizeof(int));
  out.write((char *)&numr, sizeof(int));
  out.write((char *)&rmin, sizeof(double));
  out.write((char *)&rmax, sizeof(double));

  for (int l=0; l<=lmax; l++) {

    out.write((char *)&table[l].l, sizeof(int));

    for (j=1; j<=nmax; j++)
      out.write((char *)&table[l].ev[j], sizeof(double));

    for (j=1; j<=nmax; j++)
      for (i=0; i<numr; i++)
	out.write((char *)&table[l].ef[j][i], sizeof(double));

  }


  cerr << "SphereSL::write_cached_table: done!!\n";
  return ;
}


SphereSL::~SphereSL(void)
{
  delete [] table;
}

				// Members

double SphereSL::r_to_rb(const double r)
{
  if (r<0.0) bomb("radius < 0!");
  return (r/rscale-1.0)/(r/rscale+1.0);
}
    
double SphereSL::rb_to_r(const double rb)
{
  if (rb<-1.0) bomb("rb < -1!");
  if (rb>=1.0) bomb("rb >= 1!");

  return rscale*(1.0+rb)/(1.0 - rb);
}

double SphereSL::d_r_to_rb(const double r)
{
  if (r<0.0) bomb("radius < 0!");

  double fac = r/rscale + 1.0;;
  return 2.0*(fac*fac)/rscale;
}


double SphereSL::d_rb_to_r(const double rb)
{
  if (rb<-1.0) bomb("rb < -1!");
  if (rb>=1.0) bomb("rb >= 1!");

  return 0.5*(1.0-rb)*(1.0-rb)/rscale;
}

double SphereSL::potl(const int n, const int l, const double xx)
{
  double x = xx;
  if (x<-1.0) x=-1.0;
  if (x>=1.0) x=1.0-1.0e-3;

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
  return (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
    sqrt(table[l].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
  return (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
    sqrt(table[l].ev[n]) * (*sphpot)(xi_to_r(x));
#endif
}


double SphereSL::dens(const int n, const int l, const double xx)
{
  double x = xx;
  if (x<-1.0) x=-1.0;
  if (x>=1.0) x=1.0-1.0e-3;

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
  return -(x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1]) *
    sqrt(table[l].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
  return -(x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1]) *
    sqrt(table[l].ev[n]) * (*sphdens)(xi_to_r(x));
#endif

}

double SphereSL::get_potl(const double r, const int l, const Vector& coef)
{
  double x = r_to_rb(r);

  int indx;

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  double ans = 0.0;

  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    ans += (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
      sqrt(table[l].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]) * coef[n];
#else
    ans += (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
      sqrt(table[l].ev[n]) * (*sphpot)(xi_to_r(x)) * coef[n];
#endif
  }

  return ans;
}


double SphereSL::get_dens(const double r, const int l, const Vector& coef)
{
  double x = r_to_rb(r);

  int indx;

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  double ans = 0.0;

  for (int n=1; n<=nmax; n++) {
#ifdef USE_TABLE
    ans += (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])*
      sqrt(table[l].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]) * coef[n];
#else
    ans += (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])*
      sqrt(table[l].ev[n]) * (*sphdens)(xi_to_r(x)) * coef[n];
#endif
  }

  return -ans;
}

void SphereSL::potl(const int nn, const int l, const double x,
		    Vector& t)
{
  int indx;

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=1; n<=min(nn, nmax); n++) {
#ifdef USE_TABLE
    t[n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
      sqrt(table[l].ev[n]) * (x1*p0[indx] + x2*p0[indx+1]);
#else
    t[n] = (x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])/
      sqrt(table[l].ev[n]) * (*sphpot)(xi_to_r(x));
#endif
  }
  
}


void SphereSL::dens(const int nn, const int l, const double x,
		    Vector& t)
{
  int indx;

  if (x<=xmin)
    indx = 0;
  else if (x>=xmax) 
    indx = numr - 2;
  else 
    indx = (int)( (x-xmin)/dxi );

  double x1 = (xi[indx+1] - x)/dxi;
  double x2 = (x - xi[indx])/dxi;
  

  for (int n=1; n<=min(nn, nmax); n++) {
#ifdef USE_TABLE
    t[n] = -(x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])*
      sqrt(table[l].ev[n]) * (x1*d0[indx] + x2*d0[indx+1]);
#else
    t[n] = -(x1*table[l].ef[n][indx] + x2*table[l].ef[n][indx+1])*
      sqrt(table[l].ev[n]) * (*sphdens)(xi_to_r(x));
#endif
  }

}


void SphereSL::compute_table(struct TableSph* table, int l)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol[6] = {1.0e-4,1.0e-5,  1.0e-4,1.0e-5,  1.0e-4,1.0e-5};
  int i, j, k, VERBOSE=0;
  integer NUM, N;
  logical type[8] = {0, 0, 1, 0, 0, 0, 1, 0};
  logical endfin[2] = {1, 1};
  
#ifdef DEBUG_SLEDGE
  VERBOSE = SLEDGE_VERBOSE;
#endif

  cons[6] = rmin;
  cons[7] = rmax;
  L2 = l*(l+1);
  NUM = numr;
  N = nmax;
  
  integer iflag[nmax], invec[nmax+3];
  double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
    pdef[NUM*N];
  double f;

  if (l==0) {
    f = (*sphpot)(cons[6]);
    cons[0] = (*sphdpot)(cons[6])/f;
    cons[2] = 1.0/(cons[6]*cons[6]*f*f);
  }
  else
    cons[0] = 1.0;

  f = (*sphpot)(cons[7]);
  cons[5] = 1.0/(cons[7]*cons[7]*f);

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

  sphPot = sphpot;
  sphDens = sphdens;

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
}


void SphereSL::init_table(void)
{

  int i;

  xi.setsize(0, numr-1);
  r.setsize(0, numr-1);
  p0.setsize(0, numr-1);
  d0.setsize(0, numr-1);

  xmin = (rmin/rscale - 1.0)/(rmin/rscale + 1.0);
  xmax = (rmax/rscale - 1.0)/(rmax/rscale + 1.0);
  dxi = (xmax-xmin)/(numr-1);

  for (i=0; i<numr; i++) {
    xi[i] = xmin + dxi*i;
    r[i] = rb_to_r(xi[i]);
    p0[i] = (*sphpot)(r[i]);
    d0[i] = (*sphdens)(r[i]);
  }

}


void SphereSL::compute_table_slave(void)
{

  double cons[8] = {0.0, 0.0, 0.0, 0.0,   0.0, 0.0,   0.0, 0.0};
  double tol[6] = {1.0e-4,1.0e-5,  1.0e-4,1.0e-5,  1.0e-4,1.0e-5};
  int i, j, k, VERBOSE=0;
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
    cout << "Slave " << mpi_myid << ": ordered to compute l = " << L << endl;
#endif

    cons[6] = rmin;
    cons[7] = rmax;
    L2 = L*(L+1);
    NUM = numr;
    N = nmax;

    integer iflag[N], invec[N+3];
    double ev[N], *t, *rho, store[26*(NUM+16)], xef[NUM+16], ef[NUM*N],
      pdef[NUM*N];
    double f;

    if (L==0) {
      f = (*sphpot)(cons[6]);
      cons[0] = (*sphdpot)(cons[6])/f;
      cons[2] = 1.0/(cons[6]*cons[6]*f*f);
    }
    else
      cons[0] = 1.0;
    
    f = (*sphpot)(cons[7]);
    cons[5] = 1.0/(cons[7]*cons[7]*f);

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

    sphPot = sphpot;
    sphDens = sphdens;

    sledge_(job, cons, endfin, invec, tol, type, ev, &NUM, xef, ef, pdef,
	    t, rho, iflag, store);
    //
    //     Print results:
    //
    
#ifdef DEBUG
    cout << "Slave " << mpi_myid << ": computed l = " << L << endl;

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

    table.l = L;
  
    int position = mpi_pack_table(&table, L);
    MPI_Send(mpi_buf, position, MPI_PACKED, 0, 11, MPI_COMM_WORLD);

#ifdef DEBUG
    cout << "Slave " << mpi_myid << ": sent to master l = " << L << endl;
#endif

  }

}


void SphereSL::mpi_setup(void)
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


int SphereSL::mpi_pack_table(struct TableSph* table, int l)
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


void SphereSL::mpi_unpack_table(void)
{
  int i, j, length, position = 0;
  int l;

  int retid = status.MPI_SOURCE;

  MPI_Get_count( &status, MPI_PACKED, &length);


  MPI_Unpack( mpi_buf, length, &position, &l, 1, MPI_INT,
	      MPI_COMM_WORLD);

#ifdef DEBUG    
  cout << "Process " << mpi_myid << " unpacking table entry from Process " 
       << retid << ": l=" << l << endl;
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


extern "C" int coeff_(doublereal* x, doublereal* px, doublereal* qx, 
	doublereal* rx)
{
  double f,rho;

  if (sl_dim==2) {

    f = cylPot(*x);

    *px = (*x)*f*f;
    *qx = (M2*f/(*x) + K2*f*(*x) - cylPotSL(*x)*(*x))*f;
    *rx = -cylDens(*x)*(*x)*f;
  }
  else {

    f = (*sphPot)(*x);
    rho = (*sphDens)(*x);

    *px = (*x)*(*x)*f*f;
    *qx = (L2*f - rho*(*x)*(*x))*f;
    *rx = -rho*(*x)*(*x)*f;
  }

  return 0;
}

