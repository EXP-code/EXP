// #define DEBUG 1
// #define DEBUG_NAN 1
// #define GHQL 1
// #define STURM 1
// #define EIGEN 1

#include <iostream.h>
#include <iomanip.h>
#include <String.h>
#include <strstream.h>
#include <minmax.h>

#include <numerical.h>
#include <gaussQ.h>
#include "EmpOrth6.h"

Vector Symmetric_Eigenvalues_MSRCH(Matrix& a, Matrix& ef, int M);


#undef TINY
#define TINY 1.0e-16


Logic EmpCylSL::DENS=FALSE;
int EmpCylSL::NZOF=32;
int EmpCylSL::NINTR=32;
int EmpCylSL::NINTZ=256;
int EmpCylSL::NUMX=64;
int EmpCylSL::NUMY=128;
int EmpCylSL::NOUT=2;
double EmpCylSL::RMAX=20;
double EmpCylSL::HSCALE=0.05;
String EmpCylSL::CACHEFILE = ".eof.cache.file";
String EmpCylSL::TABLEFILE = ".eof.table.file";



EmpCylSL::EmpCylSL(void)
{
  NORDER=0;
  MPIset = FALSE;
  MPIset_eof = FALSE;
  coefs_made = FALSE;
  eof_made = FALSE;
  eof_recompute = TRUE;

  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  lwR = 0;
  lwZ = 0;

  SC = 0;
  SS = 0;

  accum_cos = 0;
  accum_sin = 0;

  accum_cos0 = 0;
  accum_sin0 = 0;

  mpi_double_buf1 = 0;
  mpi_double_buf2 = 0;
  mpi_double_buf3 = 0;
}

EmpCylSL::~EmpCylSL(void)
{
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

    delete [] mpi_double_buf1;
    delete [] mpi_double_buf2;
    delete [] mpi_double_buf3;

    for (int mm=0; mm<=MMAX; mm++) {
      
      for (int ir=1; ir<=NUMR; ir++) {

	for (int j1=1; j1<=rank2; j1++) {

	  delete [] (SC[mm][ir][j1] + 1);
	  if (mm) delete [] (SS[mm][ir][j1] + 1);
	}

	delete [] (SC[mm][ir] + 1);
	if (mm) delete [] (SS[mm][ir] + 1);
      }

      delete [] (SC[mm] + 1);
      if (mm) delete [] (SS[mm] + 1);
    }

    delete [] SC;
    delete [] SS;

  }

  if (lwR) delete lwR;
  if (lwZ) delete lwZ;
  if (accum_cos) delete [] accum_cos;
  if (accum_sin) delete [] accum_sin;

  if (accum_cos0) delete [] accum_cos0;
  if (accum_sin0) delete [] accum_sin0;

  if (MPIset) {
    delete [] MPIin;
    delete [] MPIout;
  }

  if (MPIset_eof) {
    delete [] MPIin_eof;
    delete [] MPIout_eof;
  }

}


EmpCylSL::EmpCylSL(CylindricalSL *sl, double zmax, int nord)
{
  NORDER = nord;
  ZMAX = zmax;
  ortho = sl;
  MMAX = sl->get_maxM();
  NUMR = sl->get_maxNR();
  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  lwR = 0;
  lwZ = 0;

  SC = 0;
  SS = 0;

  MPIset = FALSE;
  MPIset_eof = FALSE;
  coefs_made = FALSE;
  eof_made = FALSE;
  eof_recompute = TRUE;

  accum_cos = 0;
  accum_sin = 0;

  accum_cos0 = 0;
  accum_sin0 = 0;

  mpi_double_buf1 = 0;
  mpi_double_buf2 = 0;
  mpi_double_buf3 = 0;
}


void EmpCylSL::reset(CylindricalSL *sl, double zmax, int nord)
{
  NORDER = nord;
  ZMAX = zmax;
  ortho = sl;
  MMAX = sl->get_maxM();
  NUMR = sl->get_maxNR();

  lwR = 0;
  lwZ = 0;

  SC = 0;
  SS = 0;

  MPIset = FALSE;
  MPIset_eof = FALSE;
  coefs_made = FALSE;
  eof_made = FALSE;
  eof_recompute = TRUE;

  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  accum_cos = 0;
  accum_sin = 0;

  accum_cos0 = 0;
  accum_sin0 = 0;

  mpi_double_buf1 = 0;
  mpi_double_buf2 = 0;
  mpi_double_buf3 = 0;
}


void EmpCylSL::send_eof_grid()
{
  double *MPIbuf  = new double [MPIbufsz];

#ifdef DEBUG
  cerr.form("Process %d: size=%d\n", myid, MPIbufsz);
#endif
				// Send to slaves

  if (myid==0) {

    for (int m=0; m<=MMAX; m++) {

				// Coefficients for current step

      for (int v=0; v<rank3; v++) mpi_double_buf3[v] = accum_cos[m][v];
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

      for (int v=0; v<rank3; v++) mpi_double_buf3[v] = accum_sin[m][v];
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
      for (int v=0; v<rank3; v++) accum_cos[m][v] = mpi_double_buf3[v];


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
      for (int v=0; v<rank3; v++) accum_sin[m][v] = mpi_double_buf3[v];


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


void EmpCylSL::receive_eof(int request_id, int MM, int IR)
{

  int type, icnt, off;
  int mm, ir;

  MPI_Request Req[MPItable*NORDER+1], SReq[3];
  MPI_Status Status[MPItable*NORDER+1], Stat[3];

#ifdef DEBUG
  cerr.form("master listening . . . \n");
#endif

  MPI_Recv(&type, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	   MPI_COMM_WORLD, &status);

  int current_source = status.MPI_SOURCE;

#ifdef DEBUG
  cerr.form("master beginning to receive from %d . . . \n", current_source);
#endif

  MPI_Recv(&mm, 1, MPI_INT, current_source, MPI_ANY_TAG, 
	   MPI_COMM_WORLD, &status);
  MPI_Recv(&ir, 1, MPI_INT, current_source, MPI_ANY_TAG, 
	   MPI_COMM_WORLD, &status);

#ifdef DEBUG	
  cerr.form("master receiving from %d: type=%d   M=%d  ir=%d\n", current_source, type, mm, ir);
#endif

				// Send slave new orders
  if (request_id >=0) {
    MPI_Isend(&request_id, 1, MPI_INT, current_source, 1, MPI_COMM_WORLD, 
	      &SReq[0]);
    MPI_Isend(&MM, 1, MPI_INT, current_source, 2, MPI_COMM_WORLD, &SReq[1]);
    MPI_Isend(&IR, 1, MPI_INT, current_source, 3, MPI_COMM_WORLD, &SReq[2]);
  }
  else {
    MPI_Isend(&request_id, 1, MPI_INT, current_source, 1, MPI_COMM_WORLD,
	      &SReq[0]);
  }
  
				// Receive rest of data

  MPI_Irecv(mpi_double_buf1, NORDER, MPI_DOUBLE, current_source,
	    13, MPI_COMM_WORLD, &Req[0]);

  for (int n=0; n<NORDER; n++) {
    MPI_Irecv(&mpi_double_buf2[MPIbufsz*(MPItable*n+0)], 
	      MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+1, 
	      MPI_COMM_WORLD, &Req[MPItable*n+1]);
    MPI_Irecv(&mpi_double_buf2[MPIbufsz*(MPItable*n+1)], 
	      MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+2, 
	      MPI_COMM_WORLD, &Req[MPItable*n+2]);
    MPI_Irecv(&mpi_double_buf2[MPIbufsz*(MPItable*n+2)], 
	      MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+3, 
	      MPI_COMM_WORLD, &Req[MPItable*n+3]);
    if (DENS)
      MPI_Irecv(&mpi_double_buf2[MPIbufsz*(MPItable*n+3)], 
		MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+4, 
		MPI_COMM_WORLD, &Req[MPItable*n+4]);
  }
  
  MPI_Waitall(MPItable*NORDER+1, Req, Status);


				// Read from buffers

  for (int n=0; n<NORDER; n++) {
    if (type)
      accum_cos[mm][(ir-1)*NORDER+n] = mpi_double_buf1[n];
    else
      accum_sin[mm][(ir-1)*NORDER+n] = mpi_double_buf1[n];

  }

  for (int n=0; n<NORDER; n++) {
  
    off = MPIbufsz*(MPItable*n+0);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	potC[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[off+icnt++];
      else
	potS[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[off+icnt++];
	


    off = MPIbufsz*(MPItable*n+1);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	rforceC[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[off+icnt++];
      else
	rforceS[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[off+icnt++];
	

    off = MPIbufsz*(MPItable*n+2);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	zforceC[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[off+icnt++];
      else
	zforceS[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[off+icnt++];
    

    if (DENS) {
      off = MPIbufsz*(MPItable*n+3);
      icnt = 0;
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++)
	  if (type)
	    densC[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[off+icnt++];
	  else
	    densS[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[off+icnt++];
    }
  }
  
  if (request_id >=0) {
    MPI_Waitall(3, SReq, Stat);
#ifdef DEBUG
    cerr.form("master in make_eof, continuing orders: done waiting on Slave %d\n", current_source);
#endif
  }
  else {
    MPI_Waitall(1, SReq, Stat);
#ifdef DEBUG
    cerr.form("master in make_eof, final orders: done waiting on Slave %d\n", current_source);
#endif
  }

#ifdef DEBUG
  cerr.form("master finished receiving: type=%d   M=%d  ir=%d\n", type, mm, ir);
#endif

  return;

}

void EmpCylSL::compute_eof_grid(int request_id, int m, int ir)
{

#ifdef EIGEN
  static Vector sum(ev.getlow(), ev.gethigh());
  {
    cerr.form("Node %d: m=%d  ir=%d\n", myid, m, ir);
    sum[ev.getlow()] = ev[ev.getlow()];
    for (int nk=ev.getlow()+1; nk<=ev.gethigh(); nk++) sum[nk] = sum[nk-1] +
							 ev[nk];
    for (int nk=ev.getlow(); nk<=ev.gethigh(); nk++) {
      cerr.form("       %3d,%3d,%3d>   %13.6e  %13.6e\n", m, ir, nk, 
		ev[nk], sum[nk]);
      if (sum[nk]/sum[ev.gethigh()] > 1.0 - 1.0e-3) break;
    }
  }
#endif

  static mpi_in_progress = 0;
  static MPI_Request* Req = 0;
  static MPI_Status* Status;

  if (Req == 0) {
    Req = new MPI_Request [NORDER*MPItable+1];
    Status = new MPI_Status [NORDER*MPItable+1];
  }

  //  Read in coefficient matrix or
  //  make grid if needed


  double dcos, dsin;
  double lcos, lsin, cosn, sinn;

				// Sin/cos normalization
  double x, y, r, z, k;


  int icnt, off;

  for (int v=0; v<NORDER; v++) {
    tpot[v].zero();
    trforce[v].zero();
    tzforce[v].zero();
    if (DENS) tdens[v].zero();
  }

  for (int ix=0; ix<=NUMX; ix++) {

    x = -1.0 + dX*ix;
    r = ortho->SL()->xi_to_r(x);

    ortho->SL()->get_pot(tabp, r, m);
    ortho->SL()->get_force(tabf, r, m);
    if (DENS) ortho->SL()->get_dens(tabd, r, m);
    
    for (int iy=0; iy<=NUMY; iy++) {

      y = dY*iy + YMIN;
      z = y_to_z(y);

      dcos = cos(dk*(z+ZMAX));
      dsin = sin(dk*(z+ZMAX));
      cosn=1.0;
      sinn=0.0;

      for (int nk=0; nk<NZOF; nk++) {

	k = dk*nk;

	for (int v=0; v<NORDER; v++) {

	  tpot[v][ix][iy] +=  ef[v+1][nk+1]*tabp[nk][ir] * cosn * zn[nk];
	  tzforce[v][ix][iy] += -ef[v+1][nk+1]*tabp[nk][ir] * 
	    sinn * k * zn[nk];
	  trforce[v][ix][iy] += ef[v+1][nk+1]*tabf[nk][ir] * cosn * zn[nk];
	  if (DENS) 
	    tdens[v][ix][iy] +=  ef[v+1][nk+1]*tabd[nk][ir] * cosn * zn[nk];

	  if (nk) {
	    
	    tpot[v][ix][iy] += ef[v+1][NZOF+nk]*tabp[nk][ir] * sinn * zn[nk];
	    tzforce[v][ix][iy] += ef[v+1][NZOF+nk]*tabp[nk][ir] * 
	      cosn * k * zn[nk];
	    trforce[v][ix][iy] += ef[v+1][NZOF+nk]*tabf[nk][ir] * sinn * zn[nk];
	    if (DENS)
	      tdens[v][ix][iy] += ef[v+1][NZOF+nk]*tabd[nk][ir] * sinn * zn[nk];
	  }
	    
	}

	lcos = cosn;
	lsin = sinn;
	cosn = lcos*dcos - lsin*dsin;
	sinn = lsin*dcos + lcos*dsin;
	
      }
    }
  }

				// Wait for previous delivery?

  if (mpi_in_progress)
    MPI_Waitall(NORDER*MPItable+1, Req, Status);
  mpi_in_progress = 0;

				// Send stuff back to master
      
  MPI_Send(&request_id, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
  MPI_Send(&m, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
  MPI_Send(&ir, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
      
				// First send back part of accum

  for (int n=0; n<NORDER; n++) {
    mpi_double_buf1[n] = 0.0;
    for (int i=1; i<=rank2; i++) {
      if (request_id)
	mpi_double_buf1[n] += ef[n+1][i]*accum_cos0[m][ir][i];
      else
	mpi_double_buf1[n] += ef[n+1][i]*accum_sin0[m][ir][i];
    }
  }

  MPI_Isend(mpi_double_buf1, NORDER, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, &Req[0]);

    
  for (int n=0; n<NORDER; n++) {

				// normalization factors

    tpot[n] *= facZ;
    trforce[n] *= -facZ;
    tzforce[n] *= -facZ;
    if (DENS)
      tdens[n] *= 0.25/M_PI * facZ;

				// Potential

    off = MPIbufsz*(MPItable*n+0);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = tpot[n][ix][iy];
      
    MPI_Isend(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	      13 + MPItable*n+1, MPI_COMM_WORLD, &Req[MPItable*n+1]);

				// R force

    off = MPIbufsz*(MPItable*n+1);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = trforce[n][ix][iy];
    
    MPI_Isend(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	      13 + MPItable*n+2, MPI_COMM_WORLD, &Req[MPItable*n+2]);

				// Z force

    off = MPIbufsz*(MPItable*n+2);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = tzforce[n][ix][iy];
    
    MPI_Isend(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	      13 + MPItable*n+3, MPI_COMM_WORLD, &Req[MPItable*n+3]);

				// Density

    if (DENS) {
      off = MPIbufsz*(MPItable*n+3);
      icnt = 0;
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++)
	  mpi_double_buf2[off + icnt++] = tdens[n][ix][iy];
    
      MPI_Isend(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
		13 + MPItable*n+4, MPI_COMM_WORLD, &Req[MPItable*n+4]);

    }

  }

  mpi_in_progress = 1;

}


void EmpCylSL::setup_accumulation(void)
{
  if (eof_recompute) setup_eof();

  if (!accum_cos) {
    accum_cos = new Vector [MMAX+1];
    accum_sin = new Vector [MMAX+1];
#ifdef DEBUG_NAN
    cerr.form("Slave %d: tables allocated, MMAX=%d\n", myid, MMAX);
#endif // DEBUG_NAN
  }

  for (int m=0; m<=MMAX; m++) {
    accum_cos[m].setsize(0, rank3-1);
    accum_cos[m].zero();
    if (m>0) {
      accum_sin[m].setsize(0, rank3-1);
      accum_sin[m].zero();
    }
  }

  coefs_made = FALSE;

}

void EmpCylSL::setup_eof()
{
  if (!SC) {

    int rank0 = ortho->get_maxNZ();

    if (NZOF<1) NZOF = rank0;
  
    rank2 = 2*NZOF - 1;
    rank3 = NUMR*NORDER;

    zn.setsize(0, rank0-1);
    zn.zero();
    zn += 1.0;
    zn[0] = 1.0/sqrt(2.0);

    ev.setsize(1, rank2);
    ef.setsize(1, rank2, 1, rank2);

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

    dX = (1.0+ortho->SL()->r_to_xi(RMAX))/NUMX;
    YMIN = z_to_y(-CylindricalSL::ZMAX);
    YMAX = z_to_y( CylindricalSL::ZMAX);
    dY = (YMAX - YMIN)/NUMY;
    dk = M_PI/CylindricalSL::ZMAX;
    facZ = 1.0/sqrt(ZMAX);

    SC = new double*** [MMAX+1];
    SS = new double*** [MMAX+1];

    accum_cos0 = new Vector* [MMAX+1];
    accum_sin0 = new Vector* [MMAX+1];

    vc.setsize(0, MMAX, 0, rank3-1);
    vs.setsize(0, MMAX, 0, rank3-1);

    var.setsize(1, rank2, 1, rank2);

    for (int m=0; m<=MMAX; m++) {

      SC[m] = new double** [NUMR] - 1;
      if (m) SS[m] = new double** [NUMR] - 1;

      accum_cos0[m] = new Vector [NUMR] - 1;

      if (m) accum_sin0[m] = new Vector [NUMR] - 1;

      for (int i=1; i<=NUMR; i++) {

	SC[m][i] = new double* [rank2] - 1;
	if (m) SS[m][i] = new double* [rank2] - 1;

	for (int j1=1; j1<=rank2; j1++) {
	  SC[m][i][j1] = new double [rank2] - 1;
	  if (m) SS[m][i][j1] = new double [rank2] - 1;
	}

	accum_cos0[m][i].setsize(1, rank2);

	if (m) accum_sin0[m][i].setsize(1, rank2);

      }

    }
    
    table = new Matrix [MMAX+1];
    tablef = new Matrix [MMAX+1];

    //    facC.setsize(1, rank2, 1, NUMR);
    //    facS.setsize(1, rank2, 1, NUMR);

    facC.setsize(1, NUMR, 1, rank2);
    facS.setsize(1, NUMR, 1, rank2);

    MPIbufsz = (NUMX+1)*(NUMY+1);

    mpi_double_buf1 = new double [rank2*rank2];
    mpi_double_buf2 = new double [MPIbufsz*NORDER*MPItable];
    mpi_double_buf3 = new double [rank3];

  }


  for (int m=0; m<=MMAX; m++)  {
    for (int i=1; i<=NUMR; i++)  {
      accum_cos0[m][i].zero();
      if (m>0) accum_sin0[m][i].zero();
      for (int i1=1; i1<=rank2; i1++)  {
	for (int i2=1; i2<=rank2; i2++)  {
	  SC[m][i][i1][i2] = 0.0;
	  if (m>0) SS[m][i][i1][i2] = 0.0;
	}
      }
    }
  }

  eof_made = FALSE;
}


static const double fac0 = 1.0/sqrt(2.0*M_PI);
static const double facM = 1.0/sqrt(M_PI);


void EmpCylSL::accumulate_eof(double r, double z, double phi, double mass)
{
  if (eof_made) setup_eof();

  double msin, mcos;
  int nk, mm;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double dcos = cos(dk*(z+ZMAX));
  double dsin = sin(dk*(z+ZMAX));
  double lcos, lsin, cosn, sinn;
  double norm = -4.0*M_PI;
  double fac;

  Vector tmp;

  ortho->SL()->get_pot(table, r, 0, MMAX);

  for (mm=0; mm<=MMAX; mm++) {

    mcos = cos(phi*mm);
    msin = sin(phi*mm);

    if (mm)
      fac = facM*facZ;
    else
      fac = fac0*facZ;

    cosn=1.0;
    sinn=0.0;

    for (nk=0; nk<NZOF; nk++) {

      tmp = norm * table[mm][nk] * fac * zn[nk];

      for (int ir=1; ir<=NUMR; ir++) {
	facC[ir][nk+1] = tmp[ir] * mcos * cosn;
	if (nk) facC[ir][NZOF+nk] = tmp[ir] * mcos * sinn;
      
	if (mm) {
	  facS[ir][nk+1] = tmp[ir] * msin * cosn;
	  if (nk)  facS[ir][NZOF+nk] = tmp[ir] * msin * sinn;
	}
      }

      lcos = cosn;
      lsin = sinn;
      cosn = lcos*dcos - lsin*dsin;
      sinn = lsin*dcos + lcos*dsin;
    }

    
				// Mean and Covariance
    for (int ir=1; ir<=NUMR; ir++) {

      for (int i=1; i<=rank2; i++) {
	accum_cos0[mm][ir][i]  += facC[ir][i] * mass;
	if (mm) accum_sin0[mm][ir][i] += facS[ir][i] * mass;

	for (int j=i; j<=rank2; j++) {
	  SC[mm][ir][i][j] += facC[ir][i]*facC[ir][j] * mass;
	  if (mm) SS[mm][ir][i][j] += facS[ir][i]*facS[ir][j] * mass;
	}
      }

    }

  }

}


void EmpCylSL::make_eof(void)
{
  int icnt;
  double tmp;

  if (!MPIset_eof) {
    MPIin_eof  = new double [rank2*(rank2+1)/2];
    MPIout_eof = new double [rank2*(rank2+1)/2];
    MPIset_eof = TRUE;
  }
  
  //
  //  Distribute mean and covariance to all processes
  //

  for (int mm=0; mm<=MMAX; mm++) {

    for (int ir=1; ir<=NUMR; ir++) {

				// Mean
      icnt=0;
      for (int i=1; i<=rank2; i++) 
	MPIin_eof[icnt++] = accum_cos0[mm][ir][i];

      MPI_Allreduce ( MPIin_eof, MPIout_eof, rank2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      
      icnt=0;
      for (int i=1; i<=rank2; i++)
	accum_cos0[mm][ir][i] = MPIout_eof[icnt++];


				// Covariance
      icnt=0;
      for (int i=1; i<=rank2; i++)
	for (int j=i; j<=rank2; j++)
	  MPIin_eof[icnt++] = SC[mm][ir][i][j];
  
      MPI_Allreduce ( MPIin_eof, MPIout_eof, rank2*(rank2+1)/2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=1; i<=rank2; i++)
	for (int j=i; j<=rank2; j++)
	  SC[mm][ir][i][j] = MPIout_eof[icnt++];
  
    }
  }

  for (int mm=1; mm<=MMAX; mm++) {

    for (int ir=1; ir<=NUMR; ir++) {

				// Mean
      icnt=0;
      for (int i=1; i<=rank2; i++)
	MPIin_eof[icnt++] = accum_sin0[mm][ir][i];

      MPI_Allreduce ( MPIin_eof, MPIout_eof, rank2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=1; i<=rank2; i++)
	accum_sin0[mm][ir][i] = MPIout_eof[icnt++];


				// Covariance
      icnt=0;
      for (int i=1; i<=rank2; i++)
	for (int j=i; j<=rank2; j++)
	  MPIin_eof[icnt++] = SS[mm][ir][i][j];
  
      MPI_Allreduce ( MPIin_eof, MPIout_eof, rank2*(rank2+1)/2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=1; i<=rank2; i++)
	for (int j=i; j<=rank2; j++)
	  SS[mm][ir][i][j] = MPIout_eof[icnt++];
    }
  
  }


  if (myid==0) {

    MPI_Request req[3];
    MPI_Status  stat[3];
    int slave = 0;
    int request_id = 1;		// Begin with cosine case
    int M, ir;

    M = 0;			// Initial counters
    ir = 1;
    while (M<=MMAX) {
	
      // Send request to slave
      if (slave<numprocs-1) {
	  
	slave++;
	  
	MPI_Isend(&request_id, 1, MPI_INT, slave, 1, MPI_COMM_WORLD, &req[0]);
	MPI_Isend(&M,  1, MPI_INT, slave, 2, MPI_COMM_WORLD, &req[1]);
	MPI_Isend(&ir,  1, MPI_INT, slave, 3, MPI_COMM_WORLD, &req[2]);
	  
	// Increment counters
	ir++;
	if (ir>NUMR) {
	  request_id++;
	  if (request_id>1) {
	    M++;
	    request_id = 0;
	  }
	  ir = 1;
	}

	MPI_Waitall(3, req, stat);
#ifdef DEBUG
	cerr.form("master in make_eof: done waiting on Slave %d\n", slave);
#endif
      }
	
      if (slave == numprocs-1 && M<=MMAX) {
	  
	//
	// <Wait and receive and send new request>
	//
	  
	receive_eof(request_id, M, ir);
	  
	// Increment counters

	ir++;
	if (ir>NUMR) {
	  request_id++;
	  if (request_id>1) {
	    M++;
	    request_id = 0;
	  }
	  ir = 1;
	}
      }
    }
    
    //
    // <Wait for all slaves to return and flag to continue>
    //
      
    while (slave) {
      receive_eof(-1,0,0);
      slave--;
    }
      
  } else {

    int M, request_id, ir, icnt;
    MPI_Request SReq;
    MPI_Status  SStat;

    while (1) {

      MPI_Irecv(&request_id, 1, 
	       MPI_INT, 0, 1, MPI_COMM_WORLD, &SReq);
      MPI_Wait(&SReq, &SStat);

				// Done!
      if (request_id<0) break;

      MPI_Irecv(&M, 1, 
	       MPI_INT, 0, 2, MPI_COMM_WORLD, &SReq);
      MPI_Wait(&SReq, &SStat);
	
      MPI_Irecv(&ir, 1, 
	       MPI_INT, 0, 3, MPI_COMM_WORLD, &SReq);
      MPI_Wait(&SReq, &SStat);
	

#ifdef DEBUG
      cerr.form("Slave %d: received orders type=%d  M=%d  ir=%d\n",
		myid, request_id, M, ir);
#endif

      if (request_id) {

				// Complete symmetric part
    

	for (int i=1; i<rank2; i++) {
	  for (int j=i; j<=rank2; j++)
	    var[i][j] = SC[M][ir][i][j];
	}

	for (int i=1; i<rank2; i++) {
	  for (int j=i+1; j<=rank2; j++) {
	    var[j][i] = SC[M][ir][i][j];
	  }
	}
    
	double maxV = 0.0;
	for (int i=1; i<=rank2; i++) {
	  for (int j=i; j<=rank2; j++) {
	    tmp = fabs(var[i][j]);
	    if (tmp > maxV) maxV = tmp;
	  }
	}

	var /= maxV;
    
    /*==========================*/
    /* Solve eigenvalue problem */
    /*==========================*/
    
#if defined(STURM)
	ev = Symmetric_Eigenvalues_MSRCH(var, ef, NORDER);

#elif defined(GHQL)
	ev = var.Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#else
	ev = var.Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#endif
      }
      else {

				// Complete symmetric part
    

	for (int i=1; i<rank2; i++) {
	  for (int j=i; j<=rank2; j++)
	    var[i][j] = SS[M][ir][i][j];
	}

	for (int i=1; i<rank2; i++) {
	  for (int j=i+1; j<=rank2; j++) {
	    var[j][i] = SS[M][ir][i][j];
	  }
	}
    
	double maxV = 0.0;
	for (int i=1; i<=rank2; i++) {
	  for (int j=i; j<=rank2; j++) {
	    tmp = fabs(var[i][j]);
	    if (tmp > maxV) maxV = tmp;
	  }
	}
	
	var /= maxV;
    
    /*==========================*/
    /* Solve eigenvalue problem */
    /*==========================*/
    
#if defined(STURM)
	ev = Symmetric_Eigenvalues_MSRCH(var, ef, NORDER);

#elif defined(GHQL)
	ev = var.Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#else
	ev = var.Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#endif
      }

      compute_eof_grid(request_id, M, ir);

    }

  }
				// Send grid to all processes
  send_eof_grid();

  eof_made = TRUE;
  coefs_made = TRUE;
  eof_recompute = FALSE;
}



  

void EmpCylSL::accumulate(double r, double z, double phi, double mass)
{

  if (coefs_made) setup_accumulation();

  if (eof_recompute) {
    accumulate_eof(r, z, phi, mass);
    return;
  }

  double msin, mcos;
  int mm;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double norm = -4.0*M_PI;

  get_pot(vc, vs, r,  z);

  for (mm=0; mm<=MMAX; mm++) {

				// facZ is already included during table
				// generation
    if (mm) {
      mcos = cos(phi*mm) * facM;
      msin = sin(phi*mm) * facM;
    }
    else
      mcos = fac0;

    accum_cos[mm] += norm * mass * mcos * vc[mm];
    if (mm>0) accum_sin[mm] += norm * mass * msin * vs[mm];
    
  }

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
    MPIset = TRUE;
  }
  
#ifdef MPE_PROFILE
  MPE_Log_event(7, myid, "b_distrib_c");
#endif

  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      MPIin[mm*rank3 + nn] = accum_cos[mm][nn];
  
  MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (mm=0; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      accum_cos[mm][nn] = MPIout[mm*rank3 + nn];
  



  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      MPIin[mm*rank3 + nn] = accum_sin[mm][nn];
  
  MPI_Allreduce ( MPIin, MPIout, rank3*(MMAX+1),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

  for (mm=1; mm<=MMAX; mm++)
    for (nn=0; nn<rank3; nn++)
      accum_sin[mm][nn] = MPIout[mm*rank3 + nn];
  
#ifdef MPE_PROFILE
  MPE_Log_event(8, myid, "e_distrib_c");
#endif

  coefs_made = TRUE;
}




void EmpCylSL::accumulated_eval(double r, double z, double phi,
				double& p, double& fr, double& fz, 
				double &fp)
{
  if (!coefs_made) make_coefficients();

  fr = 0.0;
  fz = 0.0;
  fp = 0.0;
  p = 0.0;


  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double xx = ortho->SL()->r_to_xi(r);

  double X = (xx+1.0)/dX;
  double Y = (z_to_y(z) - YMIN)/dY;
  int ix = (int)X;
  int iy = (int)Y;

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

  double ccos, ssin=0.0, fac, facN;
  int n, mm;

  for (mm=0; mm<=MMAX; mm++) {

    if (mm) {
      ccos = cos(phi*mm) * facM;
      ssin = sin(phi*mm) * facM;
    }
    else
      ccos = fac0;

    for (n=0; n<rank3; n++) {

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

  }

}


double EmpCylSL::accumulated_dens_eval(double r, double z, double phi)
{
  if (!DENS) return 0.0;

  if (!coefs_made) make_coefficients();

  double ans = 0.0;

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double xx = ortho->SL()->r_to_xi(r);

  double X = (xx+1.0)/dX;
  double Y = (z_to_y(z) - YMIN)/dY;
  int ix = (int)X;
  int iy = (int)Y;

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
  int n, mm, v;

  for (mm=0; mm<=MMAX; mm++) {

    if (mm) {
      ccos = cos(phi*mm) * facM;
      ssin = sin(phi*mm) * facM;
    }
    else
      ccos = fac0;

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

  }

  return ans;
}


  
void EmpCylSL::get_pot(Matrix& Vc, Matrix& Vs, double r, double z)
{
  Vc.setsize(0, MMAX, 0, rank3-1);
  Vs.setsize(0, MMAX, 0, rank3-1);

  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double xx = ortho->SL()->r_to_xi(r);

  double X = (xx+1.0)/dX;
  double Y = (z_to_y(z) - YMIN)/dY;
  int ix = (int)X;
  int iy = (int)Y;

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

  int v;
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
    
  
void EmpCylSL::get_all(int mm, int ir, int nn, 
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


  if (z > ZMAX) z = ZMAX;
  if (z <-ZMAX) z = -ZMAX;

  double xx = ortho->SL()->r_to_xi(r);

  double X = (xx+1.0)/dX;
  double Y = (z_to_y(z) - YMIN)/dY;
  int ix = (int)X;
  int iy = (int)Y;

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

  double ccos, ssin, fac;

  ccos = cos(phi*mm);
  ssin = sin(phi*mm);

  int n = (ir-1)*NORDER + nn;

  p += ccos *
    (
     potC[mm][n][ix  ][iy  ] * c00 +
     potC[mm][n][ix+1][iy  ] * c10 +
     potC[mm][n][ix  ][iy+1] * c01 +
     potC[mm][n][ix+1][iy+1] * c11 
	 );

  fr += ccos *
    (
     rforceC[mm][n][ix  ][iy  ] * c00 +
     rforceC[mm][n][ix+1][iy  ] * c10 +
     rforceC[mm][n][ix  ][iy+1] * c01 +
     rforceC[mm][n][ix+1][iy+1] * c11
	 );
  
  fz += ccos *
    (
     zforceC[mm][n][ix  ][iy  ] * c00 +
     zforceC[mm][n][ix+1][iy  ] * c10 +
     zforceC[mm][n][ix  ][iy+1] * c01 +
     zforceC[mm][n][ix+1][iy+1] * c11 
     );
  
  fp += ssin * mm *
    (
     potC[mm][n][ix  ][iy  ] * c00 +
     potC[mm][n][ix+1][iy  ] * c10 +
     potC[mm][n][ix  ][iy+1] * c01 +
     potC[mm][n][ix+1][iy+1] * c11 
     );
  
  if (DENS)
  d += ccos *
    (
     densC[mm][n][ix  ][iy  ] * c00 +
     densC[mm][n][ix+1][iy  ] * c10 +
     densC[mm][n][ix  ][iy+1] * c01 +
     densC[mm][n][ix+1][iy+1] * c11 
     );


  if (mm) {
    
    p += ssin *
      (
       potS[mm][n][ix  ][iy  ] * c00 +
       potS[mm][n][ix+1][iy  ] * c10 +
       potS[mm][n][ix  ][iy+1] * c01 +
       potS[mm][n][ix+1][iy+1] * c11 
       );

    fr += ssin *
      (
       rforceS[mm][n][ix  ][iy  ] * c00 +
       rforceS[mm][n][ix+1][iy  ] * c10 +
       rforceS[mm][n][ix  ][iy+1] * c01 +
       rforceS[mm][n][ix+1][iy+1] * c11
       );

    fz += ssin *
      (
       zforceS[mm][n][ix  ][iy  ] * c00 +
       zforceS[mm][n][ix+1][iy  ] * c10 +
       zforceS[mm][n][ix  ][iy+1] * c01 +
	   zforceS[mm][n][ix+1][iy+1] * c11 
	   );

    fp += -ccos * mm *
	  (
	   potS[mm][n][ix  ][iy  ] * c00 +
	   potS[mm][n][ix+1][iy  ] * c10 +
	   potS[mm][n][ix  ][iy+1] * c01 +
	   potS[mm][n][ix+1][iy+1] * c11 
	   );
      
    if (DENS)
    d += ssin *
      (
       densS[mm][n][ix  ][iy  ] * c00 +
       densS[mm][n][ix+1][iy  ] * c10 +
       densS[mm][n][ix  ][iy+1] * c01 +
       densS[mm][n][ix+1][iy+1] * c11 
       );

  }
  
}


#include <stdio.h>

void EmpCylSL::dump_coefs(FILE *fout)
{
  double znorm;

  for (int mm=0; mm<=MMAX; mm++) {

    if (mm)
      znorm = sqrt(1.0/M_PI/ZMAX);
    else
      znorm = sqrt(0.5/M_PI/ZMAX);

    for (int nk=0; nk<NZOF; nk++) {

      fprintf(fout, "%4d %4d\n", mm, nk);

      for (int j=1; j<=NUMR; j++)
	fprintf(fout, " %13.4e", accum_cos0[mm][j][nk]*znorm);
      fprintf(fout, "\n");

      for (int j=1; j<=NUMR; j++)
	fprintf(fout, " %13.4e", accum_cos0[mm][j][nk+NZOF]*znorm);
      fprintf(fout, "\n");

      if (mm) {

	fprintf(fout, "%4d %4d\n", mm, nk);

	for (int j=1; j<=NUMR; j++)
	  fprintf(fout, " %13.4e", accum_sin0[mm][j][nk]*znorm);
	fprintf(fout, "\n");

	for (int j=1; j<=NUMR; j++)
	  fprintf(fout, " %13.4e", accum_sin0[mm][j][nk+NZOF]*znorm);
	fprintf(fout, "\n");
      }

    }
  }
}



void EmpCylSL::dump_coefs_binary(FILE *fout, double time)
{
  double p, znorm;

  coefheader.time = time;
  coefheader.mmax = MMAX;
  coefheader.nord = NZOF;
  coefheader.nmax = NUMR;

  fwrite(&coefheader, sizeof(CoefHeader), 1, fout);

  for (int mm=0; mm<=MMAX; mm++) {

    if (mm)
      znorm = sqrt(1.0/M_PI/ZMAX);
    else
      znorm = sqrt(0.5/M_PI/ZMAX);

    for (int nk=0; nk<NZOF; nk++) {
      
      for (int j=1; j<=NUMR; j++)
	fwrite(&(p=accum_cos0[mm][j][1+nk]*znorm), sizeof(double), 1, fout);

      if (nk) {
	for (int j=1; j<=NUMR; j++)
	  fwrite(&(p=accum_cos0[mm][j][nk+NZOF]*znorm), sizeof(double), 1, fout);
      } else {
	for (int j=1; j<=NUMR; j++)
	  fwrite(&(p=0.0), sizeof(double), 1, fout);
      }


      if (mm) {

	for (int j=1; j<=NUMR; j++)
	  fwrite(&(p=accum_sin0[mm][j][1+nk]*znorm), sizeof(double), 1, fout);

	if (nk) {
	  for (int j=1; j<=NUMR; j++)
	    fwrite(&(p=accum_sin0[mm][j][nk+NZOF]*znorm), sizeof(double), 1, fout);
	} else {
	  for (int j=1; j<=NUMR; j++)
	    fwrite(&(p=0.0), sizeof(double), 1, fout);
	}

      }
      
    }
  }
}



void EmpCylSL::dump_basis(char* name, int step)
{
  static String labels [] = {"pot.", "fr.", "fz.", "dens."};
  static int numx = 60;
  static int numy = 60;
  
  double rmax = 0.33*RMAX;
  double r, dr = rmax/(numx-1);
  double z, dz = 2.0*ZMAX/(numy-1);
  double xmin = ortho->SL()->r_to_xi(0);
  double xmax = ortho->SL()->r_to_xi(rmax);
  double ymin = z_to_y(-ZMAX);
  double ymax = z_to_y( ZMAX);


  float zz;

  char strbuf[128];
  double fac;
  int n, mm;
  ofstream** outC = new ofstream* [MPItable];
  ofstream** outS = new ofstream* [MPItable];
  
  for (mm=0; mm<=MMAX; mm++) {

    if (mm)
      fac = facM;
    else
      fac = fac0;

    for (n=0; n<=min(NOUT, rank3); n++) {

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
	outC[i]->write(&(zz=-ZMAX), sizeof(float));
	outC[i]->write(&(zz= ZMAX), sizeof(float));
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
	  outS[i]->write(&(zz=-ZMAX), sizeof(float));
	  outS[i]->write(&(zz= ZMAX), sizeof(float));
	}

      }


				// Ok, write data

      for (int k=0; k<numy; k++) {

	z = -ZMAX + dz*k;

	for (int j=0; j<numx; j++) {
	  
	  r = dr*j;

	  double X = (ortho->SL()->r_to_xi(r)+1.0)/dX;
	  double Y = (z_to_y(z) - YMIN)/dY;
	  int ix = (int)X;
	  int iy = (int)Y;

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



