// #define DEBUG 1
// #define DEBUG_NAN 1
// #define GHQL 1
// #define STURM 1

#include <iostream.h>
#include <iomanip.h>
#include <String.h>

#include <numerical.h>
#include <gaussQ.h>
#include "EmpOrth4.h"

Vector Symmetric_Eigenvalues_MSRCH(Matrix& a, Matrix& ef, int M);


#undef TINY
#define TINY 1.0e-16


int EmpCylSL::NZOF=32;
int EmpCylSL::NINTR=32;
int EmpCylSL::NINTZ=256;
int EmpCylSL::NUMX=64;
int EmpCylSL::NUMY=128;
double EmpCylSL::RMAX=20;
double EmpCylSL::HSCALE=0.05;
String EmpCylSL::CACHEFILE = ".eof.cache.file";
String EmpCylSL::TABLEFILE = ".eof.table.file";

extern double dens0(int, double, double);


EmpCylSL::EmpCylSL(void)
{
  NORDER=0;
  MPIset = FALSE;
  MPIset_eof = FALSE;
  coefs_made = FALSE;
  eof_made = FALSE;
  eof_recompute = TRUE;

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
}

EmpCylSL::~EmpCylSL(void)
{
  if (SC) {

    for (int m=0; m<=MMAX; m++) {
      delete [] potC[m];
      delete [] densC[m];
      delete [] rforceC[m];
      delete [] zforceC[m];
      if (m) {
	delete [] potS[m];
	delete [] densS[m];
	delete [] rforceS[m];
	delete [] zforceS[m];
      }
    }


    delete [] potC;
    delete [] densC;
    delete [] rforceC;
    delete [] zforceC;

    delete [] potS;
    delete [] densS;
    delete [] rforceS;
    delete [] zforceS;

    delete [] mpi_double_buf1;
    delete [] mpi_double_buf2;

    for (int mm=0; mm<=MMAX; mm++) {
      delete [] (SC[mm] + 1);
      if (mm)
	delete [] (SS[mm] + 1);
    }
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

  accum_cos = 0;
  accum_sin = 0;

  accum_cos0 = 0;
  accum_sin0 = 0;

  mpi_double_buf1 = 0;
  mpi_double_buf2 = 0;
}


void EmpCylSL::send_eof_grid()
{
  int MPIbufsz = (NUMX+1)*(NUMY+1);
  double *MPIbuf  = new double [MPIbufsz];

#ifdef DEBUG
  cerr.form("Process %d: size=%d\n", myid, MPIbufsz);
#endif
				// Send to slaves

  if (myid==0) {

    for (int m=0; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = potC[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = densC[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = rforceC[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = zforceC[m][v][ix][iy];
	
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
      }

    }

    for (int m=1; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = potS[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = densS[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = rforceS[m][v][ix][iy];

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    MPIbuf[ix*(NUMY+1) + iy] = zforceS[m][v][ix][iy];
	
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
      }

    }

  }
  else {
				// Get tables from Master

    for (int m=0; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    potC[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    densC[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    rforceC[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    zforceC[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
      }
    }

    for (int m=1; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    potS[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    densS[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    rforceS[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
	MPI_Bcast(MPIbuf, MPIbufsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    zforceS[m][v][ix][iy] = MPIbuf[ix*(NUMY+1) + iy];
  
      }
    }

  }

  delete [] MPIbuf;
  
}


int EmpCylSL::receive_eof()
{

  int type, icnt;
  int mm, ir;

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

  MPI_Recv(mpi_double_buf1, NORDER, MPI_DOUBLE, current_source,
	   MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  for (int n=0; n<NORDER; n++) {
    if (type)
      accum_cos[mm][(ir-1)*NORDER+n] = mpi_double_buf1[n];
    else
      accum_sin[mm][(ir-1)*NORDER+n] = mpi_double_buf1[n];

  }

  for (int n=0; n<NORDER; n++) {
  
    MPI_Recv(mpi_double_buf2, (NUMX+1)*(NUMY+1), MPI_DOUBLE, current_source,
	     MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    icnt=0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	potC[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[icnt++];
      else
	potS[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[icnt++];
	

    MPI_Recv(mpi_double_buf2, (NUMX+1)*(NUMY+1), MPI_DOUBLE, current_source,
	     MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    icnt=0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	densC[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[icnt++];
      else
	densS[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[icnt++];
	

    MPI_Recv(mpi_double_buf2, (NUMX+1)*(NUMY+1), MPI_DOUBLE, current_source,
	     MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    icnt=0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	rforceC[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[icnt++];
      else
	rforceS[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[icnt++];
	

    MPI_Recv(mpi_double_buf2, (NUMX+1)*(NUMY+1), MPI_DOUBLE, current_source,
	     MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    icnt=0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	zforceC[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[icnt++];
      else
	zforceS[mm][(ir-1)*NORDER+n][ix][iy]  = mpi_double_buf2[icnt++];
	
    
  }

#ifdef DEBUG
  cerr.form("master finished receiving: type=%d   M=%d  ir=%d\n", type, mm, ir);
#endif

  return current_source;

}

void EmpCylSL::compute_eof_grid(int request_id, int m, int ir)
{

  //  Read in coefficient matrix or
  //  make grid if needed


  double dcos, dsin;
  double lcos, lsin, cosn, sinn;

				// Sin/cos normalization
  Matrix tabp, tabd, tabf;

  double x, y, r, z, k;


  int icnt;

  for (int v=0; v<NORDER; v++) {
    tpot[v].zero();
    tdens[v].zero();
    trforce[v].zero();
    tzforce[v].zero();
  }

  for (int ix=0; ix<=NUMX; ix++) {

    x = -1.0 + dX*ix;
    r = ortho->SL()->xi_to_r(x);

    ortho->SL()->get_pot(tabp, r, m);
    ortho->SL()->get_dens(tabd, r, m);
    ortho->SL()->get_force(tabf, r, m);
    
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
	  tdens[v][ix][iy] +=  ef[v+1][nk+1]*tabd[nk][ir] * cosn * zn[nk];
	  trforce[v][ix][iy] += ef[v+1][nk+1]*tabf[nk][ir] * cosn * zn[nk];

	  if (nk) {
	    
	    tpot[v][ix][iy] += ef[v+1][NZOF+nk]*tabp[nk][ir] * sinn * zn[nk];
	    tzforce[v][ix][iy] += ef[v+1][NZOF+nk]*tabp[nk][ir] * 
	      cosn * k * zn[nk];
	    tdens[v][ix][iy] += ef[v+1][NZOF+nk]*tabd[nk][ir] * sinn * zn[nk];
	    trforce[v][ix][iy] += ef[v+1][NZOF+nk]*tabf[nk][ir] * sinn * zn[nk];
	  }
	    
	}

	lcos = cosn;
	lsin = sinn;
	cosn = lcos*dcos - lsin*dsin;
	sinn = lsin*dcos + lcos*dsin;
	
      }
    }
  }


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

  MPI_Send(mpi_double_buf1, NORDER, MPI_DOUBLE, 0, 12, MPI_COMM_WORLD);

    
  for (int n=0; n<NORDER; n++) {

				// normalization factors

    tpot[n] *= 1.0/sqrt(ZMAX);
    tdens[n] *= 0.25/M_PI/sqrt(ZMAX);
    trforce[n] *= -1.0/sqrt(ZMAX);
    tzforce[n] *= -1.0/sqrt(ZMAX);

				// Potential

    icnt=0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[icnt++] = tpot[n][ix][iy];
      
    MPI_Send(mpi_double_buf2, (NUMX+1)*(NUMY+1), MPI_DOUBLE, 0, 12, 
	     MPI_COMM_WORLD);

				// Density

    icnt=0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[icnt++] = tdens[n][ix][iy];
    
    MPI_Send(mpi_double_buf2, (NUMX+1)*(NUMY+1), MPI_DOUBLE, 0, 12, 
	     MPI_COMM_WORLD);

				// R force

    icnt=0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[icnt++] = trforce[n][ix][iy];
    
    MPI_Send(mpi_double_buf2, (NUMX+1)*(NUMY+1), MPI_DOUBLE, 0, 12, 
	     MPI_COMM_WORLD);

				// Z force

    icnt=0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[icnt++] = tzforce[n][ix][iy];
    
    MPI_Send(mpi_double_buf2, (NUMX+1)*(NUMY+1), MPI_DOUBLE, 0, 12, 
	       MPI_COMM_WORLD);

  }

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

    potC = new Matrix* [MMAX];
    densC = new Matrix* [MMAX];
    rforceC = new Matrix* [MMAX];
    zforceC = new Matrix* [MMAX];

    potS = new Matrix* [MMAX];
    densS = new Matrix* [MMAX];
    rforceS = new Matrix* [MMAX];
    zforceS = new Matrix* [MMAX];

    for (int m=0; m<=MMAX; m++) {

      potC[m] = new Matrix [rank3];
      densC[m] = new Matrix [rank3];
      rforceC[m] = new Matrix [rank3];
      zforceC[m] = new Matrix [rank3];

      for (int v=0; v<rank3; v++) {
	potC[m][v].setsize(0, NUMX, 0, NUMY);
	densC[m][v].setsize(0, NUMX, 0, NUMY);
	rforceC[m][v].setsize(0, NUMX, 0, NUMY);
	zforceC[m][v].setsize(0, NUMX, 0, NUMY);
      }


    }

    for (int m=1; m<=MMAX; m++) {

      potS[m] = new Matrix [rank3];
      densS[m] = new Matrix [rank3];
      rforceS[m] = new Matrix [rank3];
      zforceS[m] = new Matrix [rank3];

      for (int v=0; v<rank3; v++) {
	potS[m][v].setsize(0, NUMX, 0, NUMY);
	densS[m][v].setsize(0, NUMX, 0, NUMY);
	rforceS[m][v].setsize(0, NUMX, 0, NUMY);
	zforceS[m][v].setsize(0, NUMX, 0, NUMY);
      }

    }

    tpot = new Matrix [NORDER];
    tdens = new Matrix [NORDER];
    trforce = new Matrix [NORDER];
    tzforce = new Matrix [NORDER];

    for (int n=0; n<NORDER; n++) {
      tpot[n].setsize(0, NUMX, 0, NUMY);
      tdens[n].setsize(0, NUMX, 0, NUMY);
      trforce[n].setsize(0, NUMX, 0, NUMY);
      tzforce[n].setsize(0, NUMX, 0, NUMY);
    }

    dX = (1.0+ortho->SL()->r_to_xi(RMAX))/NUMX;
    YMIN = z_to_y(-CylindricalSL::ZMAX);
    YMAX = z_to_y( CylindricalSL::ZMAX);
    dY = (YMAX - YMIN)/NUMY;
    dk = M_PI/CylindricalSL::ZMAX;
    facZ = 1.0/sqrt(ZMAX);

    SC = new Matrix* [MMAX+1];
    SS = new Matrix* [MMAX+1];

    accum_cos0 = new Vector* [MMAX+1];
    accum_sin0 = new Vector* [MMAX+1];

    vc.setsize(0, MMAX, 0, rank3-1);
    vs.setsize(0, MMAX, 0, rank3-1);

    for (int m=0; m<=MMAX; m++) {

      SC[m] = new Matrix [NUMR] - 1;
      if (m) SS[m] = new Matrix [NUMR] - 1;

      accum_cos0[m] = new Vector [NUMR] - 1;

      if (m) accum_sin0[m] = new Vector [NUMR] - 1;

      for (int i=1; i<=NUMR; i++) {

	SC[m][i].setsize(1, rank2, 1, rank2);
	if (m) SS[m][i].setsize(1, rank2, 1, rank2);

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

    mpi_double_buf1 = new double [rank2*rank2];
    mpi_double_buf2 = new double [(NUMX+1)*(NUMY+1)];

  }


  for (int m=0; m<=MMAX; m++)  {
    for (int i=1; i<=NUMR; i++)  {
      SC[m][i].zero();
      accum_cos0[m][i].zero();
      if (m>0) {
	accum_sin0[m][i].zero();
	SS[m][i].zero();
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

    int slave = 0;
    int request_id = 1;		// Begin with cosine case
    int M, ir;

    M = 0;			// Initial counters
    ir = 1;
    while (M<=MMAX) {
	
      // Send request to slave
      if (slave<numprocs-1) {
	  
	slave++;
	  
	MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	MPI_Send(&M,  1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	MPI_Send(&ir,  1, MPI_INT, slave, 11, MPI_COMM_WORLD);
	  
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
	
      if (slave == numprocs-1 && M<=MMAX) {
	  
	//
	// <Wait and receive>
	//
	  
	int retid = receive_eof();
	  
#ifdef DEBUG
	cerr.form("master in make_eof: received data and about to send orders to %d\n", retid);
#endif

	//
	// <Send new request>
	//
	  
	  
	MPI_Send(&request_id, 1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	MPI_Send(&M,  1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	MPI_Send(&ir,  1, MPI_INT, retid, 11, MPI_COMM_WORLD);
	  
#ifdef DEBUG
	cerr.form("master in make_eof: new orders sent to %d\n", retid);
#endif

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
    // <Wait for all slaves to return>
    //
      
    while (slave) {
      receive_eof();
      slave--;
    }
      
      
    //
    // <Tell slaves to continue>
    //
      
    request_id = -1;
    for (slave=1; slave < numprocs; slave++)
      MPI_Send(&request_id, 1, MPI_INT, slave, 11, MPI_COMM_WORLD);
      
  } else {

    int M, request_id, ir, icnt;

    while (1) {

      MPI_Recv(&request_id, 1, 
	       MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
				// Done!
      if (request_id<0) break;

      MPI_Recv(&M, 1, 
	       MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
      MPI_Recv(&ir, 1, 
	       MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	

#ifdef DEBUG
      cerr.form("Slave %d: received orders type=%d  M=%d  ir=%d\n",
		myid, request_id, M, ir);
#endif

      if (request_id) {

				// Complete symmetric part
    

	for (int i=1; i<rank2; i++) {
	  for (int j=i+1; j<=rank2; j++) {
	    SC[M][ir][j][i] = SC[M][ir][i][j];
	  }
	}
    
	double maxV = 0.0;
	for (int i=1; i<=rank2; i++) {
	  for (int j=i; j<=rank2; j++) {
	    tmp = fabs(SC[M][ir][i][j]);
	    if (tmp > maxV) maxV = tmp;
	  }
	}

	SC[M][ir] /= maxV;
    
    /*==========================*/
    /* Solve eigenvalue problem */
    /*==========================*/
    
#if defined(STURM)
	ev = Symmetric_Eigenvalues_MSRCH(SC[M][ir], ef, NORDER);

#elif defined(GHQL)
	ev = SC[M][ir].Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#else
	ev = SC[M][ir].Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#endif
      }
      else {

				// Complete symmetric part
    

	for (int i=1; i<rank2; i++) {
	  for (int j=i+1; j<=rank2; j++) {
	    SS[M][ir][j][i] = SS[M][ir][i][j];
	  }
	}
    
	double maxV = 0.0;
	for (int i=1; i<=rank2; i++) {
	  for (int j=i; j<=rank2; j++) {
	    tmp = fabs(SS[M][ir][i][j]);
	    if (tmp > maxV) maxV = tmp;
	  }
	}
	
	SS[M][ir] /= maxV;
    
    /*==========================*/
    /* Solve eigenvalue problem */
    /*==========================*/
    
#if defined(STURM)
	ev = Symmetric_Eigenvalues_MSRCH(SS[M][ir], ef, NORDER);

#elif defined(GHQL)
	ev = SS[M][ir].Symmetric_Eigenvalues(ef);
	ef = ef.Transpose();
#else
	ev = SS[M][ir].Symmetric_Eigenvalues(ef);
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

    if (mm) {
      mcos = cos(phi*mm) * facM * facZ;
      msin = sin(phi*mm) * facM * facZ;
    }
    else
      mcos = fac0 * facZ;

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

  double ccos, ssin, fac, facN;
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
  double fac;

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
    
double EmpCylSL::accumulated_dens_eval(double r, double z, double phi)
{
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

  double ccos, ssin, fac;
  int n, mm, v;

  for (mm=0; mm<=MMAX; mm++) {

    if (mm) {
      ccos = cos(phi*mm)/M_PI;
      ssin = sin(phi*mm)/M_PI;
    }
    else
      ccos = 0.5/M_PI;

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

  return 0.25*ans/M_PI;
}


  
void EmpCylSL::get_all(int mm, int ir, int nn, double r, double z, double phi,
		       double& p, double& d, double& fr, double& fz, double &fp)
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

  d += ccos *
    (
     densC[mm][n][ix  ][iy  ] * c00 +
     densC[mm][n][ix+1][iy  ] * c10 +
     densC[mm][n][ix  ][iy+1] * c01 +
     densC[mm][n][ix+1][iy+1] * c11 
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
  

  if (mm) {
    
    p += ssin *
      (
       potS[mm][n][ix  ][iy  ] * c00 +
       potS[mm][n][ix+1][iy  ] * c10 +
       potS[mm][n][ix  ][iy+1] * c01 +
       potS[mm][n][ix+1][iy+1] * c11 
       );

    d += ssin *
      (
       densS[mm][n][ix  ][iy  ] * c00 +
       densS[mm][n][ix+1][iy  ] * c10 +
       densS[mm][n][ix  ][iy+1] * c01 +
       densS[mm][n][ix+1][iy+1] * c11 
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
      
  }
  
  d *= 0.25/M_PI;

}


#include <stdio.h>

void EmpCylSL::dump_coefs(FILE *fout)
{

  fprintf(fout, "%4d %4d %4d\n", MMAX, NUMR, NZOF);

  for (int mm=0; mm<=MMAX; mm++) {

    fprintf(fout, "%4d\n", mm);

    for (int ir=1; ir<=NUMR; ir++)
      for (int n=0; n<rank2; n++)
	fprintf(fout, " %13.4e", accum_cos0[mm][ir][n]);
    fprintf(fout, "\n");

    for (int ir=1; ir<=NUMR; ir++)
      for (int n=0; n<rank2; n++)
	fprintf(fout, " %13.4e", accum_sin0[mm][ir][n]);
    fprintf(fout, "\n");


  }
}


void EmpCylSL::dump_coefs_binary(FILE *fout, double time)
{
  double p;

  eof_coefheader.time = time;
  eof_coefheader.mmax = MMAX;
  eof_coefheader.numr = NUMR;
  eof_coefheader.nzof = NZOF;

  fwrite(&eof_coefheader, sizeof(EOFCoefHeader), 1, fout);

  for (int mm=0; mm<=MMAX; mm++) {

    for (int ir=1; ir<=NUMR; ir++)
      for (int n=1; n<=rank2; n++)
	fwrite(&(p=accum_cos0[mm][ir][n]), sizeof(double), 1, fout);
    
    if (mm) {
      for (int ir=1; ir<=NUMR; ir++)
	for (int n=1; n<=rank2; n++)
	  fwrite(&(p=accum_sin0[mm][ir][n]), sizeof(double), 1, fout);
    }

  }
}



