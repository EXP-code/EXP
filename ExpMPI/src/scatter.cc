
/*
  void scatter_mfp()

  Turned on by setting flag "scatter".

*/

static char rcsid[] = "$Id$";

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>
#include <string>

#include <Vector.h>

#include "expand.h"


extern "C" void indexx(int n, double *arrin, int *indx);

extern "C" void scatter_mfp(void)
{
  static bool firstime = true;
  static double tau;
  const int tautab = 100;
  static double *dtau, *dtau1, dr;
  static int *indx;
  static ACG *gen;
  static Uniform *unif;
  static Normal *gaus;
  static int cntacc;

  if (!myid) return;

				// Initialization
  if (firstime) {
    cerr << "Isotropic scattering on with Tau=" << tauscat << "\n";
    firstime = false;

    dtau  = new double [tautab];
    dtau1 = new double [tautab];
    dr = rmax/tautab;

    indx = new int [nbodmax] - 1;

    gen = new ACG(11+myid, 20);
    unif = new Uniform(0.0, 1.0, gen);
    gaus = new Normal(0.0, 1.0, gen);

    cntacc = 0;
  }

				// Clean table
  for (int j=0; j<tautab; j++) dtau[j] = dtau1[j] = 0.0;

				// Accumulate table
  int ind;
  double rrr, vvv;
  for (int i=1; i<=nbodies; i++) {
    ind = (int)(rr[i]/dr);
    if (ind>=tautab) continue;
    dtau1[ind] += mass[i];
  }

				// Distribute table
  MPI_Allreduce(dtau1, dtau, tautab, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);

				// Compute density from mass per bin
  for (int j=0; j<tautab; j++)
    dtau[j] /= 4.0*M_PI/3.0*(pow(dr*(j+1), 3.0) - pow(dr*j, 3.0));


				// Compute index
  indexx(nbodies, rr, indx);

  Three_Vector vcom, vrel, vfnl;

  int i, k, cnt=0;
#ifdef DEBUG  
  double rm, rp, rtst1, rtst=0.0;
#endif
  for (int j=1; j<=nbodies; j++) {
    i = indx[j];

    if (freeze_particle(i)) continue;

    ind = (int)(rr[i]/dr);
    if (ind>=tautab) continue;

    mfp[i] += dtau[ind] * sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]) 
      * dtime;

    if (1.0 - exp(-mfp[i]/tauscat)>(*unif)()) {

				// Initialize optical depth
      mfp[i] = 0.0;

				// Choose a buddy
      if (j==1)
	k = indx[2];
#ifdef DEBUG      
      else if ((rp=fabs(rr[i]-rr[indx[j+1]])) > (rm=fabs(rr[i] - rr[indx[j-1]])))
#else

      else if (fabs(rr[i]-rr[indx[j+1]]) > fabs(rr[i] - rr[indx[j-1]]))
#endif
	k = indx[j-1];
      else
	k = indx[j+1];

#ifdef DEBUG
      // if (j>1 && rtst<(rtst1=min<double>(rm, rp))) rtst = rtst1;
      if (j>1 && rtst<(rtst1=(rm<rp ? rm : rp))) rtst = rtst1;
#endif

      vcom[1] = 0.5*(vx[i] + vx[k]);
      vcom[2] = 0.5*(vy[i] + vy[k]);
      vcom[3] = 0.5*(vz[i] + vz[k]);

      vrel[1] = vx[k] - vx[i];
      vrel[2] = vy[k] - vy[i];
      vrel[3] = vz[k] - vz[i];

				// Choose a random direction for velocity
      vfnl[1] = (*gaus)();
      vfnl[2] = (*gaus)();
      vfnl[3] = (*gaus)();

      vfnl *= sqrt(vrel*vrel)/sqrt(vfnl*vfnl);

				// To lab frame
      vx[i] = vcom[1] + 0.5*vfnl[1];
      vy[i] = vcom[2] + 0.5*vfnl[2];
      vz[i] = vcom[3] + 0.5*vfnl[3];

      vx[k] = vcom[1] - 0.5*vfnl[1];
      vy[k] = vcom[2] - 0.5*vfnl[2];
      vz[k] = vcom[3] - 0.5*vfnl[3];

      cnt++;
    }
  }

				// Keep count of mean # of interactions
  int sum=0;
  MPI_Reduce(&cnt, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_SLAVE);
  if (myid==1) {
    
    cntacc += sum;
    if (!(this_step%nscat)) {
      string file = (string)homedir + "scatter.log";
      ofstream out(file.c_str(), ios::out | ios::app);

      if (this_step == 0)
	out << "# "
	    << setw(6) << this_step
	    << setw(15) << (double)cntacc
	    << endl;
      else
	out << "# "
	    << setw(6) << this_step
	    << setw(15) << (double)cntacc/nscat
	    << endl;

      for (int j=0; j<tautab; j++)
	out << "     " 
	    << setw(15) << dr*(j+1)
	    << setw(15) << dtau[j]
	    << endl;

      cntacc = 0;
    }
  }

#ifdef DEBUG
  cerr.form("Process %d: rmax=%f\n", myid, rtst);
#endif
}
