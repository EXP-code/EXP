#include "CylDisk.h"

int   CylDisk::NUMR = 200;
bool  CylDisk::RLOG = true;

void CylDisk::setup_model(double Rmin, double Rmax)
{
  rmin = Rmin;
  rmax = Rmax;

  rr = vector<double>(NUMR);
  dd = vector<double>(NUMR);
  mm = vector<double>(NUMR);
  pp = vector<double>(NUMR);

  bool rlog = RLOG;
  if (rmin <= 0.0) rlog = false;

  if (rlog)
    dr = (log(rmax) - log(rmin))/NUMR;
  else
    dr = (rmax - rmin)/NUMR;

  for (int i=0; i<NUMR; i++) {
    rr[i] = rmin + dr*i;
    if (rlog) rr[i] = exp(rr[i]);

    double dens = 0.0, pot = 0.0;
    double p0, d0, p, fr, fz, fp;
    
    d->accumulated_eval(rr[i], 0.0, 0.0, p0, p, fr, fz, fp);
    d->accumulated_dens_eval(rr[i], 0.0, 0.0, d0);

    dd[i] = d0;
    pp[i] = p0;
  }

  mm[0] = 0.0;
  for (int i=1; i<NUMR; i++)
    mm[i] = mm[i-1] + M_PI*(dd[i-1]*rr[i-1] + dd[i]*rr[i])*(rr[i] - rr[i-1]);
}
