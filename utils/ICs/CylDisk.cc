#include "CylDisk.H"

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

  double r;

  for (int i=0; i<NUMR; i++) {
    if (rlog) {
      rr[i] = log(rmin) + dr*i;
      r = exp(rr[i]);
    } else {
      r = rr[i] = rmin + dr*i;
    }

    double dens = 0.0, pot = 0.0;
    double p0, d0, p, fr, fz, fp;
    
    d->accumulated_eval(r, 0.0, 0.0, p0, p, fr, fz, fp);
    d->accumulated_dens_eval(r, 0.0, 0.0, d0);

    dd[i] = d0;
    pp[i] = p0;
  }

  mm[0] = 0.0;
  for (int i=1; i<NUMR; i++) {
    if (rlog)
      mm[i] = mm[i-1] + 
	(dd[i-1]*exp(2.0*rr[i-1]) + dd[i]*exp(2.0*rr[i]))*dr;
    else
      mm[i] = mm[i-1] + 
	(dd[i-1]*rr[i-1] + dd[i]*rr[i])*dr;
  }

  double norm = 1.0/mm[NUMR-1];
  for (int i=0; i<NUMR; i++) mm[i] *= norm;

}
