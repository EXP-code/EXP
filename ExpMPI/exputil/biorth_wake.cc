#pragma implementation

using namespace std;

//***************************************************************************
// Wake stuff
//***************************************************************************

#include <string>
#include <iostream>
#include <math.h>
#include <Vector.h>
#include <biorth.h>
#include <biorth_wake.h>

double factrl(int n);
double plgndr(int l, int m, double x);

int BiorthWake::iterDef = 200;
double BiorthWake::tolDef = 1.0e-6;
int BiorthWake::ndim = 3;

BiorthWake::BiorthWake(AxiSymBiorth *BIO, int LMAX, int NMAX)
{
  bio = BIO;
  lmax = LMAX;
  nmax = NMAX;
  
  coefs_defined = false;
  rscl = 1.0;

  norm_grid.setsize(0, nmax-1, 0, lmax);
  for (int l=0; l<=lmax; l++) {
    for (int n=1; n<=nmax; n++) norm_grid[n-1][l] = sqrt(bio->norm(n-1, l));
  }

				// Orientation parameters
  I = Complex(0.0, 1.0);
  iter = iterDef;
  tol = tolDef;
  init_orientation = false;
}


BiorthWake::~BiorthWake(void)
{
  if (init_orientation) {
    delete [] param;
    delete [] (psum+1);
    delete [] (ptry+1);

    const int num=3;
    for (int i=1; i<=num+1; i++) delete [] (ambp[i]+1);
    delete [] (ambp+1);
    delete [] (amby+1);
  }
}


void BiorthWake::bomb(char *s)
{
  cerr << "ERROR from BiorthWake(" << bio->BiorthID << "): " << s << '\n';
  exit(-1);
}

void BiorthWake::reset_coefs(void)
{
  expcoef.zero();
}


void BiorthWake::accumulate(double x, double y, double z, double mass)
{
  switch (bio->get_dof()) {
  case 2:
    accumulate_2d(x, y, mass);
    break;
  case 3:
    accumulate_3d(x, y, z, mass);
    break;
  default:
    bomb("accumulate: only 2 and 3 accumulation");
  }
}

void BiorthWake::reconstruct(double r, double costh, double phi,
			     double& dens0, double& dens, 
			     double& potl0, double& potl,
			     int L1, int L2)
{
  switch (bio->get_dof()) {
  case 2:
    reconstruct_2d(r, phi, dens0, dens, potl0, potl, L1, L2);
    break;
  case 3:
    reconstruct_3d(r, costh, phi, dens0, dens, potl0, potl, L1, L2);
    break;
  default:
    bomb("reconstruct: only 2 and 3 reconstruction");
  }
}


void BiorthWake::accumulate_2d(double x, double y, double mass)
{
  int l,m,n,moffset;
  double r2,r,rs,fac,fac1,fac2,phi;
  double fac0=2.0*M_PI;
  const double dsmall = 1.0e-20;

  if (!coefs_defined) {

    coefs_defined = true;

    expcoef.setsize(0, 2*lmax, 1, nmax);
    expcoef.zero();		// Need this?

    used = 0;
  }

  /*======================*/
  /* Compute coefficients */
  /*======================*/

  r2 = (x*x + y*y);
  r = sqrt(r2) + dsmall;
  used++;
  phi = atan2(y,x);
  rs = r/rscl;
	
  /*		m loop */
  for (m=0; m<=lmax; m++) {
    if (m==0) {
      for (n=1; n<=nmax; n++)
	expcoef[m][n] += bio->potlR(n, m, rs)/norm_grid[n-1][l] * fac0 * mass;
	moffset++;
      }
      else {
	fac1 = fac*cos(phi*m);
	fac2 = fac*sin(phi*m);
	for (n=1; n<=nmax; n++) {
	  fac = bio->potlR(n, m, rs)/norm_grid[n-1][m] * fac0 * mass;
	  expcoef[2*(m-1)+1][n] += fac * fac1;
	  expcoef[2*(m-1)+2][n] += fac * fac2;
	}
      }
  }
}

  
void BiorthWake::reconstruct_2d(double r, double phi,
				double& dens0, double& dens, 
				double& potl0, double& potl,
				int L1, int L2)
{
  int l,m;
  double cosm,sinm;
  double fac0 = 1.0/(2.0*M_PI);

  dens = 0.0;
  dens0 = bio->get_dens(r/rscl,0,expcoef[0]) * fac0;
  potl0 = bio->get_potl(r/rscl,0,expcoef[0]) * fac0;
  potl = 0.0;

  /*		m loop */
  for (m=0; m<=lmax; m++) {
    if (l<L1 || l>L2) continue;
    if (m==0) {
      dens += bio->get_dens(r/rscl, m, expcoef[m]) * fac0;
      potl += bio->get_potl(r/rscl, m, expcoef[m]) * fac0;
    }
    else {
      cosm = cos(phi*m);
      sinm = sin(phi*m);
      dens += 2.0 * fac0 *
	  bio->get_dens(r/rscl, m, expcoef[2*(m-1)+1])*cosm + 
	  bio->get_dens(r/rscl, m, expcoef[2*(m-1)+2])*sinm ;
      potl += 2.0 * fac0 *
	  bio->get_potl(r/rscl, m, expcoef[2*(m-1)+1])*cosm + 
	  bio->get_potl(r/rscl, m, expcoef[2*(m-1)+2])*sinm ;
    }
  }

  double densfac = 1.0/(rscl*rscl*rscl) * 0.5/M_PI;
  double potlfac = 1.0/rscl;

  dens0 *= densfac;
  dens  *= densfac;
  potl0 *= potlfac;
  potl  *= potlfac;
}


void BiorthWake::accumulate_3d(double x, double y, double z, double mass)
{
  int l,loffset,moffset,m,n;
  double r2,r,rs,fac,fac1,fac2,costh,phi;
  double fac0=4.0*M_PI;
  const double dsmall = 1.0e-20;

  if (!coefs_defined) {

    coefs_defined = true;

    expcoef.setsize(0, lmax*(lmax+2), 1, nmax);
    expcoef.zero();		// Need this?

    factorial.setsize(0, lmax, 0, lmax);

    for (l=0; l<=lmax; l++) {
      for (m=0; m<=l; m++) 
	factorial[l][m] = factrl(l-m)/factrl(l+m);
    }

    used = 0;
  }

  /*======================*/
  /* Compute coefficients */
  /*======================*/

  r2 = (x*x + y*y + z*z);
  r = sqrt(r2) + dsmall;
  used++;
  costh = z/r;
  phi = atan2(y,x);
  rs = r/rscl;
	
  /*		l loop */
  for (l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

    /*		m loop */
    for (m=0, moffset=0; m<=l; m++) {
      if (m==0) {
	fac = plgndr(l,m,costh);
	for (n=1; n<=nmax; n++)
	  expcoef[loffset+moffset][n] += bio->potlR(n, l, rs)*fac*mass*
	    fac0/norm_grid[n-1][l];
	moffset++;
      }
      else {
	fac = plgndr(l,m,costh);
	fac1 = fac*cos(phi*m);
	fac2 = fac*sin(phi*m);
	for (n=1; n<=nmax; n++) {
	  expcoef[loffset+moffset][n] += bio->potlR(n, l, rs)*fac1*mass*
	    fac0/norm_grid[n-1][l];
	  expcoef[loffset+moffset+1][n] += bio->potlR(n, l, rs)*fac2*mass*
	    fac0/norm_grid[n-1][l];
	}
	moffset+=2;
      }
    }
  }
}

  
void BiorthWake::reconstruct_3d(double r, double costh, double phi,
			     double& dens0, double& dens, 
			     double& potl0, double& potl,
			     int L1, int L2)
{
  int l,m,loffset,moffset;
  double fac1,fac2,cosm,sinm;

  dens = 0.0;
  fac1 = 0.25/M_PI;
  dens0 = fac1 * bio->get_dens(r/rscl,0,expcoef[0]);
  potl0 = fac1 * bio->get_potl(r/rscl,0,expcoef[0]);
  potl = 0.0;

  /*		l loop */
  for (l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
    if (l<L1 || l>L2) continue;
    /*		m loop */
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (0.5*l+0.25)/M_PI;
      if (m==0) {
	dens += fac1*plgndr(l,m,costh)*
	  bio->get_dens(r/rscl,l,expcoef[loffset+moffset]);
	potl += fac1*plgndr(l,m,costh)*
	  bio->get_potl(r,l,expcoef[loffset+moffset]);
	moffset++;
      }
      else {
	fac2 = fac1 * factrl(l-m)/factrl(l+m);
	cosm = cos(phi*m);
	sinm = sin(phi*m);
	dens += 2.0 *fac2*plgndr(l,m,costh)*( 
		    bio->get_dens(r/rscl,l,expcoef[loffset+moffset])*cosm + 
		    bio->get_dens(r/rscl,l,expcoef[loffset+moffset+1])*sinm );
	potl += 2.0 *fac2*plgndr(l,m,costh)*( 
		    bio->get_potl(r/rscl,l,expcoef[loffset+moffset])*cosm + 
		    bio->get_potl(r/rscl,l,expcoef[loffset+moffset+1])*sinm );
	moffset +=2;
      }
    }
  }

  double densfac = 1.0/(rscl*rscl*rscl) * 0.25/M_PI;
  double potlfac = 1.0/rscl;

  dens0 *= densfac;
  dens  *= densfac;
  potl0 *= potlfac;
  potl  *= potlfac;
}
