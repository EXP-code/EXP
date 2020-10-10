#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>

#include <Vector.h>
#include <SphereSL.H>


int SphereSL::NUMR=800;
bool SphereSL::mpi=false;	// Initially off


SphereSL::SphereSL(SphericalModelTable* mod, int LMAX, int NMAX,
		   int CMAP, double Scale, bool COVAR)
{
  SphericalModelTable *model = mod;
  
  SLGridSph::mpi = mpi ? 1 : 0;

  double rmin = max<double>(mod->get_min_radius()*2.0, 
			    mod->get_max_radius()*1.0e-4);
  double rmax = mod->get_max_radius()*0.99;

  sl = new SLGridSph(LMAX, NMAX, NUMR, rmin, rmax, model, CMAP, Scale);

  lmax = LMAX;
  nmax = NMAX;
  compute_covar = COVAR;
  
  coefs_defined = false;
  rscl = 1.0;

  potd.setsize(0, LMAX, 1, NMAX);
  dpot.setsize(0, LMAX, 1, NMAX);
  dpt2.setsize(0, LMAX, 1, NMAX);
  dend.setsize(0, LMAX, 1, NMAX);
    
  legs.setsize(0, LMAX, 0, LMAX);
  dlegs.setsize(0, LMAX, 0, LMAX);
  d2legs.setsize(0, LMAX, 0, LMAX);

}

SphereSL::~SphereSL(void)
{
  delete sl;
}

void SphereSL::bomb(char *s)
{
  cerr << "ERROR from SphereSL: " << s << '\n';
  exit(-1);
}

void SphereSL::reset_coefs(void)
{
  if (expcoef.getnrows()>0 && expcoef.getncols()>0) expcoef.zero();
  if (compute_covar && cc.getnrows()>0 && cc.getncols()>0) cc.zero();
}


void SphereSL::accumulate(double x, double y, double z, double mass)
{
  double fac, fac1, fac2, fac4;
  double fac0=-4.0*M_PI;
  const double dsmall = 1.0e-20;

  if (!coefs_defined) {

    coefs_defined = true;

    expcoef.setsize(0, lmax*(lmax+2), 1, nmax);
    expcoef.zero();		// Need this?

    work1.setsize(1, nmax);

    if (compute_covar) {
      nmat = (lmax*(lmax+2)+1)*nmax;
      cc.setsize(1, nmat, 1, nmat);
      cw.setsize(1, nmat);
      cc.zero();

      work2.setsize(1, nmat);
    }

    factorial.setsize(0, lmax, 0, lmax);

    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	factorial[l][m] = sqrt( (0.5*l+0.25)/M_PI * 
				exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );
	if (m != 0) factorial[l][m] *= M_SQRT2;
      }
    }

    used = 0;
  }

  //======================
  // Compute coefficients 
  //======================

  double r2 = (x*x + y*y + z*z);
  double r = sqrt(r2) + dsmall;
  double costh = z/r;
  double phi = atan2(y,x);
  double rs = r/rscl;
	
  used++;

  sl->get_pot(potd, rs);

  legendre_R(lmax, costh, legs);

  // L loop
  for (int l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

    // M loop
    for (int m=0, moffset=0; m<=l; m++) {

      if (m==0) {
	fac = factorial[l][m] * legs[l][m];
	for (int n=1; n<=nmax; n++) {
	  fac4 = potd[l][n]*fac*fac0;
	  expcoef[loffset+moffset][n] += fac4 * mass;
	  if (compute_covar) {
	    cw[(loffset+moffset)*nmax + n] = fac4;
	  }
	}

	moffset++;
      }
      else {
	fac = factorial[l][m] * legs[l][m];
	fac1 = fac*cos(phi*m);
	fac2 = fac*sin(phi*m);
	for (int n=1; n<=nmax; n++) {
	  fac4 = potd[l][n]*fac0;
	  expcoef[loffset+moffset  ][n] += fac1 * fac4 * mass;
	  expcoef[loffset+moffset+1][n] += fac2 * fac4 * mass;
	  if (compute_covar) {
	    cw[(loffset+moffset  )*nmax + n] = fac1 * fac4;
	    cw[(loffset+moffset+1)*nmax + n] = fac2 * fac4;
	  }
	}
	moffset+=2;
      }
    }
  }

  if (compute_covar) {
    for (int n=1; n<=nmat; n++) {
      for (int nn=n; nn<=nmat; nn++) {
	cc[n][nn] += cw[n] * cw[nn] * mass;
      }
    }
  }

}

void SphereSL::make_coefs()
{
  if (mpi) {

    for (int l=0; l<=lmax*(lmax+2); l++) {
      MPI_Allreduce(&expcoef[l][1], &work1[1], nmax, MPI_DOUBLE,
		    MPI_SUM, MPI_COMM_WORLD);

      expcoef[l] = work1;
    }

    if (compute_covar) {
      for (int n=1; n<=nmat; n++) {
	MPI_Allreduce(&cc[n][1], &work2[1], nmat, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	
	cc[n] = work2;
      }
    }
  }
}


void SphereSL::dens_pot_eval(double r, double costh, double phi,
			     double& dens0, double& dens, 
			     double& potl0, double& potl,
			     int L1, int L2)
{
  double fac1, cosm, sinm;

  fac1 = factorial[0][0];

  sl->get_dens(dend, r/rscl);
  sl->get_pot (potd, r/rscl);

  legendre_R(lmax, costh, legs, dlegs);

  dens0 = fac1 * expcoef[0]*dend[0];
  potl0 = fac1 * expcoef[0]*potd[0];

  dens = 0.0;
  potl = 0.0;

  // L loop
  for (int l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
    if (l<L1 || l>L2) continue;

    // M loop
    for (int m=0, moffset=0; m<=l; m++) {
      fac1 = factorial[l][m];
      if (m==0) {
	dens += fac1*legs[l][m]* (expcoef[loffset+moffset] * dend[l]);

	potl += fac1*legs[l][m]* (expcoef[loffset+moffset] * potd[l]);

	moffset++;
      }
      else {
	cosm = cos(phi*m);
	sinm = sin(phi*m);

	dens += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * dend[l]*cosm + 
	    expcoef[loffset+moffset+1] * dend[l]*sinm );

	potl += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );

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


void SphereSL::pot_force_eval(double r, double costh, double phi,
			      double& potl,
			      double& potr, double& pott, double& potp,
			      int L1, int L2)
{
  double fac1, cosm, sinm;
  double sinth = -sqrt(fabs(1.0 - costh*costh));

  fac1 = factorial[0][0];

  sl->get_pot  (potd, r/rscl);
  sl->get_force(dpot, r/rscl);

  legendre_R(lmax, costh, legs, dlegs);

  potl = fac1 * expcoef[0]*potd[0];
  potr = fac1 * expcoef[0]*dpot[0];
  pott = 0.0;
  potp = 0.0;

  // L loop
  for (int l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
    if (l<L1 || l>L2) continue;

    // M loop
    for (int m=0, moffset=0; m<=l; m++) {
      fac1 = factorial[l][m];
      if (m==0) {
	potl += fac1*legs[l][m]* (expcoef[loffset+moffset] * potd[l]);
	dpot += fac1*legs[l][m]* (expcoef[loffset+moffset] * dpot[l]);
	pott += fac1*dlegs[l][m]* (expcoef[loffset+moffset] * potd[l]);

	moffset++;
      }
      else {
	cosm = cos(phi*m);
	sinm = sin(phi*m);

	potl += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );
	dpot += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * dpot[l]*cosm + 
	    expcoef[loffset+moffset+1] * dpot[l]*sinm );
	pott += fac1*dlegs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );
	potp += fac1*legs[l][m] * m *
	  (-expcoef[loffset+moffset]   * potd[l]*sinm + 
	    expcoef[loffset+moffset+1] * potd[l]*cosm );

	moffset +=2;
      }
    }
  }

  // double densfac = 1.0/(rscl*rscl*rscl) * 0.25/M_PI;
  double potlfac = 1.0/rscl;

  potl  *= potlfac;
  potr  *= potlfac/rscl;
  pott  *= potlfac*sinth;
  potp  *= potlfac;
}


void SphereSL::all_eval(double r, double costh, double phi,
			double& den0, double& den1,
			double& pot0, double& pot1,
			double& potr, double& pott, double& potp,
			int L1, int L2)
{
  double fac1, cosm, sinm;
  double sinth = -sqrt(fabs(1.0 - costh*costh));

  fac1 = factorial[0][0];

  sl->get_dens (dend, r/rscl);
  sl->get_pot  (potd, r/rscl);
  sl->get_force(dpot, r/rscl);

  legendre_R(lmax, costh, legs, dlegs);

  den0 = fac1 * expcoef[0]*dend[0];
  pot0 = fac1 * expcoef[0]*potd[0];
  potr = fac1 * expcoef[0]*dpot[0];
  den1 = 0.0;
  pot1 = 0.0;
  pott = 0.0;
  potp = 0.0;

  // L loop
  for (int l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
    if (l<L1 || l>L2) continue;

    // M loop
    for (int m=0, moffset=0; m<=l; m++) {
      fac1 = factorial[l][m];
      if (m==0) {
	den1 += fac1*legs[1][m] * (expcoef[loffset+moffset] * dend[l]);
	pot1 += fac1*legs[l][m] * (expcoef[loffset+moffset] * potd[l]);
	dpot += fac1*legs[l][m] * (expcoef[loffset+moffset] * dpot[l]);
	pott += fac1*dlegs[l][m]* (expcoef[loffset+moffset] * potd[l]);

	moffset++;
      }
      else {
	cosm = cos(phi*m);
	sinm = sin(phi*m);

	den1 += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * dend[l]*cosm + 
	    expcoef[loffset+moffset+1] * dend[l]*sinm );
	pot1 += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );
	dpot += fac1*legs[l][m]*
	  ( expcoef[loffset+moffset]   * dpot[l]*cosm + 
	    expcoef[loffset+moffset+1] * dpot[l]*sinm );
	pott += fac1*dlegs[l][m]*
	  ( expcoef[loffset+moffset]   * potd[l]*cosm + 
	    expcoef[loffset+moffset+1] * potd[l]*sinm );
	potp += fac1*legs[l][m] * m *
	  (-expcoef[loffset+moffset]   * potd[l]*sinm + 
	    expcoef[loffset+moffset+1] * potd[l]*cosm );

	moffset +=2;
      }
    }
  }

  double densfac = 1.0/(rscl*rscl*rscl) * 0.25/M_PI;
  double potlfac = 1.0/rscl;

  den0  *= densfac;
  den1  *= densfac;
  pot0  *= potlfac;
  pot1  *= potlfac;
  potr  *= potlfac/rscl;
  pott  *= potlfac*sinth;
  potp  *= potlfac;
}


#define MINEPS 1.0e-10

void SphereSL::legendre_R(int lmax, double x, Matrix& p)
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
      if (std::isnan(p[m][m]))
	cerr << "legendre_R: p[" << m << "][" << m << "]: pll=" << pll << "\n";
      fact += 2.0;
    }
  }
  
  for (m=0; m<lmax; m++) {
    pl2 = p[m][m];
    p[m+1][m] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[l][m] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      if (std::isnan(p[l][m]))
	cerr << "legendre_R: p[" << l << "][" << m << "]: pll=" << pll << "\n";
      
      pl2 = pl1;
      pl1 = pll;
    }
  }
  
  if (std::isnan(x))
    cerr << "legendre_R: x\n";
  for(l=0; l<=lmax; l++)
    for (m=0; m<=l; m++)
      if (std::isnan(p[l][m]))
	cerr << "legendre_R: p[" << l << "][" << m << "] lmax=" 
	     << lmax << "\n";
  
}

void SphereSL::legendre_R(int lmax, double x, Matrix &p, Matrix &dp)
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

void SphereSL::legendre_R(int lmax, double x, Matrix &p, Matrix &dp, Matrix& d2p)
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
  
  for (l=0; l<=lmax; l++) {
    for (m=0; m<=l; m++)
      d2p[l][m] = -somx2*(2.0*x*dp[l][m] - p[l][m]*(somx2*m*m + l*(l+1)));
  }
  
}


void SphereSL::install_coefs(Matrix& newcoef)
{
  if (!coefs_defined) {

    coefs_defined = true;
    
    expcoef.setsize(0, lmax*(lmax+2), 1, nmax);
    expcoef.zero();		// Need this?

    work1.setsize(1, nmax);

    if (compute_covar) {
      nmat = (lmax*(lmax+2)+1)*nmax;
      cc.setsize(1, nmat, 1, nmat);
      cw.setsize(1, nmat);
      cc.zero();

      work2.setsize(1, nmat);
    }

    factorial.setsize(0, lmax, 0, lmax);

    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	factorial[l][m] = sqrt( (0.5*l+0.25)/M_PI * 
				exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );
	if (m != 0) factorial[l][m] *= M_SQRT2;
      }
    }

    used = 0;
  }

  // Sanity check
  if (newcoef.getrhigh() != expcoef.getrhigh() ||
      newcoef.getrlow()  != expcoef.getrlow()  ||
      newcoef.getchigh() != expcoef.getchigh() ||
      newcoef.getclow()  != expcoef.getclow()  )
    {
      cerr << "SphereSL: can not install coefficients, dimension mismatch\n";
      return;
    }

  // Do the assignment
  expcoef = newcoef;
}


void SphereSL::dump_coefs(double time, ostream& out)
{
  ostringstream sout;
  sout << "SphereSL";

  char buf[64];
  for (int i=0; i<64; i++) {
    if (i<sout.str().length())  buf[i] = sout.str().c_str()[i];
    else                        buf[i] = '\0';
  }

  out.write((char *)&buf,64*sizeof(char));
  out.write((char *)&time , sizeof(double));
  out.write((char *)&rscl,  sizeof(double));
  out.write((char *)&nmax,  sizeof(int));
  out.write((char *)&lmax,  sizeof(int));

  for (int ir=1; ir<=nmax; ir++) {
    for (int l=0; l<=lmax*(lmax+2); l++)
      out.write((char *)&expcoef[l][ir], sizeof(double));
  }

}
