#include "expand.h"
#include <SphericalBasisMixtureSL.H>
#include <SphereTwoCenter.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

SphericalBasisMixtureSL::SphericalBasisMixtureSL
(string& line, SphereTwoCenter* A1, CenterType ctr1) : SphericalBasis(line)
{
  A = A1;
  ctr = ctr1;
  setup();
}


void * SphericalBasisMixtureSL::determine_coefficients_thread(void * arg)
{
  int l, i, loffset, moffset, m, n, nn;
  double r, r2, rs, fac1, fac2, costh, phi, mass;
  double facs1=0.0, facs2=0.0, fac0=4.0*M_PI;
  double xx, yy, zz;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = component->Adiabatic();
  
  use[id] = 0;

  for (i=nbeg; i<nend; i++) {

    if (cC->freeze(*(cC->Part(i)))) continue;

    mass = cC->Mass(i) * adb;

    xx = cC->Pos(i, 0, Component::Local) - A->center[0];
    yy = cC->Pos(i, 1, Component::Local) - A->center[1];
    zz = cC->Pos(i, 2, Component::Local) - A->center[2];

    r2 = (xx*xx + yy*yy + zz*zz);
    r = sqrt(r2) + DSMALL;
      
    if (r<=rmax) {
      use[id]++;
      costh = zz/r;
      phi = atan2(yy,xx);
      rs = r/scale;
	
      
      legendre_R(Lmax, costh, legs[id]);
      sinecosine_R(Lmax, phi, cosm[id], sinm[id]);

      get_potl(Lmax, nmax, rs, potd[id], id);

      /*		l loop */
      for (l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
	/*		m loop */
	for (m=0, moffset=0; m<=l; m++) {
	  if (m==0) {
	    if (selector && compute)
	      facs1 = legs[id][l][m]*legs[id][l][m]*mass;
	    for (n=1; n<=nmax; n++) {
	      expcoef0[id][loffset+moffset][n] += potd[id][l][n]*legs[id][l][m]*mass*
		fac0/normM[l][n];

	      if (selector && compute) {
		pthread_mutex_lock(&cc_lock);
		for (nn=n; nn<=nmax; nn++)
		  cc1[loffset+moffset][n][nn] += potd[id][l][n]*potd[id][l][nn]*facs1/(normM[l][n]*normM[l][nn]);
		pthread_mutex_unlock(&cc_lock);
	      }
	    }
	    moffset++;
	  }
	  else {
	    fac1 = legs[id][l][m]*cosm[id][m];
	    fac2 = legs[id][l][m]*sinm[id][m];
	    if (selector && compute) {
	      facs1 = fac1*fac1*mass;
	      facs2 = fac2*fac2*mass;
	    }
	    for (n=1; n<=nmax; n++) {
	      expcoef0[id][loffset+moffset][n] += potd[id][l][n]*fac1*mass*fac0/normM[l][n];
	      expcoef0[id][loffset+moffset+1][n] += potd[id][l][n]*fac2*mass*fac0/normM[l][n];

	      if (selector && compute) {
		pthread_mutex_lock(&cc_lock);
		for (nn=n; nn<=nmax; nn++) {
		  cc1[loffset+moffset][n][nn] += 
		    potd[id][l][n]*potd[id][l][nn]*facs1/(normM[l][n]*normM[l][nn]);
		  cc1[loffset+moffset+1][n][nn] +=
		    potd[id][l][n]*potd[id][l][nn]*facs2/(normM[l][n]*normM[l][nn]);
		}
		pthread_mutex_unlock(&cc_lock);
	      }
	      
	    }
	    moffset+=2;
	  }
	}
      }
    }
  }

  return (NULL);
}


void * SphericalBasisMixtureSL::determine_acceleration_and_potential_thread(void * arg)
{
  int l,loffset,moffset,m,ioff;
  double r,rs,r0=0.0,fac,fac1,fac2,fac3,fac4,costh,phi,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps,facp,facdp;
  double dfac=0.25/M_PI;

  Particle pt;
  double xx, yy, zz, mfactor;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  for (int i=nbeg; i<nend; i++) {

    if (use_external) {
      // Get the position in inertial coords
      cC->Pos(pt.pos, i, Component::Inertial);
      // Convert this to component coords
      component->ConvertPos(pt.pos, Component::Local);
    } else
      // Get the position in local com coords
      cc->Pos(pt.pos, i, Component::Local);

    if (ctr == ej) mfactor = 1.0 - A->mixture(pt);
    else           mfactor = A->mixture(pt);

    fac1 = dfac * mfactor;

				// Now, move to current center . . .
    xx = pt.pos[0] - A->center[0];
    yy = pt.pos[1] - A->center[1];
    zz = pt.pos[2] - A->center[2];

    r = sqrt(xx*xx + yy*yy + zz*zz) + DSMALL;
    costh = zz/r;
    rs = r/scale;
    phi = atan2(yy, xx);

    dlegendre_R(Lmax, costh, legs[id], dlegs[id]);
    sinecosine_R(Lmax, phi, cosm[id], sinm[id]);

    if (r>rmax) {
      ioff = 1;
      r0 = r;
      r = rmax;
      rs = r/scale;
    }
    else
      ioff = 0;


    get_dpotl(Lmax, nmax, rs, potd[id], dpot[id], id);
    get_pot_coefs_safe(0, expcoef[0], &p, &dp, potd[id], dpot[id]);
    if (ioff) {
      p *= rmax/r0;
      dp = -p/r0;
    }
    potl = fac1*p;
    potr = fac1*dp;
    pott = potp = 0.0;
      
    //		l loop
    //		------
    for (l=1, loffset=1; l<=Lmax; loffset+=(2*l+1), l++) {

				// Suppress L=1 terms?
      if (NO_L1  && l==1) continue;

				// Suppress odd L terms?
      if (EVEN_L && (l/2)*2 != l) continue;

      //		m loop
      //		------
      for (m=0, moffset=0; m<=l; m++) {
	fac1 = (2.0*l+1.0)/(4.0*M_PI) * mfactor;
	if (m==0) {
	  fac2 = fac1*legs[id][l][m];
	  get_pot_coefs_safe(l, expcoef[loffset+moffset], &p, &dp,
			     potd[id], dpot[id]);
	  if (ioff) {
	    p *= pow(rmax/r0,(double)(l+1));
	    dp = -p/r0 * (l+1);
	  }
	  potl += fac2*p;
	  potr += fac2*dp;
	  pott += fac1*dlegs[id][l][m]*p;
	  moffset++;
	}
	else {
	  fac2 = 2.0 * fac1 * factorial[l][m];
	  fac3 = fac2 * legs[id][l][m];
	  fac4 = fac2 * dlegs[id][l][m];
	  get_pot_coefs_safe(l, expcoef[loffset+moffset], &pc, &dpc,
			     potd[id], dpot[id]);
	  get_pot_coefs_safe(l, expcoef[loffset+moffset+1] ,&ps, &dps,
			     potd[id], dpot[id]);
	  if (ioff) {
	    facp = pow(rmax/r0,(double)(l+1));
	    facdp = -1.0/r0 * (l+1);
	    pc *= facp;
	    ps *= facp;
	    dpc = pc*facdp;
	    dps = ps*facdp;
	  }
	  potl += fac3*(pc*cosm[id][m] + ps*sinm[id][m]);
	  potr += fac3*(dpc*cosm[id][m] + dps*sinm[id][m]);
	  pott += fac4*(pc*cosm[id][m] + ps*sinm[id][m]);
	  potp += fac3*(-pc*sinm[id][m] + ps*cosm[id][m])*m;
	  moffset +=2;
	}
      }
    }

    fac = xx*xx + yy*yy;

    potr /= scale*scale;
    potl /= scale;
    pott /= scale;
    potp /= scale;

    cC->AddAcc(i, 0, -(potr*xx/r - pott*xx*zz/(r*r*r) ) );
    cC->AddAcc(i, 1, -(potr*yy/r - pott*yy*zz/(r*r*r) ) );
    cC->AddAcc(i, 2, -(potr*zz/r + pott*fac/(r*r*r))  );
    if (fac > DSMALL2) {
      cC->AddAcc(i, 0,  potp*yy/fac);
      cC->AddAcc(i, 1, -potp*xx/fac);
    }
    if (use_external)
      cC->AddPotExt(i, potl);
    else
      cC->AddPot(i, potl);

  }

  return (NULL);
}

void SphericalBasisMixtureSL::get_dpotl(int lmax, int nmax, double r, 
				Matrix& p, Matrix& dp, int tid)
{
  A->ortho->get_pot(p, r);
  A->ortho->get_force(dp, r);
}

void SphericalBasisMixtureSL::get_potl(int lmax, int nmax, double r, Matrix& p, 
			       int tid)
{
  A->ortho->get_pot(p, r);
}

void SphericalBasisMixtureSL::get_dens(int lmax, int nmax, double r, Matrix& p, 
			       int tid)
{
  A->ortho->get_dens(p, r);
}

void SphericalBasisMixtureSL::get_potl_dens(int lmax, int nmax, double r, 
				    Matrix& p, Matrix& d, int tid)
{
  A->ortho->get_pot(p, r);
  A->ortho->get_dens(d, r);
}

