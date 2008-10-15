#include "expand.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include <Particle.H>
#include <CBrockDisk.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

CBrockDisk::CBrockDisk(string& line) : AxisymmetricBasis(line)
{
  id = "Clutton-Brock two-dimensional disk";

  dof = 2;			// Two degrees of freedom

  rmax = 100.0;
  scale = 1.0;
  Lmax = 4;
  nmax = 10;
  self_consistent = true;
  selector = false;
  coef_dump = true;

  initialize();


  expcoef.setsize(0,2*Lmax+1,1,nmax);
  expcoef1.setsize(0,2*Lmax+1,1,nmax);

  expcoef0 = new Matrix [nthrds];
  if (!expcoef0) bomb("problem allocating <expcoef0>");

  for (int i=0; i<nthrds; i++)
    expcoef0[i].setsize(0, Lmax*(Lmax+2), 1, nmax);

  if (selector) {
    cc = new Matrix [Lmax*(Lmax+2)+1];
    if (!cc) bomb("problem allocating <cc>");

    for (int l=0; l<=Lmax*(Lmax+2); l++)
      cc[l].setsize(1, nmax, 1, nmax);
  
    cc1 = new Matrix [Lmax*(Lmax+2)+1];
    if (!cc1) bomb("problem allocating <cc1>");
    
    for (int l=0; l<=Lmax*(Lmax+2); l++)
      cc1[l].setsize(1, nmax, 1, nmax);
  
    pthread_mutex_init(&cc_lock, NULL);
  }
  // Allocate and compute normalization matrix

  normM.setsize(0,Lmax,1,nmax);
  dend.setsize(0,Lmax,1,nmax);
  work.setsize(0,Lmax+1,1,nmax);
  
  potd  = new Matrix [nthrds];
  if (!potd) bomb("problem allocating <potd>");

  dpot  = new Matrix [nthrds];
  if (!dpot) bomb("problem allocating <dpot>");

  for (int i=0; i<nthrds; i++) {
    potd[i].setsize(0,Lmax,1,nmax);
    dpot[i].setsize(0,Lmax,1,nmax);
  }

				// Work vectors
  u = new Vector [nthrds];
  du = new Vector [nthrds];
  if (!u) bomb("problem allocating <u>");
  if (!du) bomb("problem allocating <du>");

  for (int i=0; i<nthrds; i++) {
    u[i].setsize(0,nmax);
    du[i].setsize(0,nmax);
  }

  for (int l=0; l<=Lmax; l++) {
    for (int n=1; n<=nmax; n++) {
      normM[l][n] = norm(n-1,l);
    }
  }
    
  // Sin, cos
  
  cosm = new Vector [nthrds];
  if (!cosm) bomb("problem allocating <cosm>");

  sinm = new Vector [nthrds];
  if (!sinm) bomb("problem allocating <sinm>");

  for (int i=0; i<nthrds; i++) {
    cosm[i].setsize(0,Lmax);
    sinm[i].setsize(0,Lmax);
  }

  if (!self_consistent || initializing) determine_coefficients();

}

void CBrockDisk::initialize(void)
{
  string val;

  if (get_value("rmax", val)) rmax = atof(val.c_str());
  if (get_value("scale", val)) scale = atof(val.c_str());
  if (get_value("Lmax", val)) Lmax = atoi(val.c_str());
  if (get_value("nmax", val)) nmax = atoi(val.c_str());
  if (get_value("self_consistent", val)) {
    if (atoi(val.c_str())) self_consistent = true; 
    else self_consistent = false;
  }
  if (get_value("selector", val)) {
    if (atoi(val.c_str())) selector = true; 
    else selector = false;
  }
}

CBrockDisk::~CBrockDisk(void)
{
  delete [] expcoef0;
  if (selector) {
    delete [] cc;
    delete [] cc1;
    pthread_mutex_destroy(&cc_lock);
  }
  delete [] potd;
  delete [] dpot;
  delete [] u;
  delete [] du;
  delete [] cosm;
  delete [] sinm;
}

void CBrockDisk::get_acceleration_and_potential(Component* curComp)
{
  cC = curComp;

  /*========================================*/
  /* No coefficients for external particles */
  /*========================================*/

  if (use_external) {

    MPL_start_timer();
    determine_acceleration_and_potential();
    MPL_stop_timer();

    use_external = false;

    return;
  }

  /*======================*/
  /* Compute coefficients */
  /*======================*/

  if (self_consistent || initializing) determine_coefficients();

  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  determine_acceleration_and_potential();

  MPL_stop_timer();

  // Clear external potential flag
  use_external = false;
}

void CBrockDisk::determine_coefficients(void)
{
  int compute;
  int l, i, n, nn;

  if (selector) compute = !(this_step%npca);

				// Clean
  for (n=1; n<=nmax; n++) {
    for (l=0; l<=2*Lmax; l++) {
      expcoef[l][n] = 0.0;
      expcoef1[l][n] = 0.0;
      for (i=0; i<nthrds; i++) expcoef0[i][l][n] = 0.0;
      if (selector && compute) {
	for (nn=n; nn<=nmax; nn++) cc1[l][n][nn] = 0.0;
      }
    }
  }

  use0 = 0;
  use1 = 0;

  exp_thread_fork(true);

  for (i=0; i<nthrds; i++) {
    use1 += use[i];
    for (l=0; l<=Lmax*(Lmax+2); l++)
      for (n=1; n<=nmax; n++)
	expcoef1[l][n] += expcoef0[i][l][n];
  }

  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  used = use0;

  if (!selector) {

    MPI_Allreduce ( &expcoef1[0][1],
		    &expcoef[0][1],
		    nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	
    for (l=1;l<=Lmax; l++) {

      MPI_Allreduce ( &expcoef1[2*l - 1][1],
		      &expcoef[2*l - 1][1],
		      nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce ( &expcoef1[2*l    ][1],
		      &expcoef[2*l    ][1],
		      nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
  }

  if (selector) {

    parallel_gather_coefficients();

    if (myid == 0) pca_hall(compute);

    parallel_distribute_coefficients();

  }

}

void * CBrockDisk::determine_coefficients_thread(void * arg)
{
  double pos[3], xx, yy, zz, r, r2, phi, rs, fac1, fac2, mass;

  int nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = component->Adiabatic();

  map<unsigned long, Particle>::iterator it=cC->Particles().begin();
  unsigned long j;

  for (int i=0; i<nbeg; i++) it++;
  for (int i=nbeg; i<nend; i++) {

    j = (it++)->first;

    if (component->freeze(j)) // frozen particles do not respond
      continue;

    mass = cC->Mass(j)* adb;

    for (int k=0; k<3; k++) 
      pos[k] = cC->Pos(j, k) - component->center[k];

    xx = pos[0];
    yy = pos[1];
    zz = pos[2];

    r2 = (xx*xx + yy*yy);
    r = sqrt(r2) + DSMALL;

    if (r<=rmax) {
      use1++;
      phi = atan2(yy,xx);
      rs = r/scale;
	
      
      sinecosine_R(Lmax, phi, cosm[id], sinm[id]);
      get_potl(Lmax, nmax, rs, potd[id]);
      
				// l loop 

      for (int n=1; n<=nmax; n++) {
	expcoef0[id][0][n] += potd[id][0][n]*mass/normM[0][n];
	if (selector && compute) {
	  pthread_mutex_lock(&cc_lock);
	  for (int nn=n; nn<=nmax; nn++)
	    cc1[0][n][nn] += potd[id][0][n]*potd[id][0][nn]*mass/
	      (normM[0][n]*normM[0][nn]);
	  pthread_mutex_unlock(&cc_lock);
	}
      }
	
      for (int l=1;l<=Lmax; l++) {

	fac1 = cosm[id][l];
	fac2 = sinm[id][l];
	
	for (int n=1; n<=nmax; n++) {
	  expcoef0[id][2*l - 1][n] +=  potd[id][l][n]*fac1*mass/normM[l][n];
	  expcoef0[id][2*l    ][n] +=  potd[id][l][n]*fac2*mass/normM[l][n];
	  if (selector && compute) {
	    pthread_mutex_lock(&cc_lock);
	    for (int nn=n; nn<=nmax; nn++) {
	      cc1[2*l - 1][n][nn] += potd[id][l][n]*potd[id][l][nn]*fac1*fac1*
		mass/(normM[l][n]*normM[l][nn]);
	      cc1[2*l    ][n][nn] += potd[id][l][n]*potd[id][l][nn]*fac2*fac2*
		mass/(normM[l][n]*normM[l][nn]);
	    }
	    pthread_mutex_unlock(&cc_lock);
	  }
	}
      }
    }
  }

}



void CBrockDisk::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);
}

void * CBrockDisk::determine_acceleration_and_potential_thread(void * arg)
{
  int l;
  double r,rs,fac,phi,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps;

  double pos[3];
  double xx, yy, zz;

  int nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  map<unsigned long, Particle>::iterator it=cC->Particles().begin();
  unsigned long j;

  for (int i=0; i<nbeg; i++) it++;
  for (int i=nbeg; i<nend; i++) {

    j = it->first;

    for (int k=0; k<3; k++) 
      pos[k] = cC->Pos(j, k) - component->center[k];

    xx = pos[0];
    yy = pos[1];
    zz = pos[2];

    r = sqrt(xx*xx + yy*yy) + DSMALL;
    rs = r/scale;
    phi = atan2(yy,xx);

    sinecosine_R(Lmax, phi, cosm[id], sinm[id]);
    get_dpotl(Lmax, nmax, rs, potd[id], dpot[id]);
    get_pot_coefs_safe(0, expcoef[0], &p, &dp, potd[id], dpot[id]);

    potl = p;
    potr = dp;
    pott = potp = 0.0;
      
				// l loop
    
    for (l=1; l<=Lmax; l++) {

      get_pot_coefs_safe(l, expcoef[2*l - 1], &pc, &dpc, potd[id], dpot[id]);
      get_pot_coefs_safe(l, expcoef[2*l    ], &ps, &dps, potd[id], dpot[id]);
      potl += pc*cosm[id][l] + ps*sinm[id][l];
      potr += dpc*cosm[id][l] + dps*sinm[id][l];
      potp += (-pc*sinm[id][l] + ps*cosm[id][l])*l;
    }

    fac = xx*xx + yy*yy;
    
    potr /= scale*scale;
    potl /= scale;
    potp /= scale;

    cC->AddAcc(j, 0, -potr*xx/r);
    cC->AddAcc(j, 1, -potr*yy/r);
    cC->AddAcc(j, 2, -potr*zz/r);
    if (fac > DSMALL) {
      cC->AddAcc(j, 0,  potp*yy/fac);
      cC->AddAcc(j, 1, -potp*xx/fac);
    }

    if (use_external)
      cC->AddPotExt(j, potl);
    else
      cC->AddPot(j, potl);

    it++;
  }

}

void 
CBrockDisk::determine_fields_at_point_sph(double r, double theta, double phi,
					  double *tdens0, double *tpotl0, 
					  double *tdens, double *tpotl, 
					  double *tpotr, double *tpott, 
					  double *tpotp)

{
  determine_fields_at_point_polar(r, phi, tdens0, tpotl0, tdens, tpotl, tpotr, tpotp);
  *tpott = 0.0;
}


void 
CBrockDisk::determine_fields_at_point_cyl(double r, double z, double phi,
					  double *tdens0, double *tpotl0, 
					  double *tdens, double *tpotl, 
					  double *tpotr, double *tpott, 
					  double *tpotp)

{
  determine_fields_at_point_polar(r, phi, tdens0, tpotl0, tdens, tpotl, tpotr, tpotp);
  *tpott = 0.0;
}


void CBrockDisk::determine_fields_at_point_polar
(
 double r, double phi,
 double *tdens0, double *tpotl0,
 double *tdens, double *tpotl, double *tpotr, double *tpotp
 )
{
  int l;
  double rs,dp;
  double potr,potl,potp,p,pc,dpc,ps,dps,dens;

  rs = r/scale;

  sinecosine_R(Lmax, phi, cosm[0], sinm[0]);

  get_dens(Lmax, nmax, rs, dend);
  get_dpotl(Lmax, nmax, rs, potd[0], dpot[0]);

  get_dens_coefs(0,expcoef[0],&dens);

  get_pot_coefs(0,expcoef[0],&p,&dp);
  potl = p;
  potr = dp;
  potp = 0.0;
  
  *tdens0 = dens;
  *tpotl0 = potl;
      
  /*		l loop */
    
  for (l=1; l<=Lmax; l++) {
    
    get_dens_coefs(l,expcoef[2*l - 1],&pc);
    get_dens_coefs(l,expcoef[2*l    ],&ps);
    dens += pc*cosm[0][l] + ps*sinm[0][l];
    
    get_pot_coefs(l,expcoef[2*l - 1],&pc,&dpc);
    get_pot_coefs(l,expcoef[2*l    ],&ps,&dps);
    potl += pc*cosm[0][l] + ps*sinm[0][l];
    potr += dpc*cosm[0][l] + dps*sinm[0][l];
    potp += (-pc*sinm[0][l] + ps*cosm[0][l])*l;
  }

  *tdens0 /= scale*scale*scale;
  *tpotl0 /= scale;
      
  *tdens = dens/(scale*scale*scale);
  *tpotl = potl/scale;
  *tpotr = potr/(scale*scale);
  *tpotp = potp/scale;
}


void CBrockDisk::get_pot_coefs(int l, Vector& coef, double *p, double *dp)
{
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (i=1; i<=nmax; i++) {
    pp  += potd[0][l][i] * coef[i];
    dpp += dpot[0][l][i] * coef[i];
  }

  *p = -pp;
  *dp = -dpp;
}


void CBrockDisk::get_pot_coefs_safe(int l, Vector& coef, double *p, double *dp,
				    Matrix& potd1, Matrix& dpot1)
{
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (i=1; i<=nmax; i++) {
    pp  += potd1[l][i] * coef[i];
    dpp += dpot1[l][i] * coef[i];
  }

  *p = -pp;
  *dp = -dpp;
}


void CBrockDisk::get_dens_coefs(int l, Vector& coef, double *p)
{
  double pp;
  int i;

  pp = 0.0;

  for (i=1; i<=nmax; i++)
    pp  += dend[l][i] * coef[i];

  *p = pp;
}


				// Get potential functions by recursion

void CBrockDisk::get_dpotl(int lmax, int nmax, double r, Matrix& p, Matrix& dp)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac), rcum = 1.0;
  int l, nn;

  for (l=0; l<=lmax+1; l++) {
    cur = cur0;

    p[l][1] = cur*rcum;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      p[l][nn+1] = cur*rcum;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcum *= r;
  }


  for (l=0; l<=lmax; l++) {

    dp[l][1] = p[l][1]*l/r - p[l+1][1];
    if (nmax<1) break;
    dp[l][2] = p[l][2]*l/r - (p[l+1][2] - 2*p[l+1][1]);
    if (nmax<2) break;

    for (nn=2; nn<nmax; nn++) {
      dp[l][nn+1] = p[l][nn+1]*l/r - 
	( p[l+1][nn+1] - 2*p[l+1][nn] + p[l+1][nn-1] );
    }
  }

}

void CBrockDisk::get_potl(int lmax, int nmax, double r, Matrix& p)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac), rcum = 1.0;
  int l, nn;

  for (l=0; l<=lmax+1; l++) {
    cur = cur0;

    work[l][1]  = cur;
    p[l][1] = cur*rcum;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      p[l][nn+1] = cur*rcum;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcum *= r;
  }

}


void CBrockDisk::get_dens(int lmax, int nmax, double r, Matrix& d)
{

  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac);
  int l, nn;
  double rcum = 0.5/M_PI;

  for (l=0; l<=lmax+1; l++) {
    cur = cur0;

    work[l][1] = cur;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      work[l][nn+1] = cur;
    }
    cur0 *= fac*(2*(l+1) - 1);
  }

  for (l=0; l<=lmax; l++) {
    d[l][1] = work[l+1][1]*rcum;
    d[l][2] = work[l+1][2]*rcum;
    for (nn=2; nn<nmax; nn++)
      d[l][nn+1] = (work[l+1][nn+1] - work[l+1][nn-1])*rcum;

    rcum *= r;
  }

}

void CBrockDisk::get_potl_dens(int lmax, int nmax, double r, 
			       Matrix& p, Matrix& d)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac);
  int l, nn;
  double  rcump = 1.0, rcumd = 0.5/M_PI;

  for (l=0; l<=lmax+1; l++) {
    cur = cur0;

    work[l][1] = cur;
    curl1 = 0.0;

    for (nn=1; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/nn)*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/nn)*curl2;
      p[l][nn+1] = cur*rcump;
      work[l][nn+1] = cur;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcump *= r;
  }

  for (l=0; l<=lmax; l++) {
    d[l][1] = work[l+1][1]*rcumd;
    d[l][2] = work[l+1][2]*rcumd;
    for (nn=2; nn<nmax; nn++)
      d[l][nn+1] = (work[l+1][nn+1] - work[l+1][nn-1])*rcumd;

    rcumd *= r;
  }

}

double CBrockDisk::norm(int n, int m)
{
  double ans = 1.0;
  int i;
 
  for (i=n+1; i<=n+2*m; i++)
    ans *= i;

  return pow(0.5, 2*m+1)*ans;
}



				// Dump coefficients to a file

void CBrockDisk::dump_coefs(ostream& out)
{
  ostringstream sout;
  sout << id;

  char buf[64];
  for (int i=0; i<64; i++) {
    if (i<sout.str().length()) 
      buf[i] = sout.str().c_str()[i];
    else 
      buf[i] = ' ';
  }
  
  out.write((char *)&buf, 64*sizeof(char));

  out.write((char *)&tnow, sizeof(double));
  out.write((char *)&scale, sizeof(double));
  out.write((char *)&nmax, sizeof(int));
  out.write((char *)&Lmax, sizeof(int));

  for (int ir=1; ir<=nmax; ir++) {
    for (int l=0; l<=2*Lmax; l++)
      out.write((char *)&expcoef[l][ir], sizeof(double));
  }

}

