#include "expand.h"
#include <strstream>
#include <SphericalBasis.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

SphericalBasis::SphericalBasis(string& line) : AxisymmetricBasis(line)
{
  dof = 3;
  geometry = sphere;
  coef_dump = true;

  string val;

  if (get_value("scale", val)) 
    scale = atof(val.c_str());
  else
    scale = 1.0;

  if (get_value("rmax", val)) 
    rmax = atof(val.c_str());
  else
    rmax = 10.0;

  if (get_value("self_consistent", val)) {
    if (atoi(val.c_str())) self_consistent = true; 
    else self_consistent = false;
  } else
    self_consistent = true;

  if (get_value("selector", val)) {
    if (atoi(val.c_str())) selector = true; 
    else selector = false;
  } else
    selector = false;

  if (get_value("NO_L1", val)) {
    if (atoi(val.c_str())) NO_L1 = true; 
    else NO_L1 = false;
  } else
    NO_L1 = false;

  Lmax = Lmax<1 ? 1 : Lmax;

  if (nthrds<1) nthrds=1;


  initialize();


  //	Allocate coefficient matrix
  expcoef.setsize(0, Lmax*(Lmax+2), 1, nmax);
  expcoef1.setsize(0, Lmax*(Lmax+2), 1, nmax);
  
  expcoef0 = new Matrix [nthrds];
  if (!expcoef0) bomb("problem allocating <expcoef0>");

  for (int i=0; i<nthrds; i++)
    expcoef0[i].setsize(0, Lmax*(Lmax+2), 1, nmax);

  // Allocate normalization matrix

  normM.setsize(0,Lmax,1,nmax);
  for (int l=0; l<=Lmax; l++) {
    for (int n=1; n<=nmax; n++) {
      normM[l][n] = 1.0;
    }
  }

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

  // Potential and deriv matrices

  normM.setsize(0,Lmax,1,nmax);
  krnl.setsize(0,Lmax,1,nmax);
  dend.setsize(0,Lmax,1,nmax);

  potd  = new Matrix [nthrds];
  if (!potd) bomb("problem allocating <potd>");

  dpot  = new Matrix [nthrds];
  if (!dpot) bomb("problem allocating <dpot>");

  for (int i=0; i<nthrds; i++) {
    potd[i].setsize(0, Lmax, 1, nmax);
    dpot[i].setsize(0, Lmax, 1, nmax);
  }

  // Sin, cos, legendre

  cosm = new Vector [nthrds];
  if (!cosm) bomb("problem allocating <cosm>");

  sinm = new Vector [nthrds];
  if (!sinm) bomb("problem allocating <sinm>");

  legs = new Matrix [nthrds];
  if (!legs) bomb("problem allocating <legs>");

  dlegs = new Matrix [nthrds];
  if (!dlegs) bomb("problem allocating <dlegs>");

  for (int i=0; i<nthrds; i++) {
    cosm[i].setsize(0,Lmax);
    sinm[i].setsize(0,Lmax);
    legs[i].setsize(0,Lmax,0,Lmax);
    dlegs[i].setsize(0,Lmax,0,Lmax);
  }

				/* Work vectors */
  u = new Vector [nthrds];
  du = new Vector [nthrds];
  if (!u) bomb("problem allocating <u>");
  if (!du) bomb("problem allocating <du>");

  for (int i=0; i<nthrds; i++) {
    u[i].setsize(0,nmax);
    du[i].setsize(0,nmax);
  }

  // Factorial matrix

  factorial.setsize(0, Lmax, 0, Lmax);

  for (int l=0; l<=Lmax; l++) {
    for (int m=0; m<=l; m++) 
      factorial[l][m] = factrl(l-m)/factrl(l+m);
  }

  // Per thread counter
  use = new int [nthrds];
  if (!use) bomb("problem allocating <use>");


  firstime_coef  = true;
  firstime_accel = true;

}

void SphericalBasis::setup(void)
{				// Call normalization and kernel
  for (int l=0; l<=Lmax; l++) {	// with current binding from derived class
    for (int n=1; n<=nmax; n++) {
      normM[l][n] = norm(n-1,l);
      krnl[l][n] = knl(n-1,l);
    }
  }
}  


SphericalBasis::~SphericalBasis()
{
  delete [] expcoef0;
  if (selector) {
    delete [] cc;
    delete [] cc1;
    pthread_mutex_destroy(&cc_lock);
  }
  delete [] potd;
  delete [] dpot;
  delete [] cosm;
  delete [] sinm;
  delete [] legs;
  delete [] dlegs;
  delete [] u;
  delete [] du;
  delete [] use;
}

void SphericalBasis::initialize()
{
				// Do nothing
}

void SphericalBasis::check_range()
{
				// Do nothing
}

void SphericalBasis::get_acceleration_and_potential(vector<Particle>* P)
{
				
  particles = P;		// "Register" particles
  nbodies = particles->size();	// And compute number of bodies

  /*====================================================*/
  /* Accel & pot using previously computed coefficients */
  /*====================================================*/

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

  if (firstime_accel || self_consistent) {
    firstime_accel = false;
    determine_coefficients();
  }


  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  determine_acceleration_and_potential();

  MPL_stop_timer();

  // Clear external potential flag
  use_external = false;
}


void * SphericalBasis::determine_coefficients_thread(void * arg)
{
  int l, i, loffset, moffset, m, n, nn;
  double r, r2, rs, fac1, fac2, costh, phi, mass;
  double facs1=0.0, facs2=0.0, fac0=4.0*M_PI;
  double xx, yy, zz;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  use[id] = 0;

  for (i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;

    mass = (*particles)[i].mass;

    xx = (*particles)[i].pos[0] - component->center[0];
    yy = (*particles)[i].pos[1] - component->center[1];
    zz = (*particles)[i].pos[2] - component->center[2];

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


void SphericalBasis::determine_coefficients(void)
{
  int l, i, loffset, moffset, m, n, nn, use0, use1;

#ifdef MPE_PROFILE
  MPE_Log_event(9, myid, "b_compute_coef");
#endif

  if (selector) compute = !(this_step%npca) || firstime_coef;


  /*		Clean */
  for (n=1; n<=nmax; n++) {
      for (l=0; l<=Lmax*(Lmax+2); l++) {
	expcoef[l][n]  = 0.0;
	expcoef1[l][n] = 0.0;
	for (i=0; i<nthrds; i++) expcoef0[i][l][n] = 0.0;
	if (selector && compute) 
	  for (nn=n; nn<=nmax; nn++) cc1[l][n][nn] = 0.0;
      }
    }

  use0 = 0;
  use1 = 0;
  used = 0;

  exp_thread_fork(true);

  for (i=0; i<nthrds; i++) {
    use1 += use[i];
    for (l=0; l<=Lmax*(Lmax+2); l++)
      for (n=1; n<=nmax; n++)
	expcoef1[l][n] += expcoef0[i][l][n];
  }
    
  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  used = use0;

#ifdef MPE_PROFILE
  MPE_Log_event(10, myid, "e_compute_coef");
#endif

  if (!selector) {

#ifdef MPE_PROFILE
    MPE_Log_event(7, myid, "b_distrib_c");
#endif

    for (l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
	/*		m loop */
	for (m=0, moffset=0; m<=l; m++) {
	  if (m==0) {

	    MPI_Allreduce ( &expcoef1[loffset+moffset][1],
			    &expcoef[loffset+moffset][1],
			    nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	    moffset++;
	  }
	  else {

	    MPI_Allreduce ( &expcoef1[loffset+moffset][1],
			    &expcoef[loffset+moffset][1],
			    nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	    MPI_Allreduce ( &expcoef1[loffset+moffset+1][1],
			    &expcoef[loffset+moffset+1][1],
			    nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    moffset+=2;
	  }
	}
      }

#ifdef MPE_PROFILE
    MPE_Log_event(8, myid, "e_distrib_c");
#endif

  }
 
  if (selector) {
    
    parallel_gather_coefficients();

    if (myid == 0) pca_hall(compute);

    parallel_distribute_coefficients();
    
    firstime_coef = 0;
  }

}

void * SphericalBasis::determine_acceleration_and_potential_thread(void * arg)
{
  int l,loffset,moffset,m,ioff;
  double r,rs,r0=0.0,fac,fac1,fac2,fac3,fac4,costh,phi,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps,facp,facdp;
  double dfac=0.25/M_PI;

  double pos[3];
  double xx, yy, zz;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;

    fac1 = dfac;

    for (int k=0; k<3; k++) 
      pos[k] = (*particles)[i].pos[k] - component->center[k];

    xx = pos[0];
    yy = pos[1];
    zz = pos[2];

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
      
    /*		l loop */
    
    for (l=1, loffset=1; l<=Lmax; loffset+=(2*l+1), l++) {

      /*		m loop */
      for (m=0, moffset=0; m<=l; m++) {
	fac1 = (2.0*l+1.0)/(4.0*M_PI);
	if (m==0) {
				/* Suppress L=1 terms? */
	  if ( !NO_L1 || !(l==1) ) {
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
	  }
	  moffset++;
	}
	else {
				/* Suppress L=1 terms? */
	  if ( !NO_L1 || !(l==1) ) {
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
	  }
	  moffset +=2;
	}
      }
    }

    fac = xx*xx + yy*yy;

    potr /= scale*scale;
    potl /= scale;
    pott /= scale;
    potp /= scale;

    (*particles)[i].acc[0] += -(potr*xx/r - pott*xx*zz/(r*r*r) );
    (*particles)[i].acc[1] += -(potr*yy/r - pott*yy*zz/(r*r*r) );
    (*particles)[i].acc[2] += -(potr*zz/r + pott*fac/(r*r*r));
    if (fac > DSMALL2) {
      (*particles)[i].acc[0] +=  potp*yy/fac;
      (*particles)[i].acc[1] += -potp*xx/fac;
    }
    if (use_external)
      (*particles)[i].potext += potl;
    else
      (*particles)[i].pot += potl;

  }

  return (NULL);
}


void SphericalBasis::determine_acceleration_and_potential(void)
{
#ifdef MPE_PROFILE
  MPE_Log_event(11, myid, "b_compute_force");
#endif

  exp_thread_fork(false);

#ifdef MPE_PROFILE
  MPE_Log_event(12, myid, "e_compute_force");
#endif

}


void SphericalBasis::get_pot_coefs(int l, Vector& coef, double *p, double *dp)
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


void SphericalBasis::get_pot_coefs_safe(int l, Vector& coef, double *p, double *dp,
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


void SphericalBasis::get_dens_coefs(int l, Vector& coef, double *p)
{
  double pp;
  int i;

  pp = 0.0;

  for (i=1; i<=nmax; i++)
    pp  += dend[l][i] * coef[i];

  *p = pp;
}
				// Dump coefficients to a file

void SphericalBasis::dump_coefs(ostream& out)
{
  char buf[64];
  ostrstream sout(buf, 64);
  sout << setfill(' ') << id << '\0';

  out.write(&buf, 64*sizeof(char));
  out.write(&tnow, sizeof(double));
  out.write(&scale, sizeof(double));
  out.write(&nmax, sizeof(int));
  out.write(&Lmax, sizeof(int));

  for (int ir=1; ir<=nmax; ir++) {
    for (int l=0; l<=Lmax*(Lmax+2); l++)
      out.write(&expcoef[l][ir], sizeof(double));
  }

}


void SphericalBasis::determine_fields_at_point_cyl(double R, double z, double phi, double *tdens, double *tpotl, double *tpotR, double *tpotz, double *tpotp)
{
  double r = sqrt(R*R + z*z) + 1.0e-18;
  double theta = acos(z/R);
  double tpotr, tpott;

  determine_fields_at_point_sph(r, theta, phi, tdens, tpotl, &tpotr,
				&tpott, tpotp);
  *tpotR = tpotr*sin(theta) - tpott*cos(theta);
  *tpotz = tpotr*cos(theta) + tpott*sin(theta);
}

void SphericalBasis::determine_fields_at_point_sph(double r, double theta, double phi, double *tdens, double *tpotl, double *tpotr, double *tpott, double *tpotp)
{
  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double potr,potl,pott,potp,p,pc,dpc,ps,dps,dens;
  double dfac=0.25/M_PI;

  rs = r/scale;
  costh = cos(theta);

  fac1 = dfac;

  dlegendre_R(Lmax, costh, legs[0], dlegs[0]);
  sinecosine_R(Lmax, phi, cosm[0], sinm[0]);
  get_dens(Lmax, nmax, rs, dend, 0);
  get_dpotl(Lmax, nmax, rs, potd[0], dpot[0], 0);
  get_dens_coefs(0,expcoef[0],&dens);
  dens *= dfac*dfac;

  get_pot_coefs(0,expcoef[0],&p,&dp);
  potl = fac1*p;
  potr = fac1*dp;
  pott = potp = 0.0;
  

  // l loop
    
  for (l=1, loffset=1; l<=Lmax; loffset+=(2*l+1), l++) {
    
    // m loop
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
	fac2 = fac1*legs[0][l][m];
	get_dens_coefs(l,expcoef[loffset+moffset],&p);
	dens += dfac*fac2*p;
	get_pot_coefs(l,expcoef[loffset+moffset],&p,&dp);
	potl += fac2*p;
	potr += fac2*dp;
	pott += fac1*dlegs[0][l][m]*p;
	moffset++;
      }
      else {
	fac2 = 2.0 * fac1 * factorial[l][m];
	fac3 = fac2 * legs[0][l][m];
	fac4 = fac2 * dlegs[0][l][m];
	
	get_dens_coefs(l,expcoef[loffset+moffset],&pc);
	get_dens_coefs(l,expcoef[loffset+moffset+1],&ps);
	dens += dfac*fac3*(pc*cosm[0][m] + ps*sinm[0][m]);
	
	get_pot_coefs(l,expcoef[loffset+moffset],&pc,&dpc);
	get_pot_coefs(l,expcoef[loffset+moffset+1],&ps,&dps);
	potl += fac3*(pc*cosm[0][m] + ps*sinm[0][m]);
	potr += fac3*(dpc*cosm[0][m] + dps*sinm[0][m]);
	pott += fac4*(pc*cosm[0][m] + ps*sinm[0][m]);
	potp += fac3*(-pc*sinm[0][m] + ps*cosm[0][m])*m;
	moffset +=2;
      }
    }
  }

  *tdens = dens/(scale*scale*scale);
  *tpotl = potl/scale;
  *tpotr = potr/(scale*scale);
  *tpott = pott/scale;
  *tpotp = potp/scale;
  
}
