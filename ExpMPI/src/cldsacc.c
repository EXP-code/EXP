/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the effect of cloud-cloud interactions for biorth
 *  expansion code
 *
 *  Adjusts cloud size to preserve constant density
 *
 *  Routines:
 *  --------
 *  void new_cld(m,i,j)           creates a new cloud of mass m, index j
 *                                deriving the mass from cloud i.  The
 *                                velocity is computed from that of cloud i
 *                                plus a dispersion
 *
 *  void new_stars(m,i,j)         same as above but stars are created with
 *                                identical velocities to the parent cloud
 * 
 *  void cloud_interaction(i,j)   collide clouds i and j
 *
 *  void cloud_collision_complex(i,j,t,d)
 *  void cloud_collision(i,j,t)  --"test version"
 *                                compute the results of a cloud collision
 *                                between clouds i and j at time t with 
 *                                distance of closest approach d.  In the
 *                                "complex" version, clouds may stick,
 *                                or intersect and produce stars.
 *
 *                                In the shorter "test version" clouds stick
 *                                and form stars
 *
 *  void prune_bodies(void)       utility: remove bodies flagged as inactive
 *                                         using the eliminate array
 *
 *  void swap_bodies(i,j)         utility: swap phase-space for indices i and j
 *
 *  Returns:
 *  -------
 *  As given
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *
 ***************************************************************************/

#include "expand.h"

static char rcsid[] = "$Id$";

void new_cld(double,int,int);
void new_stars(double,int,int);
void cloud_interaction(int,int);
void cloud_collision(int,int,double);
void cloud_collision_complex(int,int,double,double);
void prune_bodies(void);
void swap_bodies(int,int);
void return_stats(double *, double *, double *, double *, double *, double *,
		  double *, double *, double *, double *);

int *eliminate,*previous,nbodies0;

double vrelcld=0.0;		/* For accumulating vrel statistics */
double vrelcld2=0.0;
int nvrel=0;

				/* This routine is called by the main
				   evolution routine */
void clouds(void)
{
  static int firstime=1;
  static int *rindx;
  static double *Rmin,*Rmax;
  double rnew,tnew,delmax,maxsize=0.0;
  int i,j,ii,jj;

  nbodies0=nbodies;

				/* Set up utility stuff on first pass */
  if (firstime) {
    firstime = 0;

				/* Set up index for sorting */
    rindx = ivector(1,nbodmax);
    Rmin = dvector(1,nbodmax);
    Rmax = dvector(1,nbodmax);
    eliminate = ivector(1,nbodmax);
    previous = ivector(1,nbodmax);
    for (i=1; i<=nbodmax; i++) previous[i] = 0;

				/* Random generator setup */
    rnd_init(seed);
  }
    
				/* Assumes radii have already been computed
				   in self-grav case */
  
				/* Sort by radius */
  indexx(nbodies,rr,rindx);

				/* Find minimum and maximum radii during
				   previous time step using rectilinear
				   orbit approximation */

  for (i=1; i<=nbodies; i++) {

    rnew = sqrt(
		(x[i] - vx[i]*dtime)*(x[i] - vx[i]*dtime) + 
		(y[i] - vy[i]*dtime)*(y[i] - vy[i]*dtime) + 
		(z[i] - vz[i]*dtime)*(z[i] - vz[i]*dtime) );
    if (rnew > rr[i]) {
      Rmin[i] = rr[i];
      Rmax[i] = rnew;
    } else {
      Rmin[i] = rnew;
      Rmax[i] = rr[i];
    }

    tnew = (x[i]*vx[i] + y[i]*vy[i] + z[i]*vz[i])/
      (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    if (tnew >=0.0 && tnew <=dtime)
      Rmin[i] = sqrt(
		(x[i] - vx[i]*tnew)*(x[i] - vx[i]*tnew) + 
		(y[i] - vy[i]*tnew)*(y[i] - vy[i]*tnew) + 
		(z[i] - vz[i]*tnew)*(z[i] - vz[i]*tnew) );

    maxsize = MAX(maxsize, size[i]);
    eliminate[i] = 0;
  }


				/* Find candidate interacting neighbors 
				   by looking for overlaps in radial
				   positions; the routine cloud_interaction()
				   will determine whether or not there has
				   actually been a collision */

				/* NB: rr[rindx[k]] is rr in ascending order
				       for k=1, 2, 3, . . . */

  for (i=1; i<=nbodies0; i++) {
    ii = rindx[i];
    if (mass[ii]<=0.0 || size[ii]<=0.0) continue;
    for (j=i+1; j<=nbodies0; j++) {
      jj = rindx[j];

				/* Objects must be within contact range */

      if ( (fabs(Rmax[jj]-Rmin[ii]) > 2.0*maxsize) &&
	   (fabs(Rmax[ii]-Rmin[jj]) > 2.0*maxsize) ) break;

				/* Objects must be clouds and not have
				   interacted with the same cloud last */

      if (mass[jj]>0.0 && size[jj]>0.0 && !eliminate[jj] &&
	  previous[ii] != jj && previous[jj] != ii ) cloud_interaction(ii,jj);

				/* Quit if fiducial object has ceased to be
				   a cloud */

      if (mass[ii]<=0.0 || size[ii]<=0.0) break;
    }
  }

  prune_bodies();

}


void cloud_interaction(int i, int j)
{
  double tclose,dx,dy,dz,dmin;

				/* Compute time of closest approach */
  tclose = (
	    (x[i]-x[j])*(vx[i]-vx[j]) + 
	    (y[i]-y[j])*(vy[i]-vy[j]) + 
	    (z[i]-z[j])*(vz[i]-vz[j]) )/
	   (
	    (vx[i]-vx[j])*(vx[i]-vx[j]) + 
	    (vy[i]-vy[j])*(vy[i]-vy[j]) + 
	    (vz[i]-vz[j])*(vz[i]-vz[j]) );

				/* Must be in during timestep */
  if (tclose > dtime) tclose = dtime;
  if (tclose < 0.0)   tclose = 0.0;
  
				/* Compute distance of closest approach */
  dx = x[i] - x[j] + (vx[j]-vx[i])*tclose;
  dy = y[i] - y[j] + (vy[j]-vy[i])*tclose;
  dz = z[i] - z[j] + (vz[j]-vz[i])*tclose;

  dmin = sqrt(dx*dx + dy*dy + dz*dz);
				/* Is there physical contact? */
  if (dmin < size[i]+size[j] ) {
    if (complexcloud)
      cloud_collision_complex(i,j,tclose,dmin);
    else
      cloud_collision(i,j,tclose);
 
				/* Add clouds to last collision list
				   to prevent interaction on next time step */
    previous[i] = j;
    previous[j] = i;
 }
}


/* test version:
   
   assume relative COM motion dissipated
   so cloud has combined mass with total
   COM motion; e.g. "perfect sticking"
 */
  
void cloud_collision(int i, int j, double t)
{

				/* Assign COM motion at impact */

  size[i] *= pow((mass[i]+mass[j])/mass[i],0.33333333333333);
  vx[i] = (mass[i]*vx[i] + mass[j]*vx[j])/(mass[i] + mass[j]);
  vy[i] = (mass[i]*vy[i] + mass[j]*vy[j])/(mass[i] + mass[j]);
  vz[i] = (mass[i]*vz[i] + mass[j]*vz[j])/(mass[i] + mass[j]);
  if (metallicity)
    metal[i] = (mass[i]*metal[i] + mass[j]*metal[j])/(mass[i]+mass[j]);
  mass[i] += mass[j];
  mass[j] = 0.0;
  ninteract++;
				/* Make new stars */

  if (stareff>0.0 && nbodies<nbodmax) {
    nbodies++;
    mass[nbodies] = mass[i]*stareff;
    mass[i] *= (1.0-stareff);
    x[nbodies] = x[i];
    y[nbodies] = y[i];
    z[nbodies] = z[i];

    vx[nbodies] = vx[i];
    vy[nbodies] = vy[i];
    vz[nbodies] = vz[i];
    size[nbodies] = 0.0;
    pot[nbodies] = pot[i];
    potext[nbodies] = potext[i];
    if (metallicity) metal[nbodies] = metal[i];
  }

  if (metallicity) metal[i] += 1.0;

}




/*
  Revised version after discussion by Kwan and Weinberg 11/22/91
*/

void cloud_collision_complex(int i, int j, double t, double d)
{
  double vrel2,ke,pe,dmass,dmassi,dmassj,area;
  int icld;


  /* 
    Check to see if clouds are bound . . . if so:
    ---------------------------------------------
    1) no star formation 
    2) update size according to constant density prescription
  */ 

  ninteract++;

  vrel2 =  (vx[i]-vx[j])*(vx[i]-vx[j]) + 
           (vy[i]-vy[j])*(vy[i]-vy[j]) + 
           (vz[i]-vz[j])*(vz[i]-vz[j]);

  vrelcld += sqrt(vrel2);
  vrelcld2 += vrel2;
  nvrel++;

  ke = mass[i]*mass[j]/(mass[i]+mass[j]) * vrel2;

  pe = -mass[i]*mass[j]/(size[i] + size[j] + DSMALL) 
    -0.6*mass[i]*mass[i]/size[i] 
    -0.6*mass[j]*mass[j]/size[j];
			
  if (ke + pe/sqrt(gravscale) < 0.0) {
    				/* Assign COM motion at impact */
  
    nmerge++;

    size[i] *= pow((mass[i]+mass[j])/mass[i],0.33333333333333);
    vx[i] = (mass[i]*vx[i] + mass[j]*vx[j])/(mass[i] + mass[j]);
    vy[i] = (mass[i]*vy[i] + mass[j]*vy[j])/(mass[i] + mass[j]);
    vz[i] = (mass[i]*vz[i] + mass[j]*vz[j])/(mass[i] + mass[j]);
    if (metallicity)
      metal[i] = (mass[i]*metal[i] + mass[j]*metal[j])/(mass[i]+mass[j]);
    mass[i] += mass[j];
    mass[j] = 0.0;
    eliminate[j] = 1;

    return;
  }


  /* 
    If not bound:
    -------------
    1) Compute geometric overlap
    2) Divide into 3 clouds, the unoverlapped segments keep their
    original velocities and overlap has the COM velocity
    3) Stars form from the overlap and the given efficiency.
    4) The remainder of the gas breaks into 3 clouds with given dispersion
  */
    
  /* Compute overlap assuming clouds are normally oriented cubes */
  
				/* Can we make more clouds? */

				/* [Must allow for creation of possibly 4
				   more bodies] */
  if (nbodies+5 > nbodmax) return;

				/* [Clouds are assumed to look like cubes
				    which sides of length 2*size */

				/* Is one cube inside the other? */

  if (size[i] > d + size[j]) {
    dmass = size[j]*size[j]/(size[i]*size[i]) * mass[i];

    size[i] *= pow((mass[i]-dmass)/mass[i],0.33333333333333);
    size[j] *= pow((dmass+mass[j])/mass[j],0.33333333333333);
    vx[j] = (dmass*vx[i] + mass[j]*vx[j])/(dmass + mass[j]);
    vy[j] = (dmass*vy[i] + mass[j]*vy[j])/(dmass + mass[j]);
    vz[j] = (dmass*vz[i] + mass[j]*vz[j])/(dmass + mass[j]);
    if (metallicity)
      metal[j] = (dmass*metal[i] + mass[j]*metal[j])/(dmass+mass[j]);

    mass[i] -= dmass;
    mass[j] += dmass;

    icld = j;
  }
  else if (size[j] > d + size[i]) {
    dmass = size[i]*size[i]/(size[j]*size[j]) * mass[j];
    size[j] *= pow((mass[j]-dmass)/mass[j],0.33333333333333);
    size[i] *= pow((dmass+mass[i])/mass[i],0.33333333333333);
    vx[i] = (dmass*vx[j] + mass[i]*vx[i])/(dmass + mass[i]);
    vy[i] = (dmass*vy[j] + mass[i]*vy[i])/(dmass + mass[i]);
    vz[i] = (dmass*vz[j] + mass[i]*vz[i])/(dmass + mass[i]);
    if (metallicity)
      metal[i] = (mass[i]*metal[i] + dmass*metal[j])/(mass[i]+dmass);


    mass[j] -= dmass;
    mass[i] += dmass;

    swap_bodies(i,j);
    icld = j;
  }				/* At least part of each cloud survies */
  else {
    area = (size[i] + size[j] - d)*2.0*MIN(size[i],size[j]);
    dmassi = area/(4.0*size[i]*size[i]) * mass[i];
    dmassj = area/(4.0*size[j]*size[j]) * mass[j];
    dmass = dmassi + dmassj;

    size[i] *= pow((mass[i]-dmassi)/mass[i],0.33333333333333);
    size[j] *= pow((mass[j]-dmassj)/mass[j],0.33333333333333);

				/* Make a new cloud out of overlap */
    nbodies++;
    
    size[nbodies] = size[i] * pow(dmass/mass[i],0.33333333333333);
    x[nbodies] = x[i];
    y[nbodies] = y[i];
    z[nbodies] = z[i];
    vx[nbodies] = (dmassi*vx[i] + dmassj*vx[j])/(dmassi + dmassj);
    vy[nbodies] = (dmassi*vy[i] + dmassj*vy[j])/(dmassi + dmassj);
    vz[nbodies] = (dmassi*vz[i] + dmassj*vz[j])/(dmassi + dmassj);
    mass[nbodies] = dmass;
    pot[nbodies] = pot[i];
    potext[nbodies] = potext[i];
    if (metallicity) metal[nbodies] = 
      (dmassi*metal[i] + dmassj*metal[j])/(dmassi + dmassj);

    mass[i] -= dmassi;
    mass[j] -= dmassj;

    icld = nbodies;
  }

				/* Make new stars and clouds out of overlap */

  if (stareff>0.0) {
    dmass = mass[icld]*(1.0-stareff)/3.0;
    new_cld(dmass,icld,++nbodies);
    if (metallicity) metal[nbodies] += 1.0;
    new_cld(dmass,icld,++nbodies);
    if (metallicity) metal[nbodies] += 1.0;
    new_cld(dmass,icld,++nbodies);
    if (metallicity) metal[nbodies] += 1.0;

    size[icld] = 0.0;		/* Remainder becomes stars */
  }
}



				/* Make new stars */
				/* Take mass away from object i and turn
				   it into stars in object j */

void new_stars(double dmass, int i, int j)
{
  mass[j] = dmass;
  mass[i] -= dmass;
  x[j] = x[i];
  y[j] = y[i];
  z[j] = z[i];

  vx[j] = vx[i];
  vy[j] = vy[i];
  vz[j] = vz[i];
  size[j] = 0.0;
  pot[j] = pot[i];
  potext[j] = potext[i];
  if (metallicity) metal[j] = metal[i];
}

/* Make new cloud with gaussian vel.disp. */
/* Take mass away from object i and turn
   it into clouds in object j */

double gasdev(void);		/* Returns a gaussian distributed random
				 variable with zero mean and unit std dev */

void new_cld(double dmass, int i, int j)
{
  double rj;

  size[j] = size[i] * pow(dmass/mass[i],0.33333333333333);
  size[i] *= pow((mass[i]-dmass)/mass[i],0.33333333333333);
  mass[j] = dmass;
  mass[i] -= dmass;
  x[j] = x[i];
  y[j] = y[i];
  z[j] = z[i];
  rj = sqrt(x[j]*x[j] + y[j]*y[j] + z[j]*z[j]);

  vx[j] = vx[i] + veld*gasdev();
  vy[j] = vy[i] + veld*gasdev();
  vz[j] = vz[i] + veld*gasdev();
  pot[j] = pot[i];
  potext[j] = potext[i];
  if (metallicity) metal[j] = metal[i];
}


void prune_bodies(void)
{
  int i,j;

  for (i=nbodies; i>0; i--) {
    if (eliminate[i]) {
      for (j=i+1; j<=nbodies; j++) {

	mass[j-1] = mass[j];

	x[j-1] = x[j];
	y[j-1] = y[j];
	z[j-1] = z[j];

	vx[j-1] = vx[j];
	vy[j-1] = vy[j];
	vz[j-1] = vz[j];
	size[j-1] = size[j];
	pot[j-1] = pot[j];
	potext[j-1] = potext[j];
	if (metallicity) metal[j-1] = metal[j];
      }
      nbodies--;
    }
  }
  
}

void swap_bodies(int i, int j)
{
  double mt,xt,yt,zt,vxt,vyt,vzt,st,pt,pet,met;

  mt = mass[i];
  xt = x[i];
  yt = y[i];
  zt = z[i];
  vxt = vx[i];
  vyt = vy[i];
  vzt = vz[i];
  st = size[i];
  pt = pot[i];
  pet = potext[i];
  if (metallicity) met = metal[i];

  mass[i] = mass[j];
  x[i] = x[j];
  y[i] = y[j];
  z[i] = z[j];
  vx[i] = vx[j];
  vy[i] = vy[j];
  vz[i] = vz[j];
  size[i] = size[j];
  pot[i] = pot[j];
  potext[i] = potext[j];
  if (metallicity) metal[i] = metal[j];

  mass[j] = mt;
  x[j] = xt;
  y[j] = yt;
  z[j] = zt;
  vx[j] = vxt;
  vy[j] = vyt;
  vz[j] = vzt;
  size[j] = st;
  pot[j] = pt;
  potext[j] = pet;
  if (metallicity) metal[j] = met;

}


void return_stats(double *vrel, double *vrel2, double *r2stars, double *r2clds, double *v2stars, double *v2clds, double *masstar, double *mcld, double *mcld1, double *mcld2)
{
  int i;

  *r2stars = *r2clds = *v2stars = *v2clds = *masstar = *mcld = *mcld1 = *mcld2 = 0.0;

  for (i=1; i<=nbodies; i++) {
    if (size[i]<=0.0) {
      *masstar += mass[i];

      *r2stars += mass[i]*(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
      *v2stars += mass[i]*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    }
    else {
      *mcld += mass[i];

      *r2clds += mass[i]*(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
      *v2clds += mass[i]*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

      if (size[i]<cloudsize) *mcld1 += mass[i];
      if (size[i]>cloudsize) *mcld2 += mass[i];
    }
  }

  if (*masstar>0.0) {
    *r2stars /= *masstar;
    *v2stars /= *masstar;
  }

  if (*mcld>0.0) {
    *r2clds /= *mcld;
    *v2clds /= *mcld;
  }

  if (nvrel) {
    *vrel = vrelcld/nvrel;
    *vrel2 = vrelcld2/nvrel;

    nvrel = 0;
    vrelcld = vrelcld2 = 0.0;
  }

  
}
