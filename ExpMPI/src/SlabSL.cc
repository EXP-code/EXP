
#include <stdlib.h>
#include <math.h>

#include "expand.h"

#include <kevin_complex.h>
#include <biorth1d.h>

#include <SlabSL.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

SlabSL::SlabSL(string& line) : PotAccel(line)
{
  id = "Slab (Sturm-Liouville)";
  NGRID = 100;
  nminx = nminy = 0;
  nmaxx = nmaxy = nmaxz = 10;
  zmax = 10.0;
  hslab = 0.2;
  coef_dump = true;

  initialize();

  SLGridSlab::mpi = 1;
  SLGridSlab::ZBEG = 0.0;
  SLGridSlab::ZEND = 0.1;
  SLGridSlab::H = hslab;
  
  int nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

  grid = new SLGridSlab(nnmax, nmaxz, NGRID, zmax);

  imx = 1+2*nmaxx;
  imy = 1+2*nmaxy;
  imz = nmaxz;
  jmax = imx*imy*imz;

  expccof = new Complex* [nthrds];
  for (int i=0; i<nthrds; i++)
    expccof[i] = new Complex[jmax];
  
  expreal = new double [jmax];
  expreal1 = new double [jmax];
  expimag = new double [jmax];
  expimag1 = new double [jmax];
    
  dfac = 2.0*M_PI;
  kfac = Complex(0.0, dfac);
    
  nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

  zpot = new Vector [nthrds];
  zfrc = new Vector [nthrds];
  for (int i=0; i<nthrds; i++) {
    zpot[i].setsize(1, nmaxz);
    zfrc[i].setsize(1, nmaxz);
  }
}

void SlabSL::initialize()
{
  string val;

  if (get_value("nmaxx", val)) nmaxx = atoi(val.c_str());
  if (get_value("nmaxy", val)) nmaxy = atoi(val.c_str());
  if (get_value("nmaxz", val)) nmaxz = atoi(val.c_str());
  if (get_value("nminx", val)) nminx = atoi(val.c_str());
  if (get_value("nminy", val)) nminy = atoi(val.c_str());
  if (get_value("hslab", val)) hslab = atof(val.c_str());
  if (get_value("zmax", val))  zmax = atof(val.c_str());
}

void SlabSL::determine_coefficients(void)
{
				//  Coefficients are ordered as follows:
				//  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax
				//  in a single array for each dimension
				//  with z dimension changing most rapidly

  // Clean 

  for (int indx=0; indx<jmax; indx++) {
    for (int i=0; i<nthrds; i++) expccof[i][indx] = 0.0;
    expreal[indx] = 0.0;
    expimag[indx] = 0.0;
    expreal1[indx] = 0.0;
    expimag1[indx] = 0.0;
  }

  exp_thread_fork(true);

  int use0, use1=0;
  for (int i=0; i<nthrds; i++) use1 += use[i];
  
  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  used = use0;

  for (int i=0; i<nthrds; i++) {
    for (int indx=0; indx<jmax; indx++) {
      expreal1[indx] = expccof[i][indx].real();
      expimag1[indx] = expccof[i][indx].imag();
    }
  }

  MPI_Allreduce( expreal1, expreal, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce( expimag1, expimag, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int indx=0; indx<jmax; indx++)
    expccof[0][indx] = Complex(expreal[indx], expimag[indx]);
}

void * SlabSL::determine_coefficients_thread(void * arg)
{
  int ix, iy, iz, iix, iiy, indx;

  Complex startx, starty, facx, facy;
  Complex stepx, stepy;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = component->Adiabatic();

  for (int i=nbeg; i<nend; i++) {
    

				// Truncate to box with sides in [0,1]
    
    if ((*particles)[i].pos[0]<0.0)
      (*particles)[i].pos[0] += 
	(double)((int)fabs((*particles)[i].pos[0])) + 1.0;
    else
      (*particles)[i].pos[0] -= (double)((int)(*particles)[i].pos[0]);
    
    if ((*particles)[i].pos[1]<0.0)
      (*particles)[i].pos[1] += 
	(double)((int)fabs((*particles)[i].pos[1])) + 1.0;
    else
      (*particles)[i].pos[1] -= (double)((int)(*particles)[i].pos[1]);
    

				/* Recursion multipliers */
    stepx = exp(-kfac*(*particles)[i].pos[0]);
    stepy = exp(-kfac*(*particles)[i].pos[1]);
   
				/* Initial values */
    startx = exp(nmaxx*kfac*(*particles)[i].pos[0]);
    starty = exp(nmaxy*kfac*(*particles)[i].pos[1]);
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
      iix = nmaxx - ix;
      if (iix<0) iix = ix - nmaxx;
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	iiy = nmaxy - iy;
	if (iiy<0) iiy = iy - nmaxy;
	
	if (iix < 0 || iix > nmaxx) {
	  cerr << "Out of bounds: iix=" << iix << endl;
	}
	if (iiy < 0 || iiy > nmaxy) {
	  cerr << "Out of bounds: iiy=" << iiy << endl;
	}
	
	if (iix>=iiy)
	  grid->get_pot(zpot[id], (*particles)[i].pos[2], iix, iiy);
	else
	  grid->get_pot(zpot[id], (*particles)[i].pos[2], iiy, iix);


	for (iz=0; iz<imz; iz++) {

	  indx = imz*(iy + imy*ix) + iz;
                             // |--- density in orthogonal series
                             // |    is 4.0*M_PI rho
			     // v
	  expccof[id][indx] += -4.0*M_PI*(*particles)[i].mass*adb*
	    facx*facy*zpot[id][iz+1];
	}
      }
    }
  }
    
}

void SlabSL::get_acceleration_and_potential(vector<Particle>* P)
{
  particles = P;

  determine_coefficients();

  exp_thread_fork(false);
}

void * SlabSL::determine_acceleration_and_potential_thread(void * arg)
{
  int ix, iy, iz, iix, iiy, ii, jj, indx;
  Complex fac, startx, starty, facx, facy, potl, facf;
  Complex stepx, stepy;
  Complex accx, accy, accz;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  for (int i=nbeg; i<nend; i++) {
    
    accx = accy = accz = potl = 0.0;
    
				/* Recursion multipliers */
    stepx = exp(kfac*(*particles)[i].pos[0]);
    stepy = exp(kfac*(*particles)[i].pos[1]);

				/* Initial values (note sign change) */
    startx = exp(-nmaxx*kfac*(*particles)[i].pos[0]);
    starty = exp(-nmaxy*kfac*(*particles)[i].pos[1]);
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
				/* Compute wavenumber; recall that the */
				/* coefficients are stored as follows: */
				/* -nmax,-nmax+1,...,0,...,nmax-1,nmax */
      ii = ix - nmaxx;
      iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	jj = iy - nmaxy;
	iiy = abs(jj);
	
	if (iix < 0 || iix > nmaxx) {
	  cerr << "Out of bounds: iix=" << iix << endl;
	}
	if (iiy < 0 || iiy > nmaxy) {
	  cerr << "Out of bounds: iiy=" << iiy << endl;
	}
	
	if (iix>=iiy) {
	  grid->get_pot  (zpot[id], (*particles)[i].pos[2], iix, iiy);
	  grid->get_force(zfrc[id], (*particles)[i].pos[2], iix, iiy);
	}
	else {
	  grid->get_pot  (zpot[id], (*particles)[i].pos[2], iiy, iix);
	  grid->get_force(zfrc[id], (*particles)[i].pos[2], iiy, iix);
	}
	
	for (iz=0; iz<imz; iz++) {
	  
	  indx = imz*(iy + imy*ix) + iz;
	  
	  fac  = facx*facy*zpot[id][iz+1]*expccof[0][indx];
	  facf = facx*facy*zfrc[id][iz+1]*expccof[0][indx];
	  
	  
				/* Don't add potential for 
				   zero wavenumber (constant term) */
	  
	  if (ii==0 && jj==0 && iz==0) fac = 0.0;
	  
	  
				/* Limit to minimum wave number */
	  
	  if (abs(ii)<nminx || abs(jj)<nminy) continue;
	  
	  potl -= fac;
	  
	  accx -= Complex(0.0,-dfac*ii)*fac;
	  accy -= Complex(0.0,-dfac*jj)*fac;
	  accz -= facf;
	  
	}
      }
    }
    
    (*particles)[i].acc[0] += Re(accx);
    (*particles)[i].acc[1] += Re(accy);
    (*particles)[i].acc[2] += Re(accz);
    if (use_external)
      (*particles)[i].potext += Re(potl);
    else
      (*particles)[i].pot += Re(potl);
  }
}

void SlabSL::dump_coefs(ostream& out)
{
  coefheader.time = tpos;
  coefheader.zmax = zmax;
  coefheader.nmaxx = nmaxx;
  coefheader.nmaxy = nmaxy;
  coefheader.nmaxz = nmaxz;
  coefheader.jmax = (1+2*nmaxx)*(1+2*nmaxy)*nmaxz;
  
  out.write(&coefheader, sizeof(SlabSLCoefHeader));
  for (int i=0; i<coefheader.jmax; i++) {
    out.write(&expccof[0][i].real(), sizeof(double));
    out.write(&expccof[0][i].imag(), sizeof(double));
  }
}

