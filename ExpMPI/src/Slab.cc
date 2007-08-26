#include <stdlib.h>
#include <math.h>

#include "expand.h"

#include <kevin_complex.h>
#include <biorth1d.h>

#include <Slab.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

static const double KEPS=1.0e-6;

Slab::Slab(string& line) : PotAccel(line)
{
  id = "Slab (trigonometric)";
  nminx = nminy = 0;
  nmaxx = nmaxy = nmaxz = 10;
  zmax = 10.0;
  coef_dump = true;

  initialize();

  imx = 1+2*nmaxx;
  imy = 1+2*nmaxy;
  imz = nmaxz;
  jmax = imx*imy*imz;

  expccof = new Complex* [nthrds];
  for (int i=0; i<nthrds; i++)
    expccof[i] = new Complex [jmax];

  expreal  = new double [jmax];
  expreal1 = new double [jmax];
  expimag  = new double [jmax];
  expimag1 = new double [jmax];
  
  dfac = 2.0*M_PI;
  kfac = Complex(0.0, dfac);
  
  nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

				// I believe that these are thread safe
  trig = new OneDTrig* [nnmax+1];
  for (int i=0; i<=nnmax; i++) {
    trig[i] = new OneDTrig [i+1];
    for (int j=0; j<=i; j++) 
      trig[i][j].reset(dfac*sqrt((double)(i*i + j*j))+KEPS, zmax);
  }

  zpot = new Vector [nthrds];
  zfrc = new Vector [nthrds];
  for (int i=0; i<nthrds; i++) {
    zpot[i].setsize(1, nmaxz);
    zfrc[i].setsize(1, nmaxz);
  }

}

Slab::~Slab()
{
  for (int i=0; i<nthrds; i++) delete [] expccof[i];
  delete [] expccof;

  delete [] expreal;
  delete [] expreal1;
  delete [] expimag;
  delete [] expimag1;
  
  for (int i=0; i<=nnmax; i++) delete [] trig[i];
  delete [] trig;

  delete [] zpot;
  delete [] zfrc;

}

void Slab::initialize()
{
  string val;

  if (get_value("nmaxx", val)) nmaxx = atoi(val.c_str());
  if (get_value("nmaxy", val)) nmaxy = atoi(val.c_str());
  if (get_value("nmaxz", val)) nmaxz = atoi(val.c_str());
  if (get_value("nminx", val)) nminx = atoi(val.c_str());
  if (get_value("nminy", val)) nminy = atoi(val.c_str());
  if (get_value("zmax",  val)) zmax  = atof(val.c_str());

}

void Slab::determine_coefficients(void)
{
				//  Coefficients are ordered as follows:
				//  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax
				//  in a single array for each dimension
				//  with z dimension changing most rapidly

  // Clean 

  for (int indx=0; indx<jmax; indx++) {
    expreal[indx] = 0.0;
    expimag[indx] = 0.0;
    expreal1[indx] = 0.0;
    expimag1[indx] = 0.0;
  }

  for (int i=0; i<nthrds; i++) {
    use[i] = 0;
    for (int indx=0; indx<jmax; indx++) expccof[i][indx] = 0.0;
  }

  exp_thread_fork(true);

  int used1 = 0;
  used = 0;
  for (int i=0; i<nthrds; i++) used1 += use[i];
  
  MPI_Allreduce ( &used1, &used,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for (int i=0; i<nthrds; i++) {
    for (int indx=0; indx<jmax; indx++) {
      expreal1[indx] += expccof[i][indx].real();
      expimag1[indx] += expccof[i][indx].imag();
    }
  }
  
  MPI_Allreduce( expreal1, expreal, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce( expimag1, expimag, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int indx=0; indx<jmax; indx++)
    expccof[0][indx] = Complex(expreal[indx], expimag[indx]);
}

void * Slab::determine_coefficients_thread(void * arg)
{
  int ix, iy, iz, iix, iiy, indx, dum;

  Complex startx, starty, facx, facy;
  Complex stepx, stepy;

  unsigned  nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = component->Adiabatic();
  double zz;

  for (int i=nbeg; i<nend; i++) {
    
				// Increment particle counter
    use[id]++;

				// Truncate to box with sides in [0,1]
    
    if (cC->Pos(i, 0)<0.0)
      cC->AddPos(i, 0, (double)((int)fabs(cC->Pos(i, 0)) + 1.0) );
    else
      cC->AddPos(i, 0, -(double)((int)cC->Pos(i, 0)));
    
    if (cC->Pos(i, 1)<0.0)
      cC->AddPos(i, 1, (double)((int)fabs(cC->Pos(i, 1))) + 1.0);
    else
      cC->AddPos(i, 1, -(double)((int)cC->Pos(i, 1)));
    

				// Recursion multipliers
    stepx = exp(-kfac*cC->Pos(i, 0));
    stepy = exp(-kfac*cC->Pos(i, 1));
    
				// Initial values
    startx = exp(nmaxx*kfac*cC->Pos(i, 0));
    starty = exp(nmaxy*kfac*cC->Pos(i, 1));
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
      iix = nmaxx - ix;
      if (iix<0) iix = ix - nmaxx;
      
      if (iix < 0 || iix > nmaxx) {
	cerr << "Out of bounds: iix=" << iix << endl;
      }

      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	iiy = nmaxy - iy;
	if (iiy<0) iiy = iy - nmaxy;
	
	if (iiy < 0 || iiy > nmaxy) {
	  cerr << "Out of bounds: iiy=" << iiy << endl;
	}
	
	zz = cC->Pos(i, 2) - cC->center[2];

	if (iix>=iiy)
	  trig[iix][iiy].potl(dum, dum, zz, zpot[id]);
	else
	  trig[iiy][iix].potl(dum, dum, zz, zpot[id]);


	for (iz=0; iz<imz; iz++) {

	  indx = imz*(iy + imy*ix) + iz;

                             // |--- density in orthogonal series
                             // |    is 4.0*M_PI rho
                             // v
	  expccof[id][indx] += -4.0*M_PI*cC->Mass(i)*adb*
	    facx*facy*zpot[id][iz+1];
	}
      }
    }
  }
    
}

void Slab::get_acceleration_and_potential(Component* C)
{
  cC = C;

  determine_coefficients();

  MPL_start_timer();
  exp_thread_fork(false);
  MPL_stop_timer();
}

void * Slab::determine_acceleration_and_potential_thread(void * arg)
{
  int ix, iy, iz, iix, iiy, ii, jj, indx, dum;

  Complex fac, startx, starty, facx, facy, potl, facf;
  Complex stepx, stepy;
  Complex accx, accy, accz;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double zz;

  for (int i=nbeg; i<nend; i++) {
    
    accx = accy = accz = potl = 0.0;
    
				// Recursion multipliers
    stepx = exp(kfac*cC->Pos(i, 0));
    stepy = exp(kfac*cC->Pos(i, 1));

				// Initial values (note sign change)
    startx = exp(-nmaxx*kfac*cC->Pos(i, 0));
    starty = exp(-nmaxy*kfac*cC->Pos(i, 1));
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
				// Compute wavenumber; recall that the
				// coefficients are stored as follows:
				// -nmax,-nmax+1,...,0,...,nmax-1,nmax
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
	
	zz = cC->Pos(i, 2) - cC->center[2];

	if (iix>=iiy) {
	  trig[iix][iiy].potl(dum, dum, zz, zpot[id]);
	  trig[iix][iiy].force(dum, dum, zz, zfrc[id]);
	}
	else {
	  trig[iiy][iix].potl(dum, dum, zz, zpot[id]);
	  trig[iiy][iix].force(dum, dum, zz, zfrc[id]);
	}

	
	for (iz=0; iz<imz; iz++) {
	  
	  indx = imz*(iy + imy*ix) + iz;
	  
	  fac  = facx*facy*zpot[id][iz+1]*expccof[0][indx];
	  facf = facx*facy*zfrc[id][iz+1]*expccof[0][indx];
	  
	  
				// Don't add potential for 
				// zero wavenumber (constant term)
	  
	  if (ii==0 && jj==0 && iz==0) fac = facf = 0.0;
	  
				// Limit to minimum wave number
	  
	  if (abs(ii)<nminx || abs(jj)<nminy) continue;
	  
	  potl += fac;
	  
	  accx += -dfac*ii*Complex(0.0,1.0)*fac;
	  accy += -dfac*jj*Complex(0.0,1.0)*fac;
	  accz += facf;
	  
	}
      }
    }
    
    cC->AddAcc(i, 0, Re(accx));
    cC->AddAcc(i, 1, Re(accy));
    cC->AddAcc(i, 2, Re(accz));
    if (use_external)
      cC->AddPotExt(i, Re(potl));
    else
      cC->AddPot(i, Re(potl));
  }
}

void Slab::dump_coefs(ostream& out)
{
  coefheader.time = tnow;
  coefheader.zmax = zmax;
  coefheader.h = 0.0;
  coefheader.type = ID;
  coefheader.nmaxx = nmaxx;
  coefheader.nmaxy = nmaxy;
  coefheader.nmaxz = nmaxz;
  coefheader.jmax = (1+2*nmaxx)*(1+2*nmaxy)*nmaxz;
  
  out.write((char *)&coefheader, sizeof(SlabCoefHeader));
  for (int i=0; i<coefheader.jmax; i++) {
    out.write((char *)&expccof[0][i].real(), sizeof(double));
    out.write((char *)&expccof[0][i].imag(), sizeof(double));
  }
}

