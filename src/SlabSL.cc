#include <stdlib.h>
#include <math.h>

#include "expand.h"

#include <kevin_complex.h>
#include <biorth1d.h>

#include <SlabSL.H>

SlabSL::SlabSL(const YAML::Node& conf) : PotAccel(conf)
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

  expccof = new KComplex* [nthrds];
  for (int i=0; i<nthrds; i++)
    expccof[i] = new KComplex[jmax];
  
  expreal  = new double [jmax];
  expreal1 = new double [jmax];
  expimag  = new double [jmax];
  expimag1 = new double [jmax];
    
  dfac = 2.0*M_PI;
  kfac = KComplex(0.0, dfac);
    
  nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

  zpot = new Vector [nthrds];
  zfrc = new Vector [nthrds];
  for (int i=0; i<nthrds; i++) {
    zpot[i].setsize(1, nmaxz);
    zfrc[i].setsize(1, nmaxz);
  }

}

SlabSL::~SlabSL()
{
  for (int i=0; i<nthrds; i++) delete [] expccof[i];
  delete [] expccof;

  delete [] expreal;
  delete [] expreal1;
  delete [] expimag;
  delete [] expimag1;
  
  delete grid;

  delete [] zpot;
  delete [] zfrc;

}

void SlabSL::initialize()
{
  try {
    if (conf["nmaxx"])          nmaxx       = conf["nmaxx"].as<int>();
    if (conf["nmaxy"])          nmaxy       = conf["nmaxy"].as<int>();
    if (conf["nmaxz"])          nmaxz       = conf["nmaxz"].as<int>();
    if (conf["nminx"])          nminx       = conf["nminx"].as<int>();
    if (conf["nminy"])          nminy       = conf["nminy"].as<int>();
    if (conf["hslab"])          hslab       = conf["hslab"].as<double>();
    if (conf["zmax"])           zmax        = conf["zmax"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in SlabSL: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}

void SlabSL::determine_coefficients(void)
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
    expccof[0][indx] = KComplex(expreal[indx], expimag[indx]);
}

void * SlabSL::determine_coefficients_thread(void * arg)
{
  int ix, iy, iz, iix, iiy, ii, jj, indx;

  KComplex startx, starty, facx, facy;
  KComplex stepx, stepy;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = cC->Adiabatic();
  double zz;

  for (int i=nbeg; i<nend; i++) {
    
				// Increment particle counter
    use[id]++;

				// Truncate to box with sides in [0,1]
    
    if (cC->Pos(i, 0)<0.0)
      cC->AddPos(i, 0, (double)((int)fabs(cC->Pos(i, 0))) + 1.0 );
    else
      cC->AddPos(i, 0, -(double)((int)cC->Pos(i, 0)) );
    
    if (cC->Pos(i, 1)<0.0)
      cC->AddPos(i, 1, (double)((int)fabs(cC->Pos(i, 1))) + 1.0 );
    else
      cC->AddPos(i, 1, -(double)((int)cC->Pos(i, 1)) );
    

				// Recursion multipliers
    stepx = exp(-kfac*cC->Pos(i, 0));
    stepy = exp(-kfac*cC->Pos(i, 1));
   
				// Initial values
    startx = exp(nmaxx*kfac*cC->Pos(i, 0));
    starty = exp(nmaxy*kfac*cC->Pos(i, 1));
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
      ii = ix - nmaxx;
      iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	jj = iy - nmaxy;
	iiy = abs(jj);
	
	if (iix > nmaxx) {
	  cerr << "Out of bounds: iix=" << ii << endl;
	}
	if (iiy > nmaxy) {
	  cerr << "Out of bounds: iiy=" << jj << endl;
	}
	
	zz = cC->Pos(i, 2, Component::Centered);

	if (iix>=iiy)
	  grid->get_pot(zpot[id], zz, iix, iiy);
	else
	  grid->get_pot(zpot[id], zz, iiy, iix);


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
    
  return (NULL);
}

void SlabSL::get_acceleration_and_potential(Component* C)
{
  cC = C;

  MPL_start_timer();
  exp_thread_fork(false);
  MPL_stop_timer();
}

void * SlabSL::determine_acceleration_and_potential_thread(void * arg)
{
  int ix, iy, iz, iix, iiy, ii, jj, indx;
  KComplex fac, startx, starty, facx, facy, potl, facf;
  KComplex stepx, stepy;
  KComplex accx, accy, accz;

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
	
	if (iix > nmaxx) {
	  cerr << "Out of bounds: ii=" << ii << endl;
	}
	if (iiy > nmaxy) {
	  cerr << "Out of bounds: jj=" << jj << endl;
	}
	
	zz = cC->Pos(i, 2, Component::Centered);

	if (iix>=iiy) {
	  grid->get_pot  (zpot[id], zz, iix, iiy);
	  grid->get_force(zfrc[id], zz, iix, iiy);
	}
	else {
	  grid->get_pot  (zpot[id], zz, iiy, iix);
	  grid->get_force(zfrc[id], zz, iiy, iix);
	}

	
	for (iz=0; iz<imz; iz++) {
	  
	  indx = imz*(iy + imy*ix) + iz;
	  
	  fac  = facx*facy*zpot[id][iz+1]*expccof[0][indx];
	  facf = facx*facy*zfrc[id][iz+1]*expccof[0][indx];
	  
				// Limit to minimum wave number
	  
	  if (abs(ii)<nminx || abs(jj)<nminy) continue;
	  
	  potl += fac;
	  
	  accx += -dfac*ii*KComplex(0.0,1.0)*fac;
	  accy += -dfac*jj*KComplex(0.0,1.0)*fac;
	  accz += -facf;
	  
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

  return (NULL);
}

void SlabSL::dump_coefs(ostream& out)
{
  coefheader.time = tnow;
  coefheader.zmax = zmax;
  coefheader.h = hslab;
  coefheader.type = ID;
  coefheader.nmaxx = nmaxx;
  coefheader.nmaxy = nmaxy;
  coefheader.nmaxz = nmaxz;
  coefheader.jmax = (1+2*nmaxx)*(1+2*nmaxy)*nmaxz;
  
  out.write((char *)&coefheader, sizeof(SlabSLCoefHeader));
  for (int i=0; i<coefheader.jmax; i++) {
    out.write((char *)&expccof[0][i].real(), sizeof(double));
    out.write((char *)&expccof[0][i].imag(), sizeof(double));
  }
}

