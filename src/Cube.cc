#include "expand.H"

#include <stdlib.h>
#include <math.h>

#include <Cube.H>

const std::set<std::string> Cube::valid_keys = {
  "nminx",
  "nminy",
  "nminz",
  "nmaxx",
  "nmaxy",
  "nmaxz"
};


Cube::Cube(Component* c0, const YAML::Node& conf) : PotAccel(c0, conf)
{
  id = "Cube";

  nminx = nminy = nminz = 0;
  nmaxx = nmaxy = nmaxz = 16;

  initialize();

  imx  = 1+2*nmaxx;
  imy  = 1+2*nmaxy;
  imz  = 1+2*nmaxz;
  jmax = imx*imy*imz;

  expccof.resize(nthrds);
  for (auto & v : expccof) v.resize(jmax);

  dfac = 2.0*M_PI;
  kfac = std::complex<double>(0.0,dfac);
}

Cube::~Cube(void)
{
}

void Cube::initialize(void)
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("Cube", "parameter", unmatched, __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["nminx"]) nminx = conf["nminx"].as<int>();
    if (conf["nminy"]) nminy = conf["nminy"].as<int>();
    if (conf["nminz"]) nminz = conf["nminz"].as<int>();
    if (conf["nmaxx"]) nmaxx = conf["nmaxx"].as<int>();
    if (conf["nmaxy"]) nmaxy = conf["nmaxy"].as<int>();
    if (conf["nmaxz"]) nmaxz = conf["nmaxz"].as<int>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Cube: "
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

void * Cube::determine_coefficients_thread(void * arg)
{

  int ix,iy,iz,indx;
  std::complex<double> startx,starty,startz,facx,facy,facz;
  std::complex<double> stepx,stepy,stepz;
  double mass;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = component->Adiabatic();

  use[id] = 0;

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {

    i = it->first;
    it++;

    use[id]++;
    mass = cC->Mass(i) * adb;

				/* Truncate to cube with sides in [0,1] */
    
    if (cC->Pos(i, 0)<0.0)
      cC->AddPos(i, 0, (double)((int)fabs(cC->Pos(i, 0))) + 1.0);
    else
      cC->AddPos(i, 0, -(double)((int)cC->Pos(i, 0)));
    
    if (cC->Pos(i, 1)<0.0)
      cC->AddPos(i, 1, (double)((int)fabs(cC->Pos(i, 1))) + 1.0);
    else
      cC->AddPos(i, 1, -(double)((int)cC->Pos(i, 1)));
    
    if (cC->Pos(i, 2)<0.0)
      cC->AddPos(i, 2, (double)((int)fabs(cC->Pos(i, 2))) + 1.0);
    else
      cC->AddPos(i, 2, -(double)((int)cC->Pos(i, 2)));
    
    
				/* Recursion multipliers */
    stepx = exp(-kfac*cC->Pos(i, 0));
    stepy = exp(-kfac*cC->Pos(i, 1));
    stepz = exp(-kfac*cC->Pos(i, 2));
    
				/* Initial values */
    startx = exp(static_cast<double>(nmaxx)*kfac*cC->Pos(i, 0));
    starty = exp(static_cast<double>(nmaxy)*kfac*cC->Pos(i, 1));
    startz = exp(static_cast<double>(nmaxz)*kfac*cC->Pos(i, 2));
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	for (facz=startz, iz=0; iz<imz; iz++, facz*=stepz) {
	  
	  indx = imz*(iy + imy*ix) + iz;
	  expccof[id][indx] += mass*facx*facy*facz;
	  
	}
      }
    }
  }
    
  return (NULL);
}

void Cube::determine_coefficients(void)
{
				//  Coefficients are ordered as follows:
				//  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax
				//  in a single array for each dimension
				//  with z dimension changing most rapidly
  // Clean 
  for (auto & v : expccof) v.setZero();

  exp_thread_fork(true);

  for (int i=0; i<nthrds; i++) use1 += use[i];
  
  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  used = use0;

  for (int i=1; i<nthrds; i++) {
    for (int indx=0; indx<jmax; indx++) {
      expccof[0][indx] += expccof[i][indx];
    }
  }

  MPI_Allreduce( MPI_IN_PLACE, expccof[0].data(), jmax,
		 MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
}

void * Cube::determine_acceleration_and_potential_thread(void * arg)
{
  int ix,iy,iz,ii,jj,kk,indx;
  std::complex<double> fac,startx,starty,startz,facx,facy,facz,dens,potl,accx,accy,accz;
  std::complex<double> stepx,stepy,stepz;
  double k2;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    
    i = it->first; it++;

    accx = accy = accz = dens = potl = 0.0;
    
				/* Recursion multipliers */
    stepx = exp(kfac*cC->Pos(i, 0));
    stepy = exp(kfac*cC->Pos(i, 1));
    stepz = exp(kfac*cC->Pos(i, 2));
    
				/* Initial values (note sign change) */
    startx = exp(-static_cast<double>(nmaxx)*kfac*cC->Pos(i, 0));
    starty = exp(-static_cast<double>(nmaxy)*kfac*cC->Pos(i, 1));
    startz = exp(-static_cast<double>(nmaxz)*kfac*cC->Pos(i, 2));
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	for (facz=startz, iz=0; iz<imz; iz++, facz*=stepz) {
	  
	  indx = imz*(iy + imy*ix) + iz;
	  
	  fac = facx*facy*facz*expccof[0][indx];
	  dens += fac;
	  
				/* Compute wavenumber; recall that the */
				/* coefficients are stored as follows: */
				/* -nmax,-nmax+1,...,0,...,nmax-1,nmax */
	  ii = ix-nmaxx;
	  jj = iy-nmaxy;
	  kk = iz-nmaxz;
	  
				/* No contribution to acceleration and */
	                        /* potential ("swindle") for zero      */
	                        /* wavenumber */
	  
	  if (ii==0 && jj==0 && kk==0) continue;
	  
				/* Limit to minimum wave number */
	  
	  if (abs(ii)<nminx || abs(jj)<nminy || abs(kk)<nminz) continue;
	  
	  k2 = 4.0*M_PI/(dfac*dfac*(ii*ii + jj*jj + kk*kk));
	  potl -= k2*fac;
	  
	  accx -= k2*std::complex<double>(0.0,-dfac*ii)*fac;
	  accy -= k2*std::complex<double>(0.0,-dfac*jj)*fac;
	  accz -= k2*std::complex<double>(0.0,-dfac*kk)*fac;
	  
	}
      }
    }
    
    cC->AddAcc(i, 0, accx.real());
    cC->AddAcc(i, 1, accy.real());
    cC->AddAcc(i, 2, accz.real());
    if (use_external)
      cC->AddPotExt(i, potl.real());
    else
      cC->AddPot(i, potl.real());
  }
  
  return (NULL);
}

void Cube::get_acceleration_and_potential(Component* C)
{
  cC = C;

  exp_thread_fork(false);

  // Clear external potential flag
  use_external = false;
}
