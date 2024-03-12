#include <complex>
#include <cstdlib>
#include <cmath>

#include "expand.H"

#include <biorth1d.H>

#include <Slab.H>

static const double KEPS=1.0e-6;

const std::set<std::string>
Slab::valid_keys = {
  "nmaxx",
  "nmaxy",
  "nmaxz",
  "nminx",
  "nminy",
  "zmax"
};

Slab::Slab(Component* c0, const YAML::Node& conf) : PotAccel(c0, conf)
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

  expccof.resize(nthrds);
  for (auto & v : expccof) v.resize(jmax);
  
  dfac = 2.0*M_PI;
  kfac = std::complex<double>(0.0, dfac);
  
  nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

				// I believe that these are thread safe
  trig.resize(nnmax+1);
  for (int i=0; i<=nnmax; i++) {
    trig[i].resize(i+1);
    for (int j=0; j<=i; j++) 
      trig[i][j].reset(dfac*sqrt((double)(i*i + j*j))+KEPS, zmax);
  }

  zpot.resize(nthrds);
  zfrc.resize(nthrds);
  
  for (auto & v : zpot) v.resize(nmaxz);
  for (auto & v : zfrc) v.resize(nmaxz);

}

Slab::~Slab()
{
  // NADA
}

void Slab::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["nmaxx"])          nmaxx     = conf["nmaxx"].as<int>();
    if (conf["nmaxy"])          nmaxy     = conf["nmaxy"].as<int>();
    if (conf["nmaxz"])          nmaxz     = conf["nmaxz"].as<int>();
    if (conf["nminx"])          nminx     = conf["nminx"].as<int>();
    if (conf["nminy"])          nminy     = conf["nminy"].as<int>();
    if (conf["zmax"])           zmax      = conf["zmax"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Slab: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("Slab::initialize: error parsing YAML");
  }
}

void Slab::determine_coefficients(void)
{
				//  Coefficients are ordered as follows:
				//  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax
				//  in a single array for each dimension
				//  with z dimension changing most rapidly

  // Clean 

  for (int i=0; i<nthrds; i++) {
    use[i] = 0;
    expccof[i].setZero();
  }

  exp_thread_fork(true);

  int used1 = 0;
  used = 0;
  for (int i=0; i<nthrds; i++) used1 += use[i];
  for (int i=1; i<nthrds; i++) {
    used1 += use[i];
    for (int n=0; n<expccof[0].size(); n++)
      expccof[0][n] += expccof[i][n];
  }

  MPI_Allreduce ( &used1, &used,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   MPI_Allreduce( MPI_IN_PLACE, expccof[0].data(), expccof[0].size(),
		  MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
}

void * Slab::determine_coefficients_thread(void * arg)
{
  int ix, iy, iz, iix, iiy, indx, dum;

  std::complex<double> startx, starty, facx, facy;
  std::complex<double> stepx, stepy;

  unsigned  nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = component->Adiabatic();
  double zz;

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    
    i = (it++)->first;
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
    startx = exp(static_cast<double>(nmaxx)*kfac*cC->Pos(i, 0));
    starty = exp(static_cast<double>(nmaxy)*kfac*cC->Pos(i, 1));
    
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
    
  return (NULL);
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

  std::complex<double> fac, startx, starty, facx, facy, potl, facf;
  std::complex<double> stepx, stepy;
  std::complex<double> accx, accy, accz;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double zz;

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    
    i = (it++)->first;

    accx = accy = accz = potl = 0.0;
    
				// Recursion multipliers
    stepx = exp(kfac*cC->Pos(i, 0));
    stepy = exp(kfac*cC->Pos(i, 1));

				// Initial values (note sign change)
    startx = exp(-static_cast<double>(nmaxx)*kfac*cC->Pos(i, 0));
    starty = exp(-static_cast<double>(nmaxy)*kfac*cC->Pos(i, 1));
    
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
	  
	  accx += -dfac*ii*std::complex<double>(0.0,1.0)*fac;
	  accy += -dfac*jj*std::complex<double>(0.0,1.0)*fac;
	  accz += facf;
	  
	}
      }
    }
    
    cC->AddAcc(i, 0, accx.real());
    cC->AddAcc(i, 1, accy.real());
    cC->AddAcc(i, 2, accz.real());
    cC->AddPot(i, potl.real());
  }

  return (NULL);
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
  out.write((char *)expccof[0].data(), jmax*sizeof(std::complex<double>));
}

