#include <filesystem>
#include <cmath>

#include "expand.H"

#include <SlabSL.H>

const std::set<std::string>
SlabSL::valid_keys = {
  "nmaxx",
  "nmaxy",
  "nmaxz",
  "nminx",
  "nminy",
  "hslab",
  "zmax",
  "ngrid",
  "type"
};

//@{
//! These are for testing exclusively (should be set false for production)
static bool cudaAccumOverride = false;
static bool cudaAccelOverride = false;
//@}

SlabSL::SlabSL(Component* c0, const YAML::Node& conf) : PotAccel(c0, conf)
{
  id = "Slab (Sturm-Liouville)";
  nminx = nminy = 0;
  nmaxx = nmaxy = nmaxz = 6;
  zmax      = 10.0;
  hslab     = 0.2;
  coef_dump = true;

#if HAVE_LIBCUDA==1
  cuda_aware = true;
#endif

  initialize();

  SLGridSlab::mpi  = 1;
  SLGridSlab::ZBEG = 0.0;
  SLGridSlab::ZEND = 0.1;
  SLGridSlab::H    = hslab;
  
  int nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

  grid = std::make_shared<SLGridSlab>(nnmax, nmaxz, ngrid, zmax, type);

  // Test for basis consistency (will generate an exception if maximum
  // error is out of tolerance)
  //
  double worst = 0.0, diff = 0.0;
  int kxw = 0, kyw = 0;
  auto test = grid->orthoCheck(10000);
  for (int kx=0, indx=0; kx<=nnmax; kx++) {
    for (int ky=0; ky<=kx; ky++, indx++) {
      for (int n1=0; n1<nmaxz; n1++) {
	for (int n2=0; n2<nmaxz; n2++) {
	  if (n1==n2) diff = fabs(test[indx](n1, n2) - 1.0);
	  else diff = fabs(test[indx](n1, n2));
	  if (diff > worst) {
	    worst = diff;
	    kxw = kx;
	    kyw = ky;
	  }
	}
      }
    }
  }
	  
  if (true) {
    std::ofstream tmp("SlabSL.ortho");
    for (int kx=0, indx=0; kx<=nnmax; kx++) {
      for (int ky=0; ky<=kx; ky++, indx++) {
	tmp << "---- kx=" << kx << "  ky=" << ky << std::endl
	    << test[indx] << std::endl;
      }
    }
  }

  if (worst > __EXP__::orthoTol) {
    if (myid==0)
      std::cout << "SlabSL: orthogonality failure, worst=" << worst
		<< " at (" << kxw << ", " << kyw << ")" << std::endl;
    throw std::runtime_error("SlabSL: biorthogonal sanity check");
  } else {
    if (myid==0)
      std::cout << "---- SlabSL: biorthogonal check passed, worst="
		<< worst << std::endl;
  }

  imx  = 1+2*nmaxx;
  imy  = 1+2*nmaxy;
  imz  = nmaxz;
  jmax = imx*imy*imz;

  expccof.resize(nthrds);
  for (auto & v : expccof) v.resize(imx, imy, imz);
    
  dfac = 2.0*M_PI;
  kfac = std::complex<double>(0.0, dfac);
    
  zpot.resize(nthrds);
  zfrc.resize(nthrds);

  for (auto & v : zpot) v.resize(nmaxz);
  for (auto & v : zfrc) v.resize(nmaxz);

  // Allocate coefficient matrix (one for each multistep level)
  // and zero-out contents
  //
  differ1 = std::vector< std::vector<coefType> >(nthrds);
  for (int n=0; n<nthrds; n++) {
    differ1[n].resize(multistep+1);
    for (int i=0; i<=multistep; i++)
      differ1[n][i].resize(imx, imy, imz);
  }

  // MPI buffer space (including multistep arrays)
  //
  unsigned sz = (multistep+1)*jmax;
  pack  .resize(sz);
  unpack.resize(sz);

  // Coefficient evaluation times
  // 
  expcoefN.resize(multistep+1);
  expcoefL.resize(multistep+1);
  for (int i=0; i<=multistep; i++) {
    expcoefN[i] = std::make_shared<coefType>(imx, imy, imz);
    expcoefL[i] = std::make_shared<coefType>(imx, imy, imz);
    expcoefN[i] -> setZero();
    expcoefL[i] -> setZero();
  }
    
}

SlabSL::~SlabSL()
{
#if HAVE_LIBCUDA==1
  if (component->cudaDevice>=0) destroy_cuda();
#endif
}

void SlabSL::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["nmaxx"])          nmaxx       = conf["nmaxx"].as<int>();
    if (conf["nmaxy"])          nmaxy       = conf["nmaxy"].as<int>();
    if (conf["nmaxz"])          nmaxz       = conf["nmaxz"].as<int>();
    if (conf["nminx"])          nminx       = conf["nminx"].as<int>();
    if (conf["nminy"])          nminy       = conf["nminy"].as<int>();
    if (conf["ngrid"])          ngrid       = conf["ngrid"].as<int>();
    if (conf["hslab"])          hslab       = conf["hslab"].as<double>();
    if (conf["zmax" ])          zmax        = conf["zmax" ].as<double>();
    if (conf["type" ])          type        = conf["type" ].as<std::string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in SlabSL: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("SlabSL::initialze: error parsing YAML");
  }
}

void SlabSL::determine_coefficients(void)
{
  //  Coefficients are ordered as follows:
  //  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax in a single array for each
  //  dimension with z dimension changing most rapidly

  // Clean 

  for (int i=0; i<nthrds; i++) {
    use[i] = 0;
    expccof[i].setZero();
  }

  // Swap interpolation arrays
  //
  if (multistep) {

    auto p = expcoefL[mlevel];
  
    expcoefL[mlevel] = expcoefN[mlevel];
    expcoefN[mlevel] = p;
  
    // Clean arrays for current level
    //
    expcoefN[mlevel]->setZero();
  }

#if HAVE_LIBCUDA==1
  (*barrier)("SlabSL::entering cuda coefficients", __FILE__, __LINE__);
  if (component->cudaDevice>=0 and use_cuda) {
    if (cudaAccumOverride) {
      component->CudaToParticles();
      exp_thread_fork(true);
    } else {
      determine_coefficients_cuda();
      DtoH_coefs(mlevel);
    }
  } else {
    exp_thread_fork(true);
  }
  (*barrier)("SlabSL::exiting cuda coefficients", __FILE__, __LINE__);
#else
  exp_thread_fork(true);
#endif

  int used1 = 0, rank = expccof[0].size();
  used = 0;
  for (int i=1; i<nthrds; i++) {
    used1 += use[i];
    
    for (int j=0; j<rank; j++) expccof[0].data()[j] += expccof[i].data()[j];
  }
  
  MPI_Allreduce ( &used1, &used,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (multistep) {

    MPI_Allreduce( expccof[0].data(), expcoefN[mlevel]->data(),
		   expccof[0].size(),
		   MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  } else {
    
    MPI_Allreduce( MPI_IN_PLACE, expccof[0].data(), expccof[0].size(),
		   MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  }

  // Last level?
  //
  if (multistep and mlevel==multistep) {
     compute_multistep_coefficients();
  }

#if HAVE_LIBCUDA==1
  cuda_initialize();
#endif
}

void * SlabSL::determine_coefficients_thread(void * arg)
{
  int ix, iy;

  std::complex<double> startx, starty, facx, facy;
  std::complex<double> stepx, stepy;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = cC->Adiabatic();

  PartMapItr it = cC->Particles().begin();

  for (int q=0; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {

    unsigned long i = it->first;
    it++;
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
    startx = exp(static_cast<double>(nmaxx)*kfac*cC->Pos(i, 0));
    starty = exp(static_cast<double>(nmaxy)*kfac*cC->Pos(i, 1));
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
      int ii = ix - nmaxx;
      int iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	int jj = iy - nmaxy;
	int iiy = abs(jj);
	
	if (iix > nmaxx) {
	  std::cerr << "Out of bounds: iix=" << ii << std::endl;
	}
	if (iiy > nmaxy) {
	  std::cerr << "Out of bounds: iiy=" << jj << std::endl;
	}
	
	double zz = cC->Pos(i, 2, Component::Centered);

	if (iix>=iiy)
	  grid->get_pot(zpot[id], zz, iix, iiy);
	else
	  grid->get_pot(zpot[id], zz, iiy, iix);


	for (int iz=0; iz<imz; iz++) {
	                           // +--- density in orthogonal series
                                   // |    is 4.0*M_PI rho
                                   // v
	  expccof[id](ix, iy, iz) += -4.0*M_PI*cC->Mass(i)*adb*
	    facx*facy*zpot[id][iz];
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

#if HAVE_LIBCUDA==1
  if (use_cuda and cC->cudaDevice>=0 and cC->force->cudaAware()) {
    if (cudaAccelOverride) {
      cC->CudaToParticles();
      exp_thread_fork(false);
      cC->ParticlesToCuda();
    } else {
      // Copy coefficients from this component to device
      //
      HtoD_coefs();
      //
      // Do the force computation
      //
      determine_acceleration_cuda();
    }
  } else {

    exp_thread_fork(false);

  }
#else

  exp_thread_fork(false);

#endif

  MPL_stop_timer();
}

void * SlabSL::determine_acceleration_and_potential_thread(void * arg)
{
  int ix, iy;
  std::complex<double> fac, facx, facy, potl, facf;
  std::complex<double> accx, accy, accz;
  const std::complex<double> I(0.0, 1.0);

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  PartMapItr it = cC->Particles().begin();

  for (int q=0; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    
    unsigned long i = it->first; it++;

    accx = accy = accz = potl = 0.0;
    
				// Recursion multipliers
    std::complex<double> stepx = exp(kfac*cC->Pos(i, 0));
    std::complex<double> stepy = exp(kfac*cC->Pos(i, 1));

				// Initial values (note sign change)
    std::complex<double> startx = exp(-static_cast<double>(nmaxx)*kfac*cC->Pos(i, 0));
    std::complex<double> starty = exp(-static_cast<double>(nmaxy)*kfac*cC->Pos(i, 1));
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
				// Compute wavenumber; recall that the
				// coefficients are stored as follows:
				// -nmax,-nmax+1,...,0,...,nmax-1,nmax
      int ii = ix - nmaxx;
      int iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	int jj = iy - nmaxy;
	int iiy = abs(jj);
	
	if (iix > nmaxx) {
	  std::cerr << "Out of bounds: ii=" << ii << std::endl;
	}
	if (iiy > nmaxy) {
	  std::cerr << "Out of bounds: jj=" << jj << std::endl;
	}
	
	double zz = cC->Pos(i, 2, Component::Centered);

	if (iix>=iiy) {
	  grid->get_pot  (zpot[id], zz, iix, iiy);
	  grid->get_force(zfrc[id], zz, iix, iiy);
	}
	else {
	  grid->get_pot  (zpot[id], zz, iiy, iix);
	  grid->get_force(zfrc[id], zz, iiy, iix);
	}

	
	for (int iz=0; iz<imz; iz++) {
	  
	  fac  = facx*facy*zpot[id][iz]*expccof[0](ix, iy, iz);
	  facf = facx*facy*zfrc[id][iz]*expccof[0](ix, iy, iz);
	  
				// Limit to minimum wave number
	  
	  if (abs(ii)<nminx || abs(jj)<nminy) continue;
	  
	  potl += fac;
	  
	  accx += -kfac*static_cast<double>(ii)*fac;
	  accy += -kfac*static_cast<double>(jj)*fac;
	  accz += -facf;
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

void SlabSL::dump_coefs_h5(const std::string& file)
{
  // Add the current coefficients
  auto cur = std::make_shared<CoefClasses::SlabStruct>();

  cur->time     = tnow;
  cur->geom     = geoname[geometry];
  cur->id       = id;
  cur->time     = tnow;
  cur->nmaxx    = nmaxx;
  cur->nmaxy    = nmaxy;
  cur->nmaxz    = nmaxz;

  cur->allocate();		// Set the storage and copy the
				// coefficients through the map
  *cur->coefs   = expccof[0];

  // Check if file exists
  //
  if (std::filesystem::exists(file)) {
    slabCoefs.clear();
    slabCoefs.add(cur);
    slabCoefs.ExtendH5Coefs(file);
  }
  // Otherwise, extend the existing HDF5 file
  //
  else {
    // Copy the YAML config.  We only need this on the first call.
    std::ostringstream sout; sout << conf;
    cur->buf = sout.str();	// Copy to CoefStruct buffer

    // Add the name attribute.  We only need this on the first call.
    slabCoefs.setName(component->name);

    // And the new coefficients and write the new HDF5
    slabCoefs.clear();
    slabCoefs.add(cur);
    slabCoefs.WriteH5Coefs(file);
  }
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
  out.write((char *)expccof[0].data(),
	    expccof[0].size()*sizeof(std::complex<double>));
}

void SlabSL::multistep_update_begin()
{
  if (play_back and not play_cnew) return;

				// Clear the update matricies
  for (int n=0; n<nthrds; n++) {
    for (int M=mfirst[mdrft]; M<=multistep; M++) differ1[n][M].setZero();
  }
}

void SlabSL::multistep_update_finish()
{
  if (play_back and not play_cnew) return;

  // Combine the update matricies from all nodes
  //
  unsigned sz = (multistep - mfirst[mdrft] + 1)*jmax;

  // Zero the buffer space
  //
  std::fill(pack.begin(), pack.begin() + sz, 0.0);

  // Pack the difference matrices
  //
  for (int M=mfirst[mdrft]; M<=multistep; M++) {

    unsigned offset = (M - mfirst[mdrft])*jmax;

    for (int i=0; i<jmax; i++)
      for (int n=0; n<nthrds; n++) {
	pack[offset + i] += differ1[n][M].data()[i];
      }
  }

  MPI_Allreduce (pack.data(), unpack.data(), sz, 
		 MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  
  // Update the local coefficients
  //
  for (int M=mfirst[mdrft]; M<=multistep; M++) {

    unsigned offset = (M - mfirst[mdrft])*jmax;

    for (int i=0; i<jmax; i++)
      expcoefN[M]->data()[i] += unpack[offset+i];
  }
}

void SlabSL::multistep_update(int from, int to, Component *c, int i, int id)
{
  if (play_back and not play_cnew) return;
  if (c->freeze(i)) return;

  double mass = c->Mass(i) * component->Adiabatic();

  double x = c->Pos(i, 0);
  double y = c->Pos(i, 1);
  double z = c->Pos(i, 2);
  
  // Only compute for points inside the unit cube
  //
  if (x<0.0 or x>1.0) return;
  if (y<0.0 or y>1.0) return;
  if (z<0.0 or z>1.0) return;
    
  // Recursion multipliers
  std::complex<double> stepx = std::exp(-kfac*x);
  std::complex<double> stepy = std::exp(-kfac*y);
  std::complex<double> stepz = std::exp(-kfac*z);
    
  // Initial values for recursion
  std::complex<double> startx = std::exp(kfac*(x*nmaxx));
  std::complex<double> starty = std::exp(kfac*(y*nmaxy));
  std::complex<double> startz = std::exp(kfac*(z*nmaxz));
  
  std::complex<double> facx, facy, facz;
  int ix, iy, iz;

  for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
    for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
      for (facz=startz, iz=0; iz<imz; iz++, facz*=stepz) {
	
	// Compute wavenumber; recall that the coefficients are
	// stored as follows: -nmax,-nmax+1,...,0,...,nmax-1,nmax
	//
	int ii = ix-nmaxx;
	int jj = iy-nmaxy;
	int kk = iz-nmaxz;
	
	// Normalization
	double norm = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));;
	
	std::complex<double> val = -mass*facx*facy*facz*norm;

	differ1[id][from](ix, iy, iz) -= val;
	differ1[id][  to](ix, iy, iz) += val;
      }
    }
  }
}

void SlabSL::compute_multistep_coefficients()
{
  if (play_back and not play_cnew) return;

  // Clean coefficient matrix
  // 
  expccof[0].setZero();
    
  // Interpolate to get coefficients above
  // 
  for (int M=0; M<mfirst[mdrft]; M++) {
    
    double numer = static_cast<double>(mdrft            - dstepL[M][mdrft]);
    double denom = static_cast<double>(dstepN[M][mdrft] - dstepL[M][mdrft]);

    double b = numer/denom;	// Interpolation weights
    double a = 1.0 - b;

    for (int i=0; i<jmax; i++) {
      expccof[0].data()[i] +=
	a*expcoefL[M]->data()[i] + b*expcoefN[M]->data()[i] ;
    }
    
    // Sanity debug check
    // 
    if (a<0.0 && a>1.0) {
      std::cout << "Process " << myid
		<< ": interpolation error in multistep [a]" 
		<< std::endl;
    }
    if (b<0.0 && b>1.0) {
      std::cout << "Process " << myid
		<< ": interpolation error in multistep [b]" 
		<< std::endl;
    }
  }

  // Add coefficients at or below this level
  // 
  for (int M=mfirst[mdrft]; M<=multistep; M++) {

    for (int i=0; i<jmax; i++) {
      expccof[0].data()[i] += expcoefN[M]->data()[i];
    }
  }
}

