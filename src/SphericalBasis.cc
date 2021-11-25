#include "expand.H"

#include <chrono>
#include <sstream>

#include <SphericalBasis.H>
#include <MixtureBasis.H>

// #define TMP_DEBUG

//@{
//! These are for testing exclusively (should be set false for production)
static bool cudaAccumOverride = false;
static bool cudaAccelOverride = false;
//@}

#ifdef DEBUG
static pthread_mutex_t io_lock;
#endif

bool SphericalBasis::NewCoefs = true;

SphericalBasis::SphericalBasis(const YAML::Node& conf, MixtureBasis *m) : 
  AxisymmetricBasis(conf)
{
#if HAVE_LIBCUDA==1
  if (m) {
    throw std::runtime_error("Error in SphericalBasis: MixtureBasis logic is not yet implemented in CUDA");
  }

  // Initialize the circular storage container 
  cuda_initialize();
  initialize_cuda_sph = true;

#endif

  dof              = 3;
  mix              = m;
  geometry         = sphere;
  coef_dump        = true;
  NO_L0            = false;
  NO_L1            = false;
  EVEN_L           = false;
  EVEN_M           = false;
  NOISE            = false;
  noiseN           = 1.0e-6;
  noise_model_file = "SLGridSph.model";
  ssfrac           = 0.0;
  subset           = false;
  setup_noise      = true;
  coefMaster       = true;
  lastPlayTime     = -std::numeric_limits<double>::max();
#if HAVE_LIBCUDA==1
  cuda_aware       = true;
#endif

  try {
    if (conf["scale"]) 
      scale = conf["scale"].as<double>();
    else
      scale = 1.0;

    if (conf["rmin"]) 
      rmin = conf["rmin"].as<double>();
    else
      rmax = 0.0;

    if (conf["rmax"]) 
      rmax = conf["rmax"].as<double>();
    else
      rmax = std::numeric_limits<double>::max();

    if (conf["self_consistent"]) {
      self_consistent = conf["self_consistent"].as<bool>();
    } else
      self_consistent = true;

    if (conf["NO_L0"])   NO_L0   = conf["NO_L0"].as<bool>();
    if (conf["NO_L1"])   NO_L1   = conf["NO_L1"].as<bool>();
    if (conf["EVEN_L"])  EVEN_L  = conf["EVEN_L"].as<bool>();
    if (conf["EVEN_M"])  EVEN_M  = conf["EVEN_M"].as<bool>();
    
    if (conf["NOISE"]) {
      if (conf["NOISE"].as<bool>()) {
	NOISE = true; 
	self_consistent = false;
      }
      else NOISE = false;
    }
    
    if (conf["noiseN"])  noiseN  = conf["noiseN"].as<bool>();

    if (conf["noise_model_file"]) noise_model_file = conf["noise_model_file"].as<std::string>();

    if (conf["seedN"])   seedN   = conf["seedN"].as<int>();
    
    if (conf["ssfrac"]) {
      ssfrac = conf["ssfrac"].as<double>();
      // Check for sane value
      if (ssfrac>0.0 && ssfrac<1.0) subset = true;
    }

    if (conf["playback"]) {
      std::string file = conf["playback"].as<std::string>();
				// Check the file exists
      {
	std::ifstream test(file);
	if (not test) {
	  std::cerr << "SphericalBasis: process " << myid << " cannot open <"
		    << file << "> for reading" << std::endl;
	  MPI_Finalize();
	  exit(-1);
	}
      }

      playback = std::make_shared<SphericalCoefs>(file);

      if (playback->nmax != nmax) {
	if (myid==0) {
	  std::cerr << "SphericalBasis: nmax for playback [" << playback->nmax
		    << "] does not match specification [" << nmax << "]"
		    << std::endl;
	}
	MPI_Finalize();
	exit(-1);
      }

      if (playback->lmax != Lmax) {
	if (myid==0) {
	  std::cerr << "SphericalBasis: Lmax for playback [" << playback->lmax
		    << "] does not match specification [" << Lmax << "]"
		    << std::endl;
	}
	MPI_Finalize();
	exit(-1);
      }

      play_back = true;
    }

    if (conf["coefMaster"]) coefMaster = conf["coefMaster"].as<bool>();

  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in SphericalBasis: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (nthrds<1) nthrds=1;

  initialize();

  // Allocate coefficient matrix (one for each multistep level)
  // and zero-out contents
  //
  differ1 = vector< vector<Eigen::MatrixXd> >(nthrds);
  for (int n=0; n<nthrds; n++) {
    differ1[n] = vector<Eigen::MatrixXd>(multistep+1);
    for (int i=0; i<=multistep; i++)
      differ1[n][i].resize((Lmax+1)*(Lmax+1), nmax);
  }

  // MPI buffer space
  //
  unsigned sz = (multistep+1)*(Lmax+1)*(Lmax+1)*nmax;
  pack   = vector<double>(sz);
  unpack = vector<double>(sz);

  // Coefficient evaluation times
  // 
  expcoefN.resize(multistep+1);
  expcoefL.resize(multistep+1);
  for (int i=0; i<=multistep; i++) {
    expcoefN[i].resize((Lmax+1)*(Lmax+1));
    expcoefL[i].resize((Lmax+1)*(Lmax+1));
    for (auto & v : expcoefN[i]) {
      v = std::make_shared<Eigen::VectorXd>(nmax);
      v->setZero();
    }
    for (auto & v : expcoefL[i]) {
      v = std::make_shared<Eigen::VectorXd>(nmax);
      v->setZero();
    }
  }
    
  expcoef .resize((Lmax+1)*(Lmax+1));
  expcoef1.resize((Lmax+1)*(Lmax+1));
  
  for (auto & v : expcoef ) v = std::make_shared<Eigen::VectorXd>(nmax);
  for (auto & v : expcoef1) v = std::make_shared<Eigen::VectorXd>(nmax);
  
  expcoef0.resize(nthrds);
  for (auto & t : expcoef0) {
    t.resize((Lmax+1)*(Lmax+1));
    for (auto & v : t) v = std::make_shared<Eigen::VectorXd>(nmax);
  }

  // Allocate normalization matrix

  normM.resize(Lmax+1, nmax);

  for (int l=0; l<=Lmax; l++) {
    for (int n=0; n<nmax; n++) {
      normM(l, n) = 1.0;
    }
  }

  if (pcavar or pcaeof) {
    muse1 = vector<double>(nthrds, 0.0);
    muse0 = 0.0;

    pthread_mutex_init(&cc_lock, NULL);
  }

  // Potential and deriv matrices
  //
  normM.resize(Lmax+1, nmax);
  krnl. resize(Lmax+1, nmax);
  dend. resize(Lmax+1, nmax);

  potd.resize(nthrds);
  dpot.resize(nthrds);

  for (auto & v : potd) v.resize(Lmax+1, nmax);
  for (auto & v : dpot) v.resize(Lmax+1, nmax);

  // Sin, cos, legendre
  //
  cosm .resize(nthrds);
  sinm .resize(nthrds);
  legs .resize(nthrds);
  dlegs.resize(nthrds);

  for (auto & v : cosm)  v.resize(Lmax+1);
  for (auto & v : sinm)  v.resize(Lmax+1);
  for (auto & v : legs)  v.resize(Lmax+1, Lmax+1);
  for (auto & v : dlegs) v.resize(Lmax+1, Lmax+1);

  // Work vectors
  //
  u. resize(nthrds);
  du.resize(nthrds);

  for (auto & v : u)  v.resize(nmax+1);
  for (auto & v : du) v.resize(nmax+1);

  // Factorial matrix
  //
  factorial.resize(Lmax+1, Lmax+1);

  for (int l=0; l<=Lmax; l++) {
    for (int m=0; m<=l; m++) 
      factorial(l, m) = factrl(l-m)/factrl(l+m);
  }

  firstime_coef  = true;
  firstime_accel = true;

#ifdef DEBUG
  pthread_mutex_init(&io_lock, NULL);
#endif
}

void SphericalBasis::setup(void)
{				// Call normalization and kernel
  for (int l=0; l<=Lmax; l++) {	// with current binding from derived class
    for (int n=0; n<nmax; n++) {
      normM (l, n) = norm(n, l);
      krnl  (l, n) = knl(n, l);
      sqnorm(l, n) = sqrt(normM(l, n));
    }
  }

  if (NOISE) compute_rms_coefs();
}  


SphericalBasis::~SphericalBasis()
{
  if (pcavar or pcaeof) {
    pthread_mutex_destroy(&cc_lock);
  }

#if HAVE_LIBCUDA==1
  if (component->cudaDevice>=0) destroy_cuda();
#endif
}

void SphericalBasis::initialize()
{
				// Do nothing
}

void SphericalBasis::check_range()
{
				// Do nothing
}

void SphericalBasis::get_acceleration_and_potential(Component* C)
{
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = std::make_shared<nvTracer>("SphericalBasis::get_acceleration");

#ifdef DEBUG
  cout << "Process " << myid 
       << ": in SphericalBasis::get_acceleration_and_potential" << endl;
#endif
				
  cC = C;			// "Register" component
  nbodies = cC->Number();	// And compute number of bodies

  if (NOISE) update_noise();

  //======================================
  // Determine potential and acceleration 
  //======================================

  MPL_start_timer();

  determine_acceleration_and_potential();

  MPL_stop_timer();

  // Clear external potential flag
  use_external = false;

  //======================================
  // Dump coefficients for debugging
  //======================================

#ifdef DEBUG
  if (myid==0) {
    if ( (multistep && mstep==0) || !multistep) {
      static int cnt = 0;
      ostringstream sout;
      sout << "SphericalBasis.debug." << runtag << "." << cnt++;
      ofstream out(sout.str().c_str());
      if (out) dump_coefs_all(out);
    }
  }
#endif

}


void * SphericalBasis::determine_coefficients_thread(void * arg)
{
  double r, r2, rs, facL, fac1, fac2, costh, phi, mass;
  double fac0=4.0*M_PI;
  double xx, yy, zz;

  unsigned nbodies = component->levlist[mlevel].size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double adb = component->Adiabatic();
  std::vector<double> wk(nmax);

#ifdef DEBUG
  pthread_mutex_lock(&io_lock);
  cout << "Process " << myid 
       << ", " << id
       << ": in determine_coefficients_thread"
       << ", rmax=" << rmax << endl;
  pthread_mutex_unlock(&io_lock);
#endif

  thread_timing_beg(id);

  vector<double> ctr;
  if (mix) mix->getCenter(ctr);

				// Compute potential using a 
				// subset of particles
  if (subset) nend = (int)floor(ssfrac*nend);

  unsigned whch = 0;		// For PCA jacknife

  for (int i=nbeg; i<nend; i++) {

    int indx = component->levlist[mlevel][i];

    if (component->freeze(indx)) continue;

    
    mass = component->Mass(indx) * adb;
				// Adjust mass for subset
    if (subset) mass /= ssfrac;
    
    if (mix) {
      xx = component->Pos(indx, 0, Component::Local) - ctr[0];
      yy = component->Pos(indx, 1, Component::Local) - ctr[1];
      zz = component->Pos(indx, 2, Component::Local) - ctr[2];
    } else {
      xx = component->Pos(indx, 0, Component::Local | Component::Centered);
      yy = component->Pos(indx, 1, Component::Local | Component::Centered);
      zz = component->Pos(indx, 2, Component::Local | Component::Centered);
    }

    r2 = (xx*xx + yy*yy + zz*zz);
    r = sqrt(r2) + DSMALL;
      
    if (r>=rmin and r<=rmax) {

      use[id]++;
      costh = zz/r;
      phi = atan2(yy,xx);
      rs = r/scale;
	
      
      legendre_R(Lmax, costh, legs[id]);
      sinecosine_R(Lmax, phi, cosm[id], sinm[id]);

      get_potl(Lmax, nmax, rs, potd[id], id);

      if (compute) {
	muse1[id] += mass;
	if (pcavar) {
	  whch = indx % sampT;
	  pthread_mutex_lock(&cc_lock);
	  massT1[whch] += mass;
	  pthread_mutex_unlock(&cc_lock);
	}
      }

      //		l loop
      for (int l=0, loffset=0, iC=0; l<=Lmax; loffset+=(2*l+1), l++) {
	//		m loop
	for (int m=0, moffset=0; m<=l; m++) {

	  if (m==0) {
	    for (int n=0; n<nmax; n++) {
	      wk[n] = potd[id](l, n)*legs[id](l, m)*mass*fac0/normM(l, n);
	      (*expcoef0[id][loffset+moffset])[n] += wk[n];
	    }

	    if (compute and pcavar) {
	      pthread_mutex_lock(&cc_lock);
	      for (int n=0; n<nmax; n++) {
		(*expcoefT1[whch][iC])[n] += wk[n];
		for (int o=0; o<nmax; o++)
		  (*expcoefM1[whch][iC])(n, o) += wk[n]*wk[o]/mass;
	      }
	      pthread_mutex_unlock(&cc_lock);
	    }

	    if (compute and pcaeof) {
	      pthread_mutex_lock(&cc_lock);
	      for (int n=0; n<nmax; n++) {
		for (int o=0; o<nmax; o++) {
		  (*tvar[iC])(n, o) += wk[n]*wk[o]/mass;
		}
	      }
	      pthread_mutex_unlock(&cc_lock);
	    }

	    iC++;
	    moffset++;
	  }
	  else {
	    facL = legs[id](l, m);
	    fac1 = facL*cosm[id][m];
	    fac2 = facL*sinm[id][m];

	    for (int n=0; n<nmax; n++) {

	      wk[n] = potd[id](l, n)*mass*fac0/normM(l, n);

	      (*expcoef0[id][loffset+moffset  ])[n] += wk[n]*fac1;
	      (*expcoef0[id][loffset+moffset+1])[n] += wk[n]*fac2;
	    }

	    if (compute and pcavar) {
	      pthread_mutex_lock(&cc_lock);
	      for (int n=0; n<nmax; n++) {
		(*expcoefT1[whch][iC])[n] += wk[n]*facL;
		for (int o=0; o<nmax; o++)
		  (*expcoefM1[whch][iC])(n, o) += wk[n]*wk[o]*facL*facL/mass;
	      }
	      pthread_mutex_unlock(&cc_lock);
	    }
	    
	    if (compute and pcaeof) {
	      pthread_mutex_lock(&cc_lock);
	      for (int n=0; n<nmax; n++) {
		for (int o=0; o<nmax; o++) {
		  (*tvar[iC])(n, o) += wk[n]*wk[o]/mass;
		}
	      }
	      pthread_mutex_unlock(&cc_lock);
	    }

	    iC++;
	    moffset+=2;
	  } // m!=0

	} // m loop

      } // l loop

    } // r < rmax

  } // particle loop

  thread_timing_end(id);

  return (NULL);
}


void SphericalBasis::determine_coefficients(void)
{
  // Playback basis coefficients
  //
  if (play_back) {
				// Do we need new coefficients?
    if (tnow <= lastPlayTime) return;
    lastPlayTime = tnow;

				// Set coefficient matrix size
    expcoef.resize((Lmax+1)*(Lmax+1));
    for (auto & v : expcoef) v = std::make_shared<Eigen::VectorXd>(1, nmax);

    if (coefMaster) {

      if (myid==0) {
	auto ret = playback->interpolate(tnow);
	for (int l=0; l<(Lmax+1)*(Lmax+1); l++) {
	  for (int n=0; n<nmax; n++) (*expcoef[l])[n] = (*ret)[l][n];
	  MPI_Bcast((*expcoef[l]).data(), nmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
      } else {
	for (int l=0; l<(Lmax+1)*(Lmax+1); l++) {
	  MPI_Bcast((*expcoef[l]).data(), nmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
      }
    } else {
      auto ret = playback->interpolate(tnow);
      for (int l=0; l<(Lmax+1)*(Lmax+1); l++) {
	for (int n=0; n<nmax; n++) (*expcoef[l])[n] = (*ret)[l][n];
      }
    }
    return;
  }

  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = std::make_shared<nvTracer>("SphericalBasis::determine_coefficients");

  std::chrono::high_resolution_clock::time_point start0, start1, finish0, finish1;

  start0 = std::chrono::high_resolution_clock::now();

  // Return if we should leave the coefficients fixed
  //
  if (!self_consistent && !firstime_coef && !initializing) return;

  if (pcavar or pcaeof) {
    if (this_step >= npca0) 
      compute = (mstep == 0) && !( (this_step-npca0) % npca);
    else
      compute = false;
  }


  int loffset, moffset, use1;

  if (compute) {

    if (massT.size() == 0) {	// Allocate storage for subsampling
      if (defSampT) sampT = defSampT;
      else          sampT = floor(sqrt(component->CurTotal()));
      massT    .resize(sampT, 0);
      massT1   .resize(sampT, 0);
      
      expcoefT .resize(sampT);
      for (auto & t : expcoefT ) {
	t.resize((Lmax+1)*(Lmax+2)/2);
	for (auto & v : t) v = std::make_shared<Eigen::VectorXd>(nmax);
      }
      
      expcoefT1.resize(sampT);
      for (auto & t : expcoefT1) {
	t.resize((Lmax+1)*(Lmax+2)/2);
	for (auto & v : t) v = std::make_shared<Eigen::VectorXd>(nmax);
      }

      expcoefM .resize(sampT);
      for (auto & t : expcoefM ) {
	t.resize((Lmax+1)*(Lmax+2)/2);
	for (auto & v : t) v = std::make_shared<Eigen::MatrixXd>(nmax, nmax);
      }
      
      expcoefM1.resize(sampT);
      for (auto & t : expcoefM1) {
	t.resize((Lmax+1)*(Lmax+2)/2);
	for (auto & v : t) v = std::make_shared<Eigen::MatrixXd>(nmax, nmax);
      }

    }

    // Zero arrays?
    //
    if (mlevel==0) {
      for (int n=0; n<nthrds; n++) muse1[n] = 0.0;
      muse0 = 0.0;
      
      for (int n=0; n<nthrds; n++) use[n] = 0.0;

      if (pcavar) {
	for (auto & t : expcoefT1) { for (auto & v : t) v->setZero(); }
	for (auto & t : expcoefM1) { for (auto & v : t) v->setZero(); }
	for (auto & v : massT1)    v = 0;
      }

      if (pcaeof) {
	for (auto & v : tvar) v->setZero();
      }
    }
  }

#ifdef DEBUG
  cout << "Process " << myid << ": in <determine_coefficients>" << endl;
#endif

  //
  // Swap interpolation arrays
  //
  auto p = expcoefL[mlevel];
  
  expcoefL[mlevel] = expcoefN[mlevel];
  expcoefN[mlevel] = p;
  
  //
  // Clean arrays for current level
  //
  for (auto & v : expcoefN[mlevel]) v->setZero();
    
  for (auto & v : expcoef0) { for (auto & u : v) u->setZero(); }
    
  use1 = 0;
  if (multistep==0) used = 0;
    
#ifdef DEBUG
  cout << "Process " << myid 
       << ": in <determine_coefficients>, about to thread, lev=" 
       << mlevel << endl;
#endif

#ifdef LEVCHECK
  MPI_Barrier(MPI_COMM_WORLD);
  for (int n=0; n<numprocs; n++) {
    if (n==myid) {
      if (myid==0) cout << "-------------------------------" << endl
			<< "Level check in Spherical Basis:" << endl 
			<< "-------------------------------" << endl;
      cout << setw(4) << myid << setw(4) << mlevel;
      if (component->levlist[mlevel].size())
	cout << setw(12) << component->levlist[mlevel].size()
	     << setw(12) << component->levlist[mlevel].front()
	     << setw(12) << component->levlist[mlevel].back() << endl;
      else
	cout << setw(12) << component->levlist[mlevel].size()
	     << setw(12) << (int)(-1)
	     << setw(12) << (int)(-1) << endl;
      
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0) cout << endl;
#endif
    
  std::fill(use.begin(), use.end(), 0);

#if HAVE_LIBCUDA==1
  if (component->cudaDevice>=0 and use_cuda) {
    if (cudaAccumOverride) {
      component->CudaToParticles();
      exp_thread_fork(true);
    } else {
      start1  = std::chrono::high_resolution_clock::now();
      determine_coefficients_cuda(compute);
      DtoH_coefs(expcoef0[0]);
      finish1 = std::chrono::high_resolution_clock::now();
    }
  } else {
    exp_thread_fork(true);
  }
#else
  exp_thread_fork(true);
#endif
  
 #ifdef DEBUG
  cout << "Process " << myid << ": in <determine_coefficients>, thread returned, lev=" << mlevel << endl;
#endif

  //
  // Sum up the results from each thread
  //
  for (int i=0; i<nthrds; i++) use1 += use[i];
  for (int i=1; i<nthrds; i++) {
    for (int l=0; l<(Lmax+1)*(Lmax+1); l++) (*expcoef0[0][l]) += (*expcoef0[i][l]);
  }
  
  if (multistep==0 or tnow==resetT) {
    used += use1;
  }
  
  for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
    
    for (int m=0, moffset=0; m<=l; m++) {

      if (m==0) {
	  
	if (multistep)
	  MPI_Allreduce ( (*expcoef0[0][loffset+moffset]).data(),
			  (*expcoefN[mlevel][loffset+moffset]).data(),
			  nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	else
	  MPI_Allreduce ( (*expcoef0[0][loffset+moffset]).data(),
			  (*expcoef[loffset+moffset]).data(),
			  nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	moffset++;
	  
      } else {
	
	if (multistep) {
	  MPI_Allreduce ( (*expcoef0[0][loffset+moffset]).data(),
			  (*expcoefN[mlevel][loffset+moffset]).data(),
			  nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  
	  MPI_Allreduce ( (*expcoef0[0][loffset+moffset+1]).data(),
			  (*expcoefN[mlevel][loffset+moffset+1]).data(),
			  nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	} else {
	  MPI_Allreduce ( (*expcoef0[0][loffset+moffset]).data(),
			  (*expcoef[loffset+moffset]).data(),
			  nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  
	  MPI_Allreduce ( (*expcoef0[0][loffset+moffset+1]).data(),
			  (*expcoef[loffset+moffset+1]).data(),
			  nmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}
	moffset+=2;
      }
    }
  }
  
  //======================================
  // Last level?
  //======================================
  
  if (mlevel==multistep) {

    //======================================
    // Multistep update
    //======================================

    if (multistep) compute_multistep_coefficients();

    //======================================
    // PCA computation
    //======================================
    
    if (compute) {
      for (int i=0; i<nthrds; i++) muse0 += muse1[i];
      MPI_Allreduce ( &muse0, &muse,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      parallel_gather_coef2();
    }

    pca_hall(compute);
  }

  print_timings("SphericalBasis: coefficient timings");

#if HAVE_LIBCUDA==1
  if (component->timers) {
    auto finish0 = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> duration0 = finish0 - start0;
    std::chrono::duration<double> duration1 = finish1 - start1;

    std::cout << std::string(60, '=') << std::endl;
    std::cout << "== Coefficient evaluation [SphericalBasis] level="
	      << mlevel << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Time in CPU: " << duration0.count()-duration1.count() << std::endl;
    if (component->cudaDevice>=0 and use_cuda) {
      std::cout << "Time in GPU: " << duration1.count() << std::endl;
    }
    std::cout << std::string(60, '=') << std::endl;
  }
#endif

  //================================
  // Dump coefficients for debugging
  //================================

  //  +--- Deep debugging. Set to 'false' for production.
  //  |
  //  v
  if (false and myid==0 and mstep==0 and mlevel==multistep) {

    std::cout << std::string(60, '-') << std::endl
	      << "-- SphericalBasis T=" << std::setw(16) << tnow << std::endl
	      << std::string(60, '-') << std::endl;

    for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
      
      for (int m=0, moffset=0; m<=l; m++) {

	if (m==0) {
	  
	  for (int n=0; n<nmax; n++) {
	    std::cout << std::setw(4)  << l
		      << std::setw(4)  << m
		      << std::setw(4)  << n
		      << std::setw(18) << (*expcoef[loffset+moffset])[n]
		      << std::endl;
	  }

	  moffset++;
	  
	} else {
	  
	  for (int n=0; n<nmax; n++) {
	    std::cout << std::setw(4)  << l
		      << std::setw(4)  << m
		      << std::setw(4)  << n
		      << std::setw(18) << (*expcoef[loffset+moffset  ])[n]
		      << std::setw(18) << (*expcoef[loffset+moffset+1])[n]
		      << std::endl;
	  }

	  moffset+=2;
	}
      }
    }
    std::cout << std::string(60, '-') << std::endl;
  }

  firstime_coef = false;
}

void SphericalBasis::multistep_reset()
{
  if (play_back) return;

  used   = 0;
  resetT = tnow;
}


void SphericalBasis::multistep_update_begin()
{
  if (play_back) return;
				// Clear the update matricies
  for (int n=0; n<nthrds; n++) {
    for (int M=mfirst[mstep]; M<=multistep; M++) {
      for (int l=0; l<=Lmax*(Lmax+2); l++) {
	for (int ir=0; ir<nmax; ir++) {
	  differ1[n][M](l, ir) = 0.0;
	}
      }
    }
  }

}

void SphericalBasis::multistep_update_finish()
{
  if (play_back) return;

				// Combine the update matricies
				// from all nodes
  unsigned sz = (multistep - mfirst[mstep]+1)*(Lmax+1)*(Lmax+1)*nmax;
  unsigned offset0, offset1;

				// Zero the buffer space
				//
  for (unsigned j=0; j<sz; j++) pack[j] = unpack[j] = 0.0;

  // Pack the difference matrices
  //
  for (int M=mfirst[mstep]; M<=multistep; M++) {
    offset0 = (M - mfirst[mstep])*(Lmax+1)*(Lmax+1)*nmax;
    for (int l=0; l<=Lmax*(Lmax+2); l++) {
      offset1 = l*nmax;
      for (int n=0; n<nthrds; n++) 
	for (int ir=0; ir<nmax; ir++) 
	  pack[offset0+offset1+ir] += differ1[n][M](l, ir);
    }
  }

  MPI_Allreduce (&pack[0], &unpack[0], sz, 
		 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  //  +--- Deep debugging
  //  |
  //  v
  if (false and myid==0) {
    std::ofstream out("test_differ.sph", ios::app);
    if (out) {
      out << std::string(10+16*nmax, '-') << std::endl;
      out << "# T=" << tnow << " mstep=" << mstep << std::endl;
      for (int M=mfirst[mstep]; M<=multistep; M++) {
	offset0 = (M - mfirst[mstep])*(Lmax+1)*(Lmax+1)*nmax;
	for (int l=0; l<=Lmax*(Lmax+2); l++) {
	  offset1 = l*nmax;
	  out << std::setw(5) << M << std::setw(5) << l;
	  for (int ir=0; ir<nmax; ir++)
	    out << std::setw(16) << unpack[offset0+offset1+ir];
	  out << std::endl;
	  break;
	}
      }
      out << std::string(10+16*nmax, '-') << std::endl;
      for (int l=0; l<=Lmax*(Lmax+2); l++) {
	out << std::setw(5) << " *** " << std::setw(5) << l;
	for (int ir=0; ir<nmax; ir++)
	  out << std::setw(16) << (*expcoef[l])[ir];
	out << std::endl;
	break;
      }
      out << std::string(10+16*nmax, '-') << std::endl;
      out << std::string(10+16*nmax, '-') << std::endl;
    } else {
      std::cout << "Error opening test file <test_differ.sph> at T=" << tnow
		<< std::endl;
    }
  }
				// Update the local coefficients
				//
  for (int M=mfirst[mstep]; M<=multistep; M++) {
    offset0 = (M - mfirst[mstep])*(Lmax+1)*(Lmax+1)*nmax;
    for (int l=0; l<=Lmax*(Lmax+2); l++) {
      offset1 = l*nmax;
      for (int ir=0; ir<nmax; ir++)
	(*expcoefN[M][l])[ir] += unpack[offset0+offset1+ir];
    }
  }

}

void SphericalBasis::multistep_update(int from, int to, Component *c, int i, int id)
{
  if (play_back)    return;
  if (c->freeze(i)) return;

  double mass = c->Mass(i) * component->Adiabatic();

				// Adjust mass for subset
  if (subset) mass /= ssfrac;

  double xx = c->Pos(i, 0, Component::Local | Component::Centered);
  double yy = c->Pos(i, 1, Component::Local | Component::Centered);
  double zz = c->Pos(i, 2, Component::Local | Component::Centered);
  
  double r2 = (xx*xx + yy*yy + zz*zz);
  double  r = sqrt(r2) + DSMALL;
      
  if (r<rmax) {

    double costh = zz/r;
    double phi   = atan2(yy,xx);
    double rs    = r/scale;
    double val, val1, val2, fac0=4.0*M_PI, fac1, fac2;
    int moffset;

    legendre_R(Lmax, costh, legs[id]);
    sinecosine_R(Lmax, phi, cosm[id], sinm[id]);

    get_potl(Lmax, nmax, rs, potd[id], 0);

    //
    // l loop
    //
    for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
      //
      // m loop
      //
      for (int m=0, moffset=0; m<=l; m++) {
	if (m==0) {
	  for (int n=0; n<nmax; n++) {
	    val = potd[id](l, n)*legs[id](l, m)*mass*fac0/normM(l, n);
	    
	    differ1[id][from](loffset+moffset, n) -= val;
	    differ1[id][  to](loffset+moffset, n) += val;
	  }
	  moffset++;

	} else {
	  fac1 = legs[id](l, m)*cosm[id][m];
	  fac2 = legs[id](l, m)*sinm[id][m];

	  for (int n=0; n<nmax; n++) {
	    val1 = potd[id](l, n)*fac1*mass*fac0/normM(l, n);
	    val2 = potd[id](l, n)*fac1*mass*fac0/normM(l, n);

	    differ1[id][from](loffset+moffset  , n) -= val1;
	    differ1[id][from](loffset+moffset+1, n) -= val2;
	    differ1[id][  to](loffset+moffset  , n) += val1;
	    differ1[id][  to](loffset+moffset+1, n) += val2;
	  }
	  moffset+=2;
	}
      }
    }
  }
  
}


void SphericalBasis::compute_multistep_coefficients()
{
  if (play_back) return;

#ifdef TMP_DEBUG
  Eigen::MatrixXd tmpcoef = expcoef;
#endif

				// Clean coefficient matrix
				// 
  for (int l=0; l<(Lmax+1)*(Lmax+1); l++) expcoef[l]->setZero();

				// Interpolate to get coefficients above
				// 
  for (int M=0; M<mfirst[mdrft]; M++) {

    double numer = static_cast<double>(mdrft            - dstepL[M][mdrft]);
    double denom = static_cast<double>(dstepN[M][mdrft] - dstepL[M][mdrft]);

    double b = numer/denom;	// Interpolation weights
    double a = 1.0 - b;

    for (int l=0; l<=Lmax*(Lmax+2); l++) {
      for (int n=0; n<nmax; n++) 
	(*expcoef[l])[n] += a*(*expcoefL[M][l])[n] + b*(*expcoefN[M][l])[n];
    }
    
    //  +--- Deep debugging
    //  |
    //  v
    if (false and myid==0) {
      std::cout << std::left << std::fixed
		<< "SPH INTERP M=" << std::setw(2) << M
		<< " mstep=" << std::setw(3) << mstep
		<< " mdrft=" << std::setw(3) << mdrft
		<< " a=" << std::setw(16) << a
		<< " b=" << std::setw(16) << b
		<< " L=" << std::setw(16) << (*expcoefL[M][0])[1]
		<< " N=" << std::setw(16) << (*expcoefN[M][0])[1]
		<< " d=" << std::setw(16) << a*(*expcoefL[M][0])[1] + b*(*expcoefN[M][0])[1]
		<< " f=" << std::setw(16) << (*expcoef[0])[1]
		<< std::endl << std::right;
    }

    if (false and myid==0) {
      std::cout << "SPH interpolate:"
		<< " M="     << std::setw( 4) << M
		<< " mstep=" << std::setw( 4) << mstep 
		<< " mstep=" << std::setw( 4) << mdrft
		<< " minS="  << std::setw( 4) << dstepL[M][mdrft]
		<< " maxS="  << std::setw( 4) << dstepN[M][mdrft]
		<< " a="     << std::setw(14) << a 
		<< " b="     << std::setw(14) << b 
		<< " L01="   << std::setw(14) << (*expcoefL[M][0])[1]
		<< " N01="   << std::setw(14) << (*expcoefN[M][0])[1]
		<< " c01="   << std::setw(14) << (*expcoef[0])[1]
		<< std::endl;
    }

				// Sanity debug check
				// 
    if (a<0.0 && a>1.0) {
      cout << "Process " << myid << ": interpolation error in multistep [a]" 
	   << endl;
    }
    if (b<0.0 && b>1.0) {
      cout << "Process " << myid << ": interpolation error in multistep [b]" 
	   << endl;
    }
  }
				// Add coefficients at or below this level
				// 
  for (int M=mfirst[mdrft]; M<=multistep; M++) {

    //  +--- Deep debugging
    //  |
    //  v
    if (false and myid==0) {
      std::cout << std::left << std::fixed
		<< "SPH FULVAL M=" << std::setw(2) << M
		<< " mstep=" << std::setw(3) << mstep
		<< " mdrft=" << std::setw(3) << mdrft
		<< std::endl << std::right;
    }

    for (int l=0; l<=Lmax*(Lmax+2); l++) {
      for (int n=0; n<nmax; n++) 
	(*expcoef[l])[n] += (*expcoefN[M][l])[n];
    }
  }

  //  +--- Deep debugging
  //  |
  //  v
  if (false and myid==0) {
    std::cout << std::left << std::fixed
	      << "SPH FULVAL mstep=" << std::setw(3) << mstep
	      << "  mdrft=" << std::setw(3) << mdrft
	      << " f=" << std::setw(16) << (*expcoef[0])[1]
	      << std::endl << std::right;
  }

  if (false and myid==0) {
    std::cout << "SPH interpolated value:"
	      << " mlev="  << std::setw( 4) << mlevel
	      << " mstep=" << std::setw( 4) << mstep
	      << " mdrft=" << std::setw( 4) << mdrft
	      << " T="     << std::setw( 4) << tnow
	      << " c01="   << std::setw(14) << (*expcoef[0])[1]
	      << std::endl;
  }

#ifdef TMP_DEBUG
  if (myid==0) {

    double maxval=0.0, maxdif=0.0, maxreldif=0.0, val;
    int irmax=1, lmax=0, irmaxdif=1, lmaxdif=0, irmaxrel=1, lmaxrel=0;
    for (int ir=0; ir<nmax; ir++) {
      for (int l=0; l<=Lmax*(Lmax+2); l++) {
	val = expcoef[l][ir] - tmpcoef[l][ir];

	if (fabs(expcoef[l][ir]) > maxval) {
	  maxval = expcoef[l][ir];
	  irmax = ir;
	  lmax = l;
	}
	if (fabs(val) > maxdif) {
	  maxdif = val;
	  irmaxdif = ir;
	  lmaxdif = l;
	}

	if (fabs(val/expcoef[l][ir]) > maxreldif) {
	  maxreldif = val/expcoef[l][ir];
	  irmaxrel = ir;
	  lmaxrel = l;
	}
      }
    }

    ofstream out("coefs.diag", ios::app);
    out << setw(15) << tnow
	<< setw(10) << mstep
	<< setw(10) << mdrft
	<< setw(18) << maxval
	<< setw(8)  << lmax
	<< setw(8)  << irmax
	<< setw(18) << maxdif
	<< setw(8)  << lmaxdif
	<< setw(8)  << irmaxdif
	<< setw(18) << maxreldif
	<< setw(8)  << lmaxrel
	<< setw(8)  << irmaxrel
	<< endl;
  }
#endif

}

void * SphericalBasis::determine_acceleration_and_potential_thread(void * arg)
{
  int l, loffset, moffset, m, ioff, indx, nbeg, nend;
  unsigned nbodies;
  double r, rs, r0=0.0, fac, fac1, fac2, fac3, fac4, costh, phi, dp;
  double potr, potl, pott, potp, p, pc, dpc, ps, dps, facp, facdp;
  double dfac=0.25/M_PI;

  double pos[3];
  double xx, yy, zz, mfactor=1.0;

  vector<double> ctr;
  if (mix) mix->getCenter(ctr);

  int id = *((int*)arg);

  thread_timing_beg(id);

  // If we are multistepping, compute accel only at or above <mlevel>
  //
  for (int lev=mlevel; lev<=multistep; lev++) {

    nbodies = cC->levlist[lev].size();

    if (nbodies==0) continue;

    nbeg = nbodies*(id  )/nthrds;
    nend = nbodies*(id+1)/nthrds;

#ifdef DEBUG
  pthread_mutex_lock(&io_lock);
  cout << "Process " << myid << ": in thread"
       << " id=" << id 
       << " level=" << lev
       << " nbeg=" << nbeg << " nend=" << nend << endl;
  pthread_mutex_unlock(&io_lock);
#endif

  for (int i=nbeg; i<nend; i++) {

    indx = cC->levlist[lev][i];

      if (cC->freeze(indx)) continue;

      if (mix) {
	if (use_external) {
	  cC->Pos(pos, indx, Component::Inertial);
	  component->ConvertPos(pos, Component::Local);
	} else
	  cC->Pos(pos, indx, Component::Local);

	mfactor = mix->Mixture(pos);
	xx = pos[0] - ctr[0];
	yy = pos[1] - ctr[1];
	zz = pos[2] - ctr[2];
      } else {
	if (use_external) {
	  cC->Pos(pos, indx, Component::Inertial);
	  component->ConvertPos(pos, Component::Local | Component::Centered);
	} else
	  cC->Pos(pos, indx, Component::Local | Component::Centered);
	
	xx = pos[0];
	yy = pos[1];
	zz = pos[2];
      }	

      fac1 = dfac * mfactor;

      r = sqrt(xx*xx + yy*yy + zz*zz) + DSMALL;
      costh = zz/r;
      rs = r/scale;
      phi = atan2(yy, xx);

      dlegendre_R (Lmax, costh, legs[id], dlegs[id]);
      sinecosine_R(Lmax, phi,   cosm[id], sinm [id]);

      if (r>rmax) {
	ioff = 1;
	r0 = r;
	r = rmax;
	rs = r/scale;
      }
      else
	ioff = 0;


      potl = potr = pott = potp = 0.0;
      
      get_dpotl(Lmax, nmax, rs, potd[id], dpot[id], id);

      if (!NO_L0) {
	get_pot_coefs_safe(0, *expcoef[0], p, dp, potd[id], dpot[id]);
	if (ioff) {
	  p *= rmax/r0;
	  dp = -p/r0;
	}
	potl = fac1*p;
	potr = fac1*dp;
      }
      
      //		l loop
      //		------
      for (l=1, loffset=1; l<=Lmax; loffset+=(2*l+1), l++) {

				// Suppress L=1 terms?
	if (NO_L1 && l==1) continue;
	
				// Suppress odd L terms?
	if (EVEN_L && (l/2)*2 != l) continue;

	//		m loop
	//		------
	for (m=0, moffset=0; m<=l; m++) {
	  
				// Suppress odd M terms?
	  if (EVEN_M && (m/2)*2 != m) continue;

	  fac1 = (2.0*l+1.0)/(4.0*M_PI) * mfactor;
	  if (m==0) {
	    fac2 = fac1*legs[id](l, m);
	    get_pot_coefs_safe(l, *expcoef[loffset+moffset], p, dp,
			       potd[id], dpot[id]);
	    if (ioff) {
	      p *= pow(rmax/r0,(double)(l+1));
	      dp = -p/r0 * (l+1);
	    }
	    potl += fac2*p;
	    potr += fac2*dp;
	    pott += fac1*dlegs[id](l, m)*p;
	    moffset++;
	  }
	  else {
	    fac2 = 2.0 * fac1 * factorial(l, m);
	    fac3 = fac2 * legs[id](l, m);
	    fac4 = fac2 * dlegs[id](l, m);
	    get_pot_coefs_safe(l, *expcoef[loffset+moffset], pc, dpc,
			       potd[id], dpot[id]);
	    get_pot_coefs_safe(l, *expcoef[loffset+moffset+1], ps, dps,
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

      cC->AddAcc(indx, 0, -(potr*xx/r - pott*xx*zz/(r*r*r)) );
      cC->AddAcc(indx, 1, -(potr*yy/r - pott*yy*zz/(r*r*r)) );
      cC->AddAcc(indx, 2, -(potr*zz/r + pott*fac/(r*r*r))   );
      if (fac > DSMALL) {
	cC->AddAcc(indx, 0,  potp*yy/fac );
	cC->AddAcc(indx, 1, -potp*xx/fac );
      }
      if (use_external)
	cC->AddPotExt(indx, potl);
      else
	cC->AddPot(indx, potl);
    }

  }

  thread_timing_end(id);

  return (NULL);
}


void SphericalBasis::determine_acceleration_and_potential(void)
{
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = std::make_shared<nvTracer>("SphericalBasis::determine_acceleration");

  std::chrono::high_resolution_clock::time_point start0, start1, finish0, finish1;
  start0 = std::chrono::high_resolution_clock::now();

#ifdef DEBUG
  cout << "Process " << myid << ": in determine_acceleration_and_potential\n";
#endif

  if (use_external == false) {

    if (multistep && (self_consistent || initializing)) {
      compute_multistep_coefficients();
    }

  }

#if HAVE_LIBCUDA==1
  if (use_cuda and cC->cudaDevice>=0 and cC->force->cudaAware()) {
    if (cudaAccelOverride) {
      cC->CudaToParticles();
      exp_thread_fork(false);
      cC->ParticlesToCuda();
    } else {
      start1 = std::chrono::high_resolution_clock::now();
      //
      // Copy coefficients from this component to device
      //
      HtoD_coefs(expcoef);
      //
      // Do the force computation
      //
      determine_acceleration_cuda();
      
      finish1 = std::chrono::high_resolution_clock::now();
    }
  } else {

    exp_thread_fork(false);

  }
#else

  exp_thread_fork(false);

#endif

#ifdef DEBUG
  if (not use_cuda) {
    cout << "SphericalBasis: process " << myid << " returned from fork" << endl;
    cout << "SphericalBasis: process " << myid << " name=<" << cC->name << ">";

    if (cC->Particles().size()) {

      unsigned long imin = std::numeric_limits<unsigned long>::max();
      unsigned long imax = 0, kmin = kmin, kmax = 0;
      for (auto p : cC->Particles()) {
	imin = std::min<unsigned long>(imin, p.first);
	imax = std::max<unsigned long>(imax, p.first);
	kmin = std::min<unsigned long>(kmin, p.second->indx);
	kmax = std::max<unsigned long>(kmax, p.second->indx);
      }
      
      cout << " bodies ["
	   << kmin << ", " << kmax << "], ["
	   << imin << ", " << imax << "]"
	   << " #=" << cC->Particles().size() << endl;
      
    } else
      cout << " zero bodies!" << endl;
  }
#endif

  print_timings("SphericalBasis: acceleration timings");

#if HAVE_LIBCUDA==1
  if (component->timers) {
    finish0 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration0 = finish0 - start0;
    std::chrono::duration<double> duration1 = finish1 - start1;
  
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "== Force evaluation [SphericalBasis::" << cC->name
	      << "] level=" << mlevel << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Time in CPU: " << duration0.count()-duration1.count() << std::endl;
    if (cC->cudaDevice>=0) {
      std::cout << "Time in GPU: " << duration1.count() << std::endl;
    }
    std::cout << std::string(60, '=') << std::endl;
  }
#endif
}


void SphericalBasis::get_pot_coefs(int l, const Eigen::VectorXd& coef,
				   double& p, double& dp)
{
  double pp, dpp;

  pp = dpp = 0.0;

  for (int i=0; i<nmax; i++) {
    pp  += potd[0](l, i) * coef[i];
    dpp += dpot[0](l, i) * coef[i];
  }

  p = -pp;
  dp = -dpp;
}


void SphericalBasis::get_pot_coefs_safe(int l, const Eigen::VectorXd& coef, 
					double& p, double& dp,
					Eigen::MatrixXd& potd1, Eigen::MatrixXd& dpot1)
{
  double pp, dpp;
  int i;

  pp = dpp = 0.0;

  for (int i=0; i<nmax; i++) {
    pp  += potd1(l, i) * coef[i];
    dpp += dpot1(l, i) * coef[i];
  }

  p = -pp;
  dp = -dpp;
}


void SphericalBasis::get_dens_coefs(int l, Eigen::VectorXd& coef, double& p)
{
  double pp;

  pp = 0.0;

  for (int i=0; i<nmax; i++)
    pp  += dend(l, i) * coef[i];

  p = pp;
}
				// Dump coefficients to a file

void SphericalBasis::dump_coefs(ostream& out)
{
  if (NewCoefs) {

    // This is a node of simple {key: value} pairs.  More general
    // content can be added as needed.
    //
    YAML::Node node;

    node["id"    ] = id;
    node["time"  ] = tnow;
    node["scale" ] = scale;
    node["nmax"  ] = nmax;
    node["lmax"  ] = Lmax;
    node["normed"] = true;

    // Serialize the node
    //
    YAML::Emitter y; y << node;
    
    // Get the size of the string
    //
    unsigned int hsize = strlen(y.c_str());
    
    // Write magic #
    //
    out.write(reinterpret_cast<const char *>(&cmagic),   sizeof(unsigned int));

    // Write YAML string size
    //
    out.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));
    
    // Write YAML string
    //
    out.write(reinterpret_cast<const char *>(y.c_str()), hsize);

    double z;

    for (int ir=0; ir<nmax; ir++) {
      for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
	double fac1 = (2.0*l+1.0)/(4.0*M_PI);
	for (int m=0, moffset=0; m<=l; m++) {
	  double fac2 = sqrt(fac1*factorial(l, m));
	  if (m==0) {
	    out.write((char *)&(z=fac2*(*expcoef[loffset+moffset+0])[ir]), sizeof(double));
	    moffset += 1;
	  } else {
	    out.write((char *)&(z=fac2*(*expcoef[loffset+moffset+0])[ir]), sizeof(double));
	    out.write((char *)&(z=fac2*(*expcoef[loffset+moffset+1])[ir]), sizeof(double));
	    moffset += 2;
	  }
	}
      }
    }

  } else {

    std::ostringstream sout;
    sout << id;

    char buf[64];
    for (int i=0; i<64; i++) {
      if (i<sout.str().length())  buf[i] = sout.str().c_str()[i];
      else                        buf[i] = '\0';
    }

    out.write((char *)&buf, 64*sizeof(char));
    out.write((char *)&tnow, sizeof(double));
    out.write((char *)&scale, sizeof(double));
    out.write((char *)&nmax, sizeof(int));
    out.write((char *)&Lmax, sizeof(int));

    for (int ir=0; ir<nmax; ir++) {
      for (int l=0; l<=Lmax*(Lmax+2); l++)
	out.write((char *)&(*expcoef[l])[ir], sizeof(double));
    }
  }

}


void SphericalBasis::determine_fields_at_point
(double x, double y, double z, 
 double *tdens0, double *tpotl0, 
 double *tdens,  double *tpotl, 
 double *tpotX,  double *tpotY, double *tpotZ)
{
  double r = sqrt(x*x + y*y + z*z) + 1.0e-18;
  double theta = acos(z/r);
  double phi   = atan2(y, x);
  double cth   = cos(theta), sth = sin(theta);
  double cph   = cos(phi),   sph = sin(phi);
  double tpotr, tpott, tpotp;

  determine_fields_at_point_sph(r, theta, phi,
				tdens0, tpotl0, 
				tdens, tpotl, 
				&tpotr,	&tpott, &tpotp);

  *tpotX = tpotr*sth*cph + tpott*cth*cph - tpotp*sph;
  *tpotY = tpotr*sth*sph + tpott*cth*sph + tpotp*cph;
  *tpotZ = tpotr*cth     - tpott*sth;
}


void SphericalBasis::determine_fields_at_point_cyl
(double R, double z, double phi, 
 double *tdens0, double *tpotl0, 
 double *tdens, double *tpotl, 
 double *tpotR, double *tpotz, double *tpotp)
{
  double r = sqrt(R*R + z*z) + 1.0e-18;
  double theta = acos(z/r);
  double tpotr, tpott;

  determine_fields_at_point_sph(r, theta, phi, 
				tdens0, tpotl0, 
				tdens, tpotl, 
				&tpotr,	&tpott, tpotp);

  *tpotR = tpotr*sin(theta) - tpott*cos(theta);
  *tpotz = tpotr*cos(theta) + tpott*sin(theta);
}

void SphericalBasis::determine_fields_at_point_sph
(double r, double theta, double phi, 
 double *tdens0, double *tpotl0, 
 double *tdens, double *tpotl, 
 double *tpotr, double *tpott, double *tpotp)
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
  get_dens_coefs(0, *expcoef[0], dens);
  dens *= dfac*dfac;

  get_pot_coefs(0, *expcoef[0], p, dp);
  potl = fac1*p;
  potr = fac1*dp;
  pott = potp = 0.0;
  
  *tdens0 = dens;
  *tpotl0 = potl;

  // l loop
    
  for (l=1, loffset=1; l<=Lmax; loffset+=(2*l+1), l++) {
    
    // m loop
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
	fac2 = fac1*legs[0](l, m);
	get_dens_coefs(l, *expcoef[loffset+moffset], p);
	dens += dfac*fac2*p;
	get_pot_coefs(l, *expcoef[loffset+moffset], p, dp);
	potl += fac2*p;
	potr += fac2*dp;
	pott += fac1*dlegs[0](l, m)*p;
	moffset++;
      }
      else {
	fac2 = 2.0 * fac1 * factorial(l, m);
	fac3 = fac2 * legs[0](l, m);
	fac4 = fac2 * dlegs[0](l, m);
	
	get_dens_coefs(l, *expcoef[loffset+moffset], pc);
	get_dens_coefs(l, *expcoef[loffset+moffset+1], ps);
	dens += dfac*fac3*(pc*cosm[0][m] + ps*sinm[0][m]);
	
	get_pot_coefs(l, *expcoef[loffset+moffset], pc, dpc);
	get_pot_coefs(l, *expcoef[loffset+moffset+1], ps, dps);
	potl += fac3*(pc*cosm[0][m] + ps*sinm[0][m]);
	potr += fac3*(dpc*cosm[0][m] + dps*sinm[0][m]);
	pott += fac4*(pc*cosm[0][m] + ps*sinm[0][m]);
	potp += fac3*(-pc*sinm[0][m] + ps*cosm[0][m])*m;
	moffset +=2;
      }
    }
  }

  *tdens0 /= scale*scale*scale;
  *tpotl0 /= scale;

  *tdens = dens/(scale*scale*scale);
  *tpotl = potl/scale;
  *tpotr = potr/(scale*scale);
  *tpott = pott/scale;
  *tpotp = potp/scale;
  
}


void SphericalBasis::compute_rms_coefs(void)
{
  meanC.resize(nmax);
  rmsC.resize(Lmax+1, nmax);

  meanC.setZero();
  rmsC.setZero();

  const int numg = 100;
  LegeQuad qe(numg);

  SphericalModelTable modl(noise_model_file);
  double rmin = modl.get_min_radius();
  double rmax = modl.get_max_radius();
  double del = rmax - rmin;
  double r, rs;

  for (int i=0; i<numg; i++) {
    r = rmin + del*qe.knot(i);
    rs = r / scale;

    get_potl(Lmax, nmax, rs, potd[0], 0);

    for(int l=0; l<=Lmax; l++) {
      
      for (int n=0; n<nmax; n++) {

	if (l==0)
	  meanC[n] += del * qe.weight(i) * r * r * potd[0](l, n)/scale *
	    modl.get_density(r);

	rmsC(l, n) += del * qe.weight(i) * r * r * potd[0](l, n)/scale *
	  potd[0](l, n)/scale * modl.get_density(r);
      }
    }
  }

  double fac, fac1;
  double mtot = modl.get_mass(rmax);

  for(int l=0; l<=Lmax; l++) {

    fac1 = (4.0*M_PI)/(2.0*l+1.0);

    for (int n=0; n<nmax; n++) {
      fac = normM(l, n);
      if (l==0) meanC[n] *= 4.0*M_PI*fac1/fac;
      rmsC(l, n) *= mtot*4.0*M_PI*fac1*fac1/(fac*fac);
    }
  }

}


void SphericalBasis::update_noise(void)
{

  if (setup_noise) {
    setup_noise = false;	// Only do this initialization once
    
    rgen.seed(seedN);		// Want the same seed on each process

    if (myid==0) {		// Diagnostic output
      ofstream out("rmscoef.dat");

      for(int l=0; l<=Lmax; l++) {

	out << "# L=" << l << endl;
	for (int n=0; n<nmax; n++) {
	  if (l==0)
	    out << setw(5)  << n
		<< setw(16) << (*expcoef[l])[n]
		<< setw(16) << meanC[n]
		<< setw(16) << rmsC(l, n)
		<< setw(16) << rmsC(l, n) - meanC[n]*meanC[n]
		<< endl;
	  else
	    out << setw(5)  << n
		<< setw(16) << (*expcoef[l])[n]
		<< setw(16) << 0.0
		<< setw(16) << rmsC(l, n)
		<< setw(16) << rmsC(l, n)
		<< endl;
	}
	out << endl;
      }
    }
  }


				// l loop
  for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
				// m loop
    for (int m=0, moffset=0; m<=l; m++) {

      if (m==0) {
	for (int n=0; n<nmax; n++) {
	  (*expcoef[loffset+moffset])[n] = 
	    sqrt(fabs(rmsC(l, n) - meanC[n]*meanC[n])/factorial(l, m)/noiseN)*nrand(rgen);
	  if (l==0) (*expcoef[l])[n] += meanC[n];
	}
	moffset++;
      }
      else {
	for (int n=0; n<nmax; n++) {
	  (*expcoef[loffset+moffset+0])[n] = 
	    sqrt(0.5*fabs(rmsC(l, n) - meanC[n]*meanC[n])/factorial(l, m)/noiseN)*nrand(rgen);
	  (*expcoef[loffset+moffset+1])[n] = 
	    sqrt(0.5*fabs(rmsC(l, n) - meanC[n]*meanC[n])/factorial(l, m)/noiseN)*nrand(rgen);
	}
	moffset+=2;
      }
    }
  }

}


void SphericalBasis::dump_coefs_all(ostream& out)
{
  out << setw(10) << "Time:"   << setw(18) << tnow      << endl
      << setw(10) << "Scale:"  << setw(18) << scale     << endl
      << setw(10) << "Nmax:"   << setw(8)  << nmax      << endl
      << setw(10) << "Lmax:"   << setw(8)  << Lmax      << endl
      << setw(10) << "Levels:" << setw(8)  << multistep << endl;
  
  out << setw(70) << setfill('=') << "=" << endl << setfill(' ');
  out << "Total" << endl;
  out << setw(70) << setfill('=') << "=" << endl << setfill(' ');

  for (int ir=0; ir<nmax; ir++) {
    for (int l=0; l<=Lmax*(Lmax+2); l++)
      out << setw(18) << (*expcoef[l])[ir];
    out << endl;
  }
  
  out << endl;

  if (multistep) {

    for (int M=0; M<=multistep; M++) {
      out << setw(70) << setfill('=') << "=" << endl << setfill(' ');
      out << "Level " << M << endl;
      out << setw(70) << setfill('=') << "=" << endl << setfill(' ');

      for (int ir=0; ir<nmax; ir++) {
	for (int l=0; l<=Lmax*(Lmax+2); l++)
	  out << setw(18) << (*expcoefN[M][l])[ir];
	out << endl;
      }
    }
  }

}

