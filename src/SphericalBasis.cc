#include "expand.H"

#include <filesystem>
#include <sstream>
#include <chrono>
#include <string>
#include <vector>
#include <set>

#include "SphericalBasis.H"
#include "MixtureBasis.H"

// #define TMP_DEBUG
// #define MULTI_DEBUG

//@{
//! These are for testing exclusively (should be set false for production)
static bool cudaAccumOverride = false;
static bool cudaAccelOverride = false;
//@}

#ifdef DEBUG
static pthread_mutex_t io_lock;
#endif

bool SphericalBasis::NewCoefs = true;

const std::set<std::string>
SphericalBasis::valid_keys = {
  "scale",
  "rmin",
  "rmax",
  "self_consistent",
  "FIX_L0",
  "NO_L0",
  "NO_L1",
  "EVEN_L",
  "EVEN_M",
  "M0_ONLY",
  "NOISE",
  "noiseN",
  "noise_model_file",
  "seedN",
  "ssfrac",
  "playback",
  "coefCompute",
  "coefMaster",
  "orthocheck",
  "subsampleFloat",
  "totalCovar",
  "fullCovar"
};

SphericalBasis::SphericalBasis(Component* c0, const YAML::Node& conf, MixtureBasis *m) : 
  AxisymmetricBasis(c0, conf)
{
#if HAVE_LIBCUDA==1
  if (m) {
    throw GenericError("Error in SphericalBasis: MixtureBasis logic is not yet "
		       "implemented in CUDA", __FILE__, __LINE__, 1030, false);
  }

  // Initialize the circular storage container 
  cuda_initialize();
  initialize_cuda_sph = true;

#endif

  dof              = 3;
  mix              = m;
  geometry         = sphere;
  coef_dump        = true;
  FIX_L0           = false;
  NO_L0            = false;
  NO_L1            = false;
  EVEN_L           = false;
  EVEN_M           = false;
  M0_only          = false;
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
  ortho_check      = false;

  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);

  // Assign values from YAML
  //
  try {
    if (conf["scale"]) 
      scale = conf["scale"].as<double>();
    else
      scale = 1.0;

    if (conf["rmin"]) 
      rmin = conf["rmin"].as<double>();
    else
      rmin = 0.0;

    if (conf["rmax"]) 
      rmax = conf["rmax"].as<double>();
    else
      rmax = std::numeric_limits<double>::max();

    if (conf["self_consistent"]) {
      self_consistent = conf["self_consistent"].as<bool>();
    } else
      self_consistent = true;

    if (conf["FIX_L0"])  FIX_L0  = conf["FIX_L0"].as<bool>();
    if (conf["NO_L0"])   NO_L0   = conf["NO_L0"].as<bool>();
    if (conf["NO_L1"])   NO_L1   = conf["NO_L1"].as<bool>();
    if (conf["EVEN_L"])  EVEN_L  = conf["EVEN_L"].as<bool>();
    if (conf["EVEN_M"])  EVEN_M  = conf["EVEN_M"].as<bool>();
    if (conf["M0_ONLY"]) M0_only = conf["M0_ONLY"].as<bool>();

    if (conf["fullCovar"]) {
      fullCovar = conf["fullCovar"].as<bool>();
    }
    
    if (conf["totalCovar"]) {
      totlCovar = conf["totalCovar"].as<bool>();
      if (totlCovar) fullCovar = true;
    }

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
	  throw std::runtime_error("SphericalBasis: file error");
	}
      }

      // This creates the Coefs instance
      playback = std::dynamic_pointer_cast<CoefClasses::SphCoefs>(CoefClasses::Coefs::factory(file));

      // Check to make sure that has been created
      if (not playback) {
	throw GenericError("SphericalBasis: failure in downcasting",
			   __FILE__, __LINE__, 1031, false);
      }

      // Set tolerance to 2 master time steps
      playback->setDeltaT(dtime*2);

      if (playback->nmax() != nmax) {
	if (myid==0) {
	  std::cerr << "SphericalBasis: nmax for playback [" << playback->nmax()
		    << "] does not match specification [" << nmax << "]"
		    << std::endl;
	}
	throw std::runtime_error("SphericalBasis: parameter mismatch");
      }

      if (playback->lmax() != Lmax) {
	if (myid==0) {
	  std::cerr << "SphericalBasis: Lmax for playback [" << playback->lmax()
		    << "] does not match specification [" << Lmax << "]"
		    << std::endl;
	}
	throw std::runtime_error("SphericalBasis: parameter mismatch");
      }

      play_back = true;

      if (conf["coefCompute"]) play_cnew = conf["coefCompute"].as<bool>();

      if (conf["coefMaster"]) coefMaster = conf["coefMaster"].as<bool>();

      if (myid==0) {
	std::cout << "---- Playback is ON for Component " << component->name
		  << " using Force " << component->id << std::endl;
	if (coefMaster)
	  std::cout << "---- Playback will use MPI master" << std::endl;

	if (play_cnew)
	  std::cout << "---- New coefficients will be computed from particles on playback" << std::endl;
      }
    }
    // END: playback config

    if (conf["orthocheck"]) ortho_check = conf["orthocheck"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in SphericalBasis: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("SphericalBasis: error parsing YAML");
  }

  if (nthrds<1) nthrds=1;

  initialize();

  // Allocate coefficient matrix (one for each multistep level)
  // and zero-out contents
  //
  differ1 = std::vector< std::vector<Eigen::MatrixXd> >(nthrds);
  for (int n=0; n<nthrds; n++) {
    differ1[n].resize(multistep+1);
    for (int i=0; i<=multistep; i++)
      differ1[n][i].resize((Lmax+1)*(Lmax+1), nmax);
  }

  // MPI buffer space
  //
  unsigned sz = (multistep+1)*(Lmax+1)*(Lmax+1)*nmax;
  pack  .resize(sz);
  unpack.resize(sz);

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
    for (int m=0; m<=l; m++) {
      factorial(l, m) = sqrt( (2.0*l+1.0)/(4.0*M_PI) * factrl(l-m)/factrl(l+m) );
      if (m) factorial(l, m) *= M_SQRT2;
    }
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
      krnl  (l, n) = knl (n, l);
      sqnorm(l, n) = sqrt(normM(l, n));
    }
  }

  if (NOISE) compute_rms_coefs();
  if (ortho_check) biorthogonality_check();
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
  // Biorthgonal normalization factor
  //
  double fac0=-4.0*M_PI;

  // Partition bodies and get features
  //
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

    
    double mass = component->Mass(indx) * adb;
				// Adjust mass for subset
    if (subset) mass /= ssfrac;
    
    double xx, yy, zz;
    if (mix) {
      xx = component->Pos(indx, 0, Component::Local) - ctr[0];
      yy = component->Pos(indx, 1, Component::Local) - ctr[1];
      zz = component->Pos(indx, 2, Component::Local) - ctr[2];
    } else {
      xx = component->Pos(indx, 0, Component::Local | Component::Centered);
      yy = component->Pos(indx, 1, Component::Local | Component::Centered);
      zz = component->Pos(indx, 2, Component::Local | Component::Centered);
    }

    double r2 = (xx*xx + yy*yy + zz*zz);
    double r  = sqrt(r2) + DSMALL;
      
    if (r>=rmin and r<=rmax) {

      use[id]++;
      double costh = zz/r;
      double phi = atan2(yy,xx);
      double rs = r/scale;
      
      legendre_R(Lmax, costh, legs[id]);
      sinecosine_R(Lmax, phi, cosm[id], sinm[id]);

      get_potl(Lmax, nmax, rs, potd[id], id);

      if (compute) {
	muse1[id] += mass;
	if (pcavar) {
	  whch = indx % sampT;
	  pthread_mutex_lock(&cc_lock);
	  countT1[whch] += 1;
	  massT1 [whch] += mass;
	  pthread_mutex_unlock(&cc_lock);
	}
      }

      //		l loop
      for (int l=0, loffset=0, iC=0; l<=Lmax; loffset+=(2*l+1), l++) {
	//		m loop
	for (int m=0, moffset=0; m<=l; m++) {

	  double facL = factorial(l, m) * legs[id](l, m);

	  if (m==0) {
	    for (int n=0; n<nmax; n++) {
	      wk[n] = potd[id](l, n)*facL*mass*fac0/sqnorm(l, n);
	      (*expcoef0[id][loffset+moffset])[n] += wk[n];
	    }

	    if (compute and pcavar) {
	      pthread_mutex_lock(&cc_lock);
	      for (int n=0; n<nmax; n++) {
		if (std::fabs(wk[n]) < 1.0e-40) {
		  std::cout << "Coef contribution is " << wk[n] << " for "
			    << l << " " << m << " " << n << std::endl;
		}
		(*expcoefT1[whch][iC])[n] += wk[n];
		for (int o=0; o<nmax; o++) {
		  if (std::fabs(wk[n]*wk[o]/mass) < 1.0e-40) {
		    std::cout << "Covar contribution is " << wk[n]*wk[o]/mass << " for "
			      << l << " " << m << " "
			      << n << " " << o << std::endl;
		  }
		  (*expcoefM1[whch][iC])(n, o) += wk[n]*wk[o]/mass;
		}
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
	    if (not M0_only) {

	      double fac1 = facL*cosm[id][m];
	      double fac2 = facL*sinm[id][m];

	      for (int n=0; n<nmax; n++) {

		wk[n] = potd[id](l, n)*mass*fac0/sqnorm(l, n);

		(*expcoef0[id][loffset+moffset  ])[n] += wk[n]*fac1;
		(*expcoef0[id][loffset+moffset+1])[n] += wk[n]*fac2;
	      }

	      if (compute and pcavar) {
		pthread_mutex_lock(&cc_lock);
		for (int n=0; n<nmax; n++) {
		  if (std::fabs(wk[n]*facL) < 1.0e-40) {
		    std::cout << "Contribution is " << wk[n]*facL << " for "
			      << l << " " << m << " " << n << std::endl;
		  }
		  (*expcoefT1[whch][iC])[n] += wk[n]*std::complex<double>(fac1, fac2);
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
  if (play_back) {
    determine_coefficients_playback();
    if (play_cnew) determine_coefficients_particles();
  } else {
    determine_coefficients_particles();
  }
}

void SphericalBasis::determine_coefficients_playback(void)
{
  // Do we need new coefficients?
  if (tnow <= lastPlayTime) return;
  lastPlayTime = tnow;

  // Set coefficient matrix size (only do it once)
  if (expcoefP.size()==0) {
    expcoefP.resize((Lmax+1)*(Lmax+1));
    for (auto & v : expcoefP) v = std::make_shared<Eigen::VectorXd>(nmax);
  }

  if (coefMaster) {

    if (myid==0) {
      auto ret = playback->interpolate(tnow);

      // Get the matrix
      auto mat = std::get<0>(ret);

      // Get the error signal
      if (not std::get<1>(ret)) stop_signal = 1;

      //            +--------- Counter in real array (cosine and sine arrays
      //            |          are interleaved)
      //            |    +---- Counter in complex array (cosing and sine
      //            |    |     components are the real and imag parts)
      //            v    v
      for (int l=0, L=0, M=0; l<=Lmax; l++) {
	for (int m=0; m<=l; m++, M++) {
	  *expcoefP[L] = mat.row(M).real();
	  MPI_Bcast((*expcoefP[L++]).data(), nmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  if (m) {
	    *expcoefP[L] = mat.row(M).imag();
	    MPI_Bcast((*expcoefP[L++]).data(), nmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  }
	}
      }
    } else {
      for (int l=0, L=0; l<=Lmax; l++) {
	for (int m=0; m<=l; m++) {
	  MPI_Bcast((*expcoefP[L++]).data(), nmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  if (m) {
	    MPI_Bcast((*expcoefP[L++]).data(), nmax, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  }
	}
      }
    }
    
  } else {

    auto ret = playback->interpolate(tnow);

    // Get the matrix
    auto mat = std::get<0>(ret);

    // Get the error signal
    if (not std::get<1>(ret)) stop_signal = 1;

    for (int l=0, L=0, M=0; l<=Lmax; l++) {
      for (int m=0; m<=l; m++, M++) {
	*expcoefP[L++] = mat.row(M).real();
	if (m) {
	  *expcoefP[L++] = mat.row(M).imag();
	}
      }
    }
  }
}

void SphericalBasis::determine_coefficients_particles(void)
{
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = std::make_shared<nvTracer>("SphericalBasis::determine_coefficients");

  std::chrono::high_resolution_clock::time_point start0, start1, finish0, finish1;

  start0 = std::chrono::high_resolution_clock::now();

  // Return if we should leave the coefficients fixed
  //
  if (!self_consistent && !firstime_coef && !initializing) return;

  // No covarance by default
  //
  compute = false;

  // Subsample computation flag
  //
  if (nint) {
    // Computed only for mstep==0 and every nint steps
    if (mstep==0 and this_step % nint == 0) {
      compute = true;
      requestSubsample = true;
    }
    else {
      requestSubsample  = false;
      subsampleComputed = false;

    }
  }

  if (pcavar or pcaeof) {
    if (this_step >= npca0) {
      compute = (mstep == 0) && !( (this_step-npca0) % npca);
    }
  }

  int loffset, moffset, use1;

  if (compute) {

    if (massT.size() == 0) {	// Allocate storage for subsampling
      if (defSampT) sampT = defSampT;
      else          sampT = floor(sqrt(component->CurTotal()));
      massT    .resize(sampT, 0);
      massT1   .resize(sampT, 0);
      countT   .resize(sampT, 0);
      countT1  .resize(sampT, 0);
      
      expcoefT .resize(sampT);
      for (auto & t : expcoefT ) {
	t.resize((Lmax+1)*(Lmax+2)/2);
	for (auto & v : t) v = std::make_shared<Eigen::VectorXcd>(nmax);
      }
      
      expcoefT1.resize(sampT);
      for (auto & t : expcoefT1) {
	t.resize((Lmax+1)*(Lmax+2)/2);
	for (auto & v : t) v = std::make_shared<Eigen::VectorXcd>(nmax);
      }

      expcoefM .resize(sampT);
      for (auto & t : expcoefM ) {
	t.resize((Lmax+1)*(Lmax+2)/2);
	for (auto & v : t) v = std::make_shared<Eigen::MatrixXcd>(nmax, nmax);
      }
      
      expcoefM1.resize(sampT);
      for (auto & t : expcoefM1) {
	t.resize((Lmax+1)*(Lmax+2)/2);
	for (auto & v : t) v = std::make_shared<Eigen::MatrixXcd>(nmax, nmax);
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
	for (auto & v : countT1)   v = 0;
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

  // Swap interpolation arrays
  //
  auto p = expcoefL[mlevel];
  
  expcoefL[mlevel] = expcoefN[mlevel];
  expcoefN[mlevel] = p;
  
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

      requestSubsample  = false;
      subsampleComputed = true;
    }

    if (npca0 < std::numeric_limits<int>::max()) pca_hall(compute);
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
  if (play_back and not play_cnew) return;

  used   = 0;
  resetT = tnow;
}


void SphericalBasis::multistep_update_begin()
{
  if (play_back and not play_cnew) return;

				// Clear the update matricies
  for (int n=0; n<nthrds; n++) {
    for (int M=mfirst[mdrft]; M<=multistep; M++) differ1[n][M].setZero();
  }

#ifdef SPH_UPDATE_TABLE
  if ((not component->NoSwitch() or mdrft==1) or (this_step==0 and mdrft==0)) {
    occt.resize(multistep+1);
    for (auto & o : occt) {
      o.resize(multistep+1);
      std::fill(o.begin(), o.end(), 0);
    }
  }
#endif
}

void SphericalBasis::multistep_update_finish()
{
  if (play_back and not play_cnew) return;

				// Combine the update matricies
				// from all nodes
  unsigned sz = (multistep - mfirst[mdrft] + 1)*(Lmax+1)*(Lmax+1)*nmax;

				// Zero the buffer space
				//
  std::fill(pack.begin(), pack.begin() + sz, 0.0);

  // Pack the difference matrices
  //
  for (int M=mfirst[mdrft]; M<=multistep; M++) {

    unsigned offset0 = (M - mfirst[mdrft])*(Lmax+1)*(Lmax+1)*nmax;

    for (int l=0; l<=Lmax*(Lmax+2); l++) {

      unsigned offset1 = l*nmax;

      for (int n=0; n<nthrds; n++) {
	for (int ir=0; ir<nmax; ir++) {
	  pack[offset0+offset1+ir] += differ1[n][M](l, ir);
	}
      }
    }
  }

  MPI_Allreduce (pack.data(), unpack.data(), sz, 
		 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  // Update the local coefficients
  //
  for (int M=mfirst[mdrft]; M<=multistep; M++) {

    unsigned offset0 = (M - mfirst[mdrft])*(Lmax+1)*(Lmax+1)*nmax;

    for (int l=0; l<=Lmax*(Lmax+2); l++) {

      unsigned offset1 = l*nmax;

      for (int ir=0; ir<nmax; ir++)
	(*expcoefN[M][l])[ir] += unpack[offset0+offset1+ir];
    }
  }

  //  +--- Deep debugging
  //  |
  //  v
  if (false and myid==0) {
    std::string filename = runtag + ".differ_sph_" + component->name;
    std::ofstream out(filename, ios::app);

    if (out) {
      
      out << std::string(14+16*nmax, '-') << std::endl;
      out << "# T=" << tnow << " mstep=" << mstep
	  << " mfirst=" << mfirst[mdrft] << std::endl;
      
      for (int M=mfirst[mdrft]; M<=multistep; M++) {
	unsigned offset0 = (M - mfirst[mdrft])*(Lmax+1)*(Lmax+1)*nmax;
	for (int l=0, ll=0; l<=Lmax; l++) {
	  for (int m=0; m<=l; m++) {
	    unsigned offset1 = ll*nmax;
	    out << std::setw(3) << M << std::setw(3) << l
		<< std::setw(3) << m << std::setw(3) << 'c';
	    for (int ir=0; ir<nmax; ir++)
	      out << std::setw(16) << unpack[offset0+offset1+ir];
	    out << std::endl;
	    ll++;
	    if (m) {
	      offset1 = ll*nmax;
	      out << std::setw(3) << M << std::setw(3) << l
		  << std::setw(3) << m << std::setw(3) << 's';
	      for (int ir=0; ir<nmax; ir++)
		out << std::setw(16) << unpack[offset0+offset1+ir];
	      out << std::endl;
	      ll++;
	    }
	  }
	}
      }
      
      out << std::string(14+16*nmax, '-') << std::endl
	  << "---- Updated coefficients"  << std::endl
	  << std::string(14+16*nmax, '-') << std::endl;

      for (int M=mfirst[mdrft]; M<=multistep; M++) {
	for (int l=0, ll=0; l<=Lmax; l++) {
	  for (int m=0; m<=l; m++) {
	    out << std::setw(3) << M << std::setw(3) << l
		<< std::setw(3) << m << std::setw(3) << 'c';
	    for (int ir=0; ir<nmax; ir++)
	      out << std::setw(16) << (*expcoefN[M][ll])[ir];
	    out << std::endl;
	    ll++;
	    if (m) {
	      out << std::setw(3) << M << std::setw(3) << l
		  << std::setw(3) << m << std::setw(3) << 's';
	      for (int ir=0; ir<nmax; ir++)
		out << std::setw(16) << (*expcoef[ll])[ir];
	      out << std::endl;
	      ll++;
	    }
	  }
	}
      }
      out << std::string(14+16*nmax, '-') << std::endl;

    } else {
	std::cout << "Error opening test file <" << filename << "> at T=" << tnow
		  << std::endl;
    }
  }
  // END: deep debug

#ifdef SPH_UPDATE_TABLE
  occt_output();
#endif
}

void SphericalBasis::multistep_update(int from, int to, Component *c, int i, int id)
{
  if (play_back and not play_cnew) return;
  if (c->freeze(i)) return;

  double mass = c->Mass(i) * component->Adiabatic();

#ifdef SPH_UPDATE_TABLE
  occt[from][to]++;
#endif
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
    double fac0  = -4.0*M_PI;
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

	double facL = factorial(l, m)*legs[id](l, m);
	
	if (m==0) {
	  for (int n=0; n<nmax; n++) {
	    double val = potd[id](l, n)*facL*mass*fac0/sqnorm(l, n);
	    
	    differ1[id][from](loffset+moffset, n) -= val;
	    differ1[id][  to](loffset+moffset, n) += val;
	  }
	  moffset++;

	} else {
	  double fac1 = facL*cosm[id][m];
	  double fac2 = facL*sinm[id][m];

	  for (int n=0; n<nmax; n++) {
	    double val1 = potd[id](l, n)*fac1*mass*fac0/sqnorm(l, n);
	    double val2 = potd[id](l, n)*fac2*mass*fac0/sqnorm(l, n);

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
  if (play_back and not play_cnew) return;

#ifdef TMP_DEBUG
  Eigen::MatrixXd tmpcoef((Lmax+1)*(Lmax+1), nmax);
  for (int l=0; l<(Lmax+1)*(Lmax+1); l++) {
    for (int j=0; j<nmax; j++) tmpcoef(l, j) = (*expcoef[l])[j];
  }
#endif

				// Clean coefficient matrix
				// 
  for (int l=0; l<(Lmax+1)*(Lmax+1); l++) expcoef[l]->setZero();

				// For debugging only
				//
  double saveF = 0.0, saveG = 0.0;
    
				// Interpolate to get coefficients above
				// 
  for (int M=0; M<mfirst[mdrft]; M++) {
    
    double numer = static_cast<double>(mdrft            - dstepL[M][mdrft]);
    double denom = static_cast<double>(dstepN[M][mdrft] - dstepL[M][mdrft]);

    double b = numer/denom;	// Interpolation weights
    double a = 1.0 - b;

    for (int l=0; l<=Lmax*(Lmax+2); l++) {
      for (int n=0; n<nmax; n++) {
	(*expcoef[l])[n] += a*(*expcoefL[M][l])[n] + b*(*expcoefN[M][l])[n];
	if (l==0 and n==1) saveF += a*(*expcoefL[M][l])[n] + b*(*expcoefN[M][l])[n];
      }
    }
    
    //  +--- Deep debugging
    //  |
    //  v
    if (false and myid==0) {
      std::cout << std::left << std::fixed
		<< "SPH INTERP M=" << std::setw(2) << M
		<< " mstep=" << std::setw(3) << mstep
		<< " mdrft=" << std::setw(3) << mdrft
		<< " T=" << std::setw(16) << tnow
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
		<< " M="     << std::setw( 3) << M
		<< " mstep=" << std::setw( 3) << mstep 
		<< " mstep=" << std::setw( 3) << mdrft
		<< " minS="  << std::setw( 3) << dstepL[M][mdrft]
		<< " maxS="  << std::setw( 3) << dstepN[M][mdrft]
		<< " T="     << std::setw(12) << tnow
		<< " a="     << std::setw(12) << a 
		<< " b="     << std::setw(12) << b 
		<< " L01="   << std::setw(12) << (*expcoefL[M][0])[1]
		<< " N01="   << std::setw(12) << (*expcoefN[M][0])[1]
		<< " c01="   << std::setw(12) << (*expcoef[0])[1]
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
      for (int n=0; n<nmax; n++) {
	(*expcoef[l])[n] += (*expcoefN[M][l])[n];
	if (l==0 and n==1)saveG += (*expcoefN[M][l])[n];
      }
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
	      << " T="     << std::setw( 8) << tnow
	      << " c01="   << std::setw(14) << (*expcoef[0])[1]
	      << " u01="   << std::setw(14) << saveF
	      << " v01="   << std::setw(14) << saveG
	      << std::endl;
  }

#ifdef TMP_DEBUG
  if (myid==0) {

    static bool first = true;
    double maxval=0.0, maxdif=0.0, maxreldif=0.0, val;
    int irmax=1, lmax=0, irmaxdif=1, lmaxdif=0, irmaxrel=1, lmaxrel=0;
    for (int ir=0; ir<nmax; ir++) {
      for (int l=0; l<=Lmax*(Lmax+2); l++) {
	val = (*expcoef[l])[ir] - tmpcoef(l, ir);

	if (fabs((*expcoef[l])[ir]) > maxval) {
	  maxval = (*expcoef[l])[ir];
	  irmax = ir;
	  lmax = l;
	}
	if (fabs(val) > maxdif) {
	  maxdif = val;
	  irmaxdif = ir;
	  lmaxdif = l;
	}

	if (fabs(val/(*expcoef[l])[ir]) > maxreldif) {
	  maxreldif = val/(*expcoef[l])[ir];
	  irmaxrel = ir;
	  lmaxrel = l;
	}
      }
    }

    std::ofstream out(runtag+".coefs.diag", std::ios::app);
    if (first) {
      out << std::left << std::setw(15) << "# Time"
	  << std::left << std::setw(10) << "| mstep"
	  << std::left << std::setw(10) << "| mdrft"
	  << std::left << std::setw(18) << "| max(val)"
	  << std::left << std::setw(8)  << "| l"
	  << std::left << std::setw(8)  << "| n"
	  << std::left << std::setw(18) << "| max(dif)"
	  << std::left << std::setw(8)  << "| l"
	  << std::left << std::setw(8)  << "| n"
	  << std::left << std::setw(18) << "| max(reldif)"
	  << std::left << std::setw(8)  << "| l"
	  << std::left << std::setw(8)  << "| n"
	  << std::endl;
      first = false;
    }
    out << std::left << std::setw(15) << tnow
	<< std::left << std::setw(10) << mstep
	<< std::left << std::setw(10) << mdrft
	<< std::left << std::setw(18) << maxval
	<< std::left << std::setw(8)  << lmax
	<< std::left << std::setw(8)  << irmax
	<< std::left << std::setw(18) << maxdif
	<< std::left << std::setw(8)  << lmaxdif
	<< std::left << std::setw(8)  << irmaxdif
	<< std::left << std::setw(18) << maxreldif
	<< std::left << std::setw(8)  << lmaxrel
	<< std::left << std::setw(8)  << irmaxrel
	<< std::endl;
  }
#endif

#ifdef MULTI_DEBUG
  if (myid==0) {

    std::set<int> L = {0, 1};
    std::set<int> N = {0, 1};

    static bool first = true;

    std::ofstream out(runtag+".multi_diag", std::ios::app);

    if (first) {
      out << "#" << std::setw(9) << std::right << "Time "
	  << std::setw(4) << "l" << std::setw(4) << "n";
      for (int M=0; M<=multistep; M++)  out << std::setw(14) << M << " ";
      out << std::setw(14) << "Sum" << std::endl;
      first = false;
    }

    for (int l : L) {

      for (int n : N) {

	out << std::setw(10) << std::fixed << tnow
	    << std::setw(4) << l << std::setw(4) << n;

	double sum = 0.0, val;
	for (int M=0; M<mfirst[mdrft]; M++) {
      
	  double numer = static_cast<double>(mdrft            - dstepL[M][mdrft]);
	  double denom = static_cast<double>(dstepN[M][mdrft] - dstepL[M][mdrft]);

	  double b = numer/denom;	// Interpolation weights
	  double a = 1.0 - b;
      
	  val = a*(*expcoefL[M][l])[n] + b*(*expcoefN[M][l])[n];
	  sum += val;
	  out << std::setw(14) << val;
	}
	
	for (int M=mfirst[mdrft]; M<=multistep; M++) {
	  val = (*expcoefN[M][l])[n];
	  sum += val;
	  out << std::setw(14) << val;
	}

	out << std::setw(14) << sum << std::endl; // Next line

      }
      // END: N loop
    }
    // END: L loop
  }
  // END: MULTI_DEBUG
#endif
}

void * SphericalBasis::determine_acceleration_and_potential_thread(void * arg)
{
  // Density basis is 4piG times physical density
  //
  double dfac=0.25/M_PI;

  // Local variables
  //
  double r0=0.0, dp;
  double potr, potl, pott, potp, p, pc, dpc, ps, dps, facp, facdp;

  double pos[3];
  double xx, yy, zz, mfactor=1.0;

  vector<double> ctr;
  if (mix) mix->getCenter(ctr);

  int id = *((int*)arg);

  thread_timing_beg(id);

  // If we are multistepping, compute accel only at or above <mlevel>
  //
  for (int lev=mlevel; lev<=multistep; lev++) {

    unsigned nbodies = cC->levlist[lev].size();

    if (nbodies==0) continue;

    int nbeg = nbodies*(id  )/nthrds;
    int nend = nbodies*(id+1)/nthrds;

#ifdef DEBUG
    pthread_mutex_lock(&io_lock);
    std::cout << "Process " << myid << ": in thread"
	      << " id=" << id 
	      << " level=" << lev
	      << " nbeg=" << nbeg << " nend=" << nend << std::endl;
    pthread_mutex_unlock(&io_lock);
#endif

    for (int i=nbeg; i<nend; i++) {

      int indx = cC->levlist[lev][i];

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

      double r = sqrt(xx*xx + yy*yy + zz*zz) + DSMALL;
      double costh = zz/r;
      double rs = r/scale;
      double phi = atan2(yy, xx);

      dlegendre_R (Lmax, costh, legs[id], dlegs[id]);
      sinecosine_R(Lmax, phi,   cosm[id], sinm [id]);

      int ioff = 0;
      if (r>rmax) {
	ioff = 1;
	r0   = r;
	r    = rmax;
	rs   = r/scale;
      }

      // Zero coefficient accumulated field values
      //
      potl = potr = pott = potp = 0.0;
      
      get_dpotl(Lmax, nmax, rs, potd[id], dpot[id], id);

      if (!NO_L0) {
	get_pot_coefs_safe(0, *expcoef[0], p, dp, potd[id], dpot[id]);
	if (ioff) {
	  p *= rmax/r0;
	  dp = -p/r0;
	}
	double facL = mfactor * factorial(0, 0);
	potl = facL * p;
	potr = facL * dp;
      }
      
      //		l loop
      //		------
      for (int l=1, loffset=1; l<=Lmax; loffset+=(2*l+1), l++) {

				// Suppress L=1 terms?
	if (NO_L1 && l==1) continue;
	
				// Suppress odd L terms?
	if (EVEN_L && (l/2)*2 != l) continue;

	//		m loop
	//		------
	for (int m=0, moffset=0; m<=l; m++) {
	  
	  double facL = factorial(l, m) *  legs[id](l, m) * mfactor;
	  double facD = factorial(l, m) * dlegs[id](l, m) * mfactor;

				// Suppress odd M terms?
	  if (EVEN_M && (m/2)*2 != m) continue;

				// Suppress all asymmetric terms
	  if (M0_only and m!=0) continue;

	  if (m==0) {
	    get_pot_coefs_safe(l, *expcoef[loffset+moffset], p, dp,
			       potd[id], dpot[id]);
	    if (ioff) {
	      p *= pow(rmax/r0,(double)(l+1));
	      dp = -p/r0 * (l+1);
	    }
	    potl += facL * p;
	    potr += facL * dp;
	    pott += facD * p;
	    moffset++;
	  }
	  else {
	    get_pot_coefs_safe(l, *expcoef[loffset+moffset], pc, dpc,
			       potd[id], dpot[id]);

	    get_pot_coefs_safe(l, *expcoef[loffset+moffset+1], ps, dps,
			       potd[id], dpot[id]);
	    if (ioff) {		// Factors for external multipole solution
	      facp  = pow(rmax/r0,(double)(l+1));
	      facdp = -1.0/r0 * (l+1);
				// Apply the factors
	      pc   *= facp;
	      ps   *= facp;
	      dpc   = pc * facdp;
	      dps   = ps * facdp;
	    }
	    potl += facL * (pc *cosm[id][m] + ps *sinm[id][m] );
	    potr += facL * (dpc*cosm[id][m] + dps*sinm[id][m] );
	    pott += facD * (pc *cosm[id][m] + ps *sinm[id][m] );
	    potp += facL * (-pc*sinm[id][m] + ps *cosm[id][m] )*m;
	    moffset +=2;
	  }
	}
      }

      double fac = xx*xx + yy*yy;

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

  if (play_back) {
    swap_coefs(expcoefP, expcoef);
  }

  if (use_external == false) {

    if (multistep && (self_consistent || initializing)) {
      compute_multistep_coefficients();
    }

  }

  // This should suffice for both CPU and GPU
  if (FIX_L0) {
    // Save the monopole coefficients on the first evaluation
    if (C0.size() == 0) C0 = *expcoef[0];
    // Copy the saved coefficients to the active array
    else *expcoef[0] = C0;
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

  if (play_back) {
    swap_coefs(expcoef, expcoefP);
  }


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

  p  = pp;
  dp = dpp;
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

  p  = pp;
  dp = dpp;
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
	for (int m=0, moffset=0; m<=l; m++) {
	  if (m==0) {
	    out.write((char *)&(z=(*expcoef[loffset+moffset+0])[ir]), sizeof(double));
	    moffset += 1;
	  } else {
	    out.write((char *)&(z=(*expcoef[loffset+moffset+0])[ir]), sizeof(double));
	    out.write((char *)&(z=(*expcoef[loffset+moffset+1])[ir]), sizeof(double));
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

// Dump coefficients to an HDF5 file

void SphericalBasis::dump_coefs_h5(const std::string& file)
{
  // Add the current coefficients
  auto cur = std::make_shared<CoefClasses::SphStruct>();

  cur->time   = tnow;
  cur->geom   = geoname[geometry];
  cur->id     = id;
  cur->time   = tnow;
  cur->lmax   = Lmax;
  cur->nmax   = nmax;
  cur->scale  = scale;
  cur->normed = true;

  cur->allocate();

  auto & cof = *cur->coefs;	// Reference for convenience

  for (int ir=0; ir<nmax; ir++) {
    for (int l=0, L=0, offset=0; l<=Lmax; l++) {
      for (int m=0; m<=l; m++, L++) {
	if (m==0) {
	  cof(L, ir) = {(*expcoef[offset])[ir], 0.0};
	  offset += 1;
	} else {
	  cof(L, ir) = {(*expcoef[offset])[ir], (*expcoef[offset+1])[ir]};
	  offset += 2;
	}
      }
    }
  }

  // Add center
  //
  cur->ctr = component->getCenter(Component::Local | Component::Centered);

  // Add the orientation
  cur->rot = component->getRotation();

  // Check if file exists
  //
  if (std::filesystem::exists(file)) {
    sphCoefs.clear();
    sphCoefs.add(cur);
    sphCoefs.ExtendH5Coefs(file);
  }
  // Otherwise, extend the existing HDF5 file
  //
  else {
    // Copy the YAML config.  We only need this on the first call.
    std::ostringstream sout; sout << conf;
    cur->buf = sout.str();	// Copy to CoefStruct buffer

    // Add the name attribute.  We only need this on the first call.
    sphCoefs.setName(component->name);

    // Add the default units
    sphCoefs.setUnits({{"length", "none", 1.0},
		       {"mass",   "none", 1.0},
		       {"time",   "none", 1.0}});

    // And the new coefficients and write the new HDF5
    sphCoefs.clear();
    sphCoefs.add(cur);
    sphCoefs.WriteH5Coefs(file);
  }
}


void SphericalBasis::determine_fields_at_point
(double x, double y, double z, 
 double *tdens0, double *tpotl0, 
 double *tdens,  double *tpotl, 
 double *tpotX,  double *tpotY, double *tpotZ)
{
  double R2 = x*x + y*y;
  double R  = sqrt(R) + DSMALL;
  double r  = sqrt(R2 + z*z) + DSMALL;
  double r2 = R2 + z*z + DSMALL;
  double r3 = r2*r;
  
  double theta = acos(z/r);
  double phi   = atan2(y, x);
  double cth   = cos(theta), sth = sin(theta);
  double cph   = cos(phi),   sph = sin(phi);
  double tpotr, tpott, tpotp;

  determine_fields_at_point_sph(r, theta, phi,
				tdens0, tpotl0, 
				tdens, tpotl, 
				&tpotr,	&tpott, &tpotp);

  *tpotX = tpotr*x/r - tpott*x*z/r3;
  *tpotY = tpotr*y/r - tpott*y*z/r3;
  *tpotZ = tpotr*z/r + tpott*R2/r3;
      
  if (R > DSMALL) {
    *tpotX +=  tpotp*y/R;
    *tpotY += -tpotp*x/R;
  }
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

  *tpotR = tpotr*sin(theta) - tpott*cos(theta)/r;
  *tpotz = tpotr*cos(theta) + tpott*sin(theta)/r;
}

void SphericalBasis::determine_fields_at_point_sph
(double r, double theta, double phi, 
 double *tdens0, double *tpotl0, 
 double *tdens, double *tpotl, 
 double *tpotr, double *tpott, double *tpotp)
{
  double dp, potr, potl, pott, potp, p, pc, dpc, ps, dps, dens;
  double dfac=0.25/M_PI;

  double rs = r/scale;
  double costh = cos(theta);

  dlegendre_R(Lmax, costh, legs[0], dlegs[0]);
  sinecosine_R(Lmax, phi, cosm[0], sinm[0]);

  double facL = factorial(0, 0) * legs[0](0, 0);

  get_dens(Lmax, nmax, rs, dend, 0);
  get_dpotl(Lmax, nmax, rs, potd[0], dpot[0], 0);
  get_dens_coefs(0, *expcoef[0], dens);
  dens *= dfac * facL;

  get_pot_coefs(0, *expcoef[0], p, dp);
  potl = facL * p;
  potr = facL * dp;
  pott = potp = 0.0;
  
  *tdens0 = dens;
  *tpotl0 = potl;

  // l loop
  //
  for (int l=1, loffset=1; l<=Lmax; loffset+=(2*l+1), l++) {
    
    // m loop
    //
    for (int m=0, moffset=0; m<=l; m++) {

      double facL = factorial(l, m) * legs[0](l, m);
      double facD = factorial(l, m) * dlegs[0](l, m);

      if (m==0) {
	get_dens_coefs(l, *expcoef[loffset+moffset], p);
	dens += dfac*facL*p;
	get_pot_coefs(l, *expcoef[loffset+moffset], p, dp);
	potl += facL * p;
	potr += facL * dp;
	pott += facD * p;
	moffset++;
      }
      else {
	get_dens_coefs(l, *expcoef[loffset+moffset], pc);
	get_dens_coefs(l, *expcoef[loffset+moffset+1], ps);
	dens += dfac * facL * (pc*cosm[0][m] + ps*sinm[0][m]);
	
	get_pot_coefs(l, *expcoef[loffset+moffset], pc, dpc);
	get_pot_coefs(l, *expcoef[loffset+moffset+1], ps, dps);
	potl += facL * (pc*cosm[0][m] + ps*sinm[0][m]);
	potr += facL * (dpc*cosm[0][m] + dps*sinm[0][m]);
	pott += facD * (pc*cosm[0][m] + ps*sinm[0][m]);
	potp += facL * (-pc*sinm[0][m] + ps*cosm[0][m])*m;
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

  for (int i=0; i<numg; i++) {
    double r = rmin + del*qe.knot(i);
    double rs = r / scale;

    get_potl(Lmax, nmax, rs, potd[0], 0);

    for(int l=0; l<=Lmax; l++) {
      
      for (int n=0; n<nmax; n++) {

	double pot = potd[0](l, n)/sqnorm(l, n)/scale;

	if (l==0)
	  meanC[n] += del * qe.weight(i) * r * r * pot *
	    4.0*M_PI * modl.get_density(r);

	rmsC(l, n) += del * qe.weight(i) * r * r * pot * pot *
	  4.0*M_PI * modl.get_density(r);
      }
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

	double fac = sqrt((2.0*l + 1.0)/(4.0*M_PI));

	for (int n=0; n<nmax; n++) {
	  if (l==0)
	    out << setw(5)  << n
		<< setw(16) << (*expcoef[l])[n]
		<< setw(16) << meanC[n] * fac
		<< setw(16) << rmsC(l, n) * fac*fac
		<< setw(16) << (rmsC(l, n) - meanC[n]*meanC[n])*fac*fac
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
	    sqrt(fabs(rmsC(l, n) - meanC[n]*meanC[n])*factorial(l, m)/noiseN)*nrand(rgen);
	  if (l==0) (*expcoef[l])[n] += meanC[n];
	}
	moffset++;
      }
      else {
	for (int n=0; n<nmax; n++) {
	  (*expcoef[loffset+moffset+0])[n] = sqrt(fabs(rmsC(l, n) - meanC[n]*meanC[n])*factorial(l, m)/noiseN)*nrand(rgen);
	  (*expcoef[loffset+moffset+1])[n] = 
	    sqrt(fabs(rmsC(l, n) - meanC[n]*meanC[n])*factorial(l, m)/noiseN)*nrand(rgen);
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

void SphericalBasis::occt_output()
{
  for (int m=0; m<=multistep; m++) {
    if (myid) {
      MPI_Reduce(occt[m].data(), NULL, multistep+1,
		 MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(MPI_IN_PLACE, occt[m].data(), multistep+1,
		 MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  }
    
  if (myid==0) {

    unsigned total = 0;
    for (int m0=0; m0<=multistep; m0++) {
      for (int m1=0; m1<=multistep; m1++) total += occt[m0][m1];
    }
    
    if (total) {
      std::cout << std::string(8*(multistep+1), '-') << std::endl
		<< "SphericalBasis update matrix at T=" << tnow
		<< " mdrft=" << mdrft << std::endl
		<< std::string(8*(multistep+1), '-') << std::endl;
      
      for (int m0=0; m0<=multistep; m0++) {
	for (int m1=0; m1<=multistep; m1++) {
	  std::cout << std::setw(8) << occt[m0][m1];
	}
	std::cout << std::endl;
      }
      std::cout << std::string(8*(multistep+1), '-') << std::endl;
    }
  }
}

void SphericalBasis::biorthogonality_check()
{
  // Allocate storage:
  // - 1st index is the node
  // - 2nd index is L
  //
  std::vector<std::vector<Eigen::MatrixXd>> one(numprocs);
  for (auto & v : one) {
    v.resize(Lmax+1);
    for (auto & u : v) {
      u.resize(nmax, nmax);
      u.setZero();
    }
  }

  // Number of knots
  //
  const int num = 400;
  LegeQuad wk(num);
  
  // Mapped radial range
  //
  double xmin = r_to_xi(rmin), xmax = r_to_xi(rmax);
  double dx = xmax - xmin;

  Eigen::MatrixXd p(Lmax+1, nmax), d(Lmax+1, nmax);


  // Biorthogonal integral loop
  //
  for (int i=0; i<num; i++) {

    // Each node contributes individual radii
    //
    if (i % numprocs == myid) {

      double x = xmin + dx*wk.knot(i);
      double r = xi_to_r(x);
      double w = dx*wk.weight(i) * d_r_to_xi(r) * r * r;

      // Evaluate basis at radius r
      //
      get_potl(Lmax, nmax, r, p, 0);
      get_dens(Lmax, nmax, r, d, 0);

      // Contribution to the integrand
      //
      for (int L=0, cnt=0; L<=Lmax; L++) {
	for (int n1=0; n1<nmax; n1++) {
	  for (int n2=0; n2<nmax; n2++, cnt++) {
	    one[myid][L](n1, n2) += w * p(L, n1) * d(L, n2);
	  }
	}
      }
    }
  }

  // Reduce from workers to root and print to file
  //
  if (myid==0) {

    // Get contributions from worker nodes
    //
    for (int L=0; L<=Lmax; L++)
      MPI_Reduce(MPI_IN_PLACE, one[0][L].data(), one[0][L].size(),
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Output file
    //
    std::ostringstream filename;
    filename << "ortho_check." << component->name << "." << runtag;
    std::ofstream out(filename.str());
    if (out) {
      // Header
      //
      out << "# component=" << component->name << std::endl
	  << "# force ID =" << this->id << std::endl
	  << "#" << std::endl;

      // Print biorthogonality matrices
      //
      for (int L=0, cnt=0; L<=Lmax; L++) {
	out << "#" << std::string(72, '-') << std::endl;
	out << "# L=" << L << std::endl;
	out << "#" << std::string(72, '-') << std::endl;
	out << std::scientific << std::setprecision(4)
	    << one[0][L] << std::endl;
      }
    }
  }
  else {
    // Send contributions to root node
    //
    for (int L=0; L<=Lmax; L++)
      MPI_Reduce(one[myid][L].data(), 0, one[myid][L].size(),
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
}

PotAccel::CovarData SphericalBasis::getSubsample()
{
  using covTuple = std::tuple<Eigen::VectorXcd,
			      Eigen::MatrixXcd>;
  // Prepare the covariance structure
  PotAccel::CovarData ret;

  // Number of real angular terms l in [0, Lmax], m in [0, l]
  int totL = (Lmax+1)*(Lmax+2)/2;

  // Sanity check
  if (expcoefT[0].size() != totL) {
    std::ostringstream ss;
    ss << "SphericalBasis::getSubsample() "
       << "internal error: unexpected number of absolute angular terms, "
       << expcoefT[0].size() << " != " << totL;
    throw std::runtime_error(ss.str());
  }

  std::get<0>(ret).resize(sampT);	// Sample counts
  std::get<1>(ret).resize(sampT);	// Sample masses
  std::get<2>(ret) = Eigen::Tensor<std::complex<double>, 3>(sampT, totL, nmax);

  if (fullCovar)
    std::get<3>(ret) = Eigen::Tensor<std::complex<double>, 4>(sampT, totL, nmax, nmax);

  // Fill the covariance structure with subsamples
  for (int T=0; T<sampT; T++) {

    std::get<0>(ret)(T) = countT[T];
    std::get<1>(ret)(T) =  massT[T];

    // l loop
    for (int l=0, k=0; l<=Lmax; l++) {
      // m loop
      for (int m=0; m<=l; m++) {
	// outer n loop
	for (int n1=0; n1<nmax; n1++) {
	  auto z = (*expcoefT[T][k])(n1);
	  if (isnan(z.real()) or isnan(z.imag())) {
	    std::cout << "Pack coef NaN at " << l << ", " << m << ", "
		      << n1 << std::endl;
	  }
	  if (std::abs(z.real()) < 1.0e-12) {
	    std::cout << "Pack coef zero at " << l << ", " << m << ", "
		      << n1 << std::endl;
	  }
	  std::get<2>(ret)(T, k, n1) = (*expcoefT[T][k])(n1);
	  // inner n loop
	  for (int n2=0; n2<nmax; n2++) {
	    auto z = (*expcoefM[T][k])(n1, n2);
	    if (isnan(z.real()) or isnan(z.imag())) {
	      std::cout << "Pack covar NaN at " << l << ", " << m << ", "
			<< n1 << ", " << n1 << std::endl;
	    }
	    if (std::fabs(z) < 1.0e-20) {
	      std::cout << "Pack covar zero at " << l << ", " << m << ", "
			<< n1 << ", " << n2 << std::endl;
	    }
	    if (fullCovar)
	      std::get<3>(ret)(T, k, n1, n2) = (*expcoefM[T][k])(n1, n2);
	  }
	  // END: n2 loop
	}
	// END: n1 loop
	k++;
      }
      // END: m loop
    }
    // END: l loop
  }
  // END: T loop
    
  return ret;
}


void SphericalBasis::writeCovarH5Params(HighFive::File& file)
{
  file.createAttribute<int>("lmax", HighFive::DataSpace::From(Lmax)).write(Lmax);
  file.createAttribute<int>("nmax", HighFive::DataSpace::From(nmax)).write(nmax);
  file.createAttribute<double>("scale", HighFive::DataSpace::From(scale)).write(scale);
  file.createAttribute<double>("rmin", HighFive::DataSpace::From(rmin)).write(rmin);
  file.createAttribute<double>("rmax", HighFive::DataSpace::From(rmax)).write(rmax);
}

