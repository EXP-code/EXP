#include "expand.H"

#include <filesystem>
#include <cstdlib>
#include <cmath>

#include "Cube.H"

const std::set<std::string>
Cube::valid_keys = {
  "nminx",
  "nminy",
  "nminz",
  "nmaxx",
  "nmaxy",
  "nmaxz",
  "method",
  "wrap",
  "nint",
  "samplesz",
  "subsampleFloat"
};

//@{
//! These are for testing exclusively (should be set false for production)
static bool cudaAccumOverride = false;
static bool cudaAccelOverride = false;
static bool deepDebug = false;
static bool coefDebug = false;
//@}

Cube::Cube(Component* c0, const YAML::Node& conf) : PotAccel(c0, conf)
{
  // ID parameters
  //
  id         = "Cube";
  geometry   = cube;
  coef_dump  = true;
  byPlanes   = true;
  cuMethod   = "planes";
  sampT      = 100;

  // Default parameter values
  //
  nminx = nminy = nminz = 0;
  nmaxx = nmaxy = nmaxz = 16;

#if HAVE_LIBCUDA==1
  cuda_aware = true;
#endif

  // Parse parameters
  //
  initialize();

  // Cache computed paramters
  //
  imx   = 1 + 2*nmaxx;		// number of x wave numbers
  imy   = 1 + 2*nmaxy;		// number of x wave numbers
  imz   = 1 + 2*nmaxz;		// number of x wave numbers
  osize = imx * imy * imz;	// total number of coefficients

  // Allocate storage
  //
  expcoef.resize(imx, imy, imz);
  expcoef0.resize(nthrds);
  for (auto & v : expcoef0) v.resize(imx, imy, imz);

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
  unsigned sz = (multistep+1)*osize;
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
    
  // Constant factors
  //
  dfac = 2.0*M_PI;
  kfac = std::complex<double>(0.0, dfac);

  // Initialize covariance
  //
  init_covariance();

#if HAVE_LIBCUDA==1
  cuda_initialize();
#endif

}

Cube::~Cube(void)
{
#if HAVE_LIBCUDA==1
  if (component->cudaDevice>=0) destroy_cuda();
#endif
}

void Cube::initialize(void)
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["nminx" ])  nminx      = conf["nminx" ].as<int>();
    if (conf["nminy" ])  nminy      = conf["nminy" ].as<int>();
    if (conf["nminz" ])  nminz      = conf["nminz" ].as<int>();
    if (conf["nmaxx" ])  nmaxx      = conf["nmaxx" ].as<int>();
    if (conf["nmaxy" ])  nmaxy      = conf["nmaxy" ].as<int>();
    if (conf["nmaxz" ])  nmaxz      = conf["nmaxz" ].as<int>();
    if (conf["method"])  cuMethod   = conf["method"].as<std::string>();
    if (conf["wrap"  ])  wrap       = conf["wrap"  ].as<bool>();

    if (conf["nint"]) {
      nint = conf["nint"].as<int>();
      if (nint>0) computeSubsample = true;
    }

    if (conf["samplesz"]) {
      sampT = conf["samplesz"].as<int>();
    }

  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Cube: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("Cube::initialize: error parsing YAML");
  }

  if (myid==0)
    std::cout << "---- Cube::initialize: wrap="
	      << std::boolalpha << wrap << std::endl;
}

void Cube::init_covariance()
{
  if (computeSubsample) {

    meanV.resize(sampT);
    for (auto& v : meanV) {
      v.resize(osize);
    }

    workV1.resize(nthrds);
    for (auto& v : workV1) v.resize(osize);

    if (fullCovar) {
      covrV.resize(sampT);
      for (auto& v : covrV) {
	v.resize(osize, osize);
      }
    } else {
      covrV.clear();
    }

    sampleCounts.resize(sampT);
    sampleMasses.resize(sampT);
      
    meanV1.resize(nthrds);
    covrV1.resize(nthrds);
    countV1.resize(nthrds);
    massV1.resize(nthrds);

    for (int n=0; n<nthrds; n++) {
      meanV1[n].resize(sampT);
      covrV1[n].resize(sampT);
      for (int T=0; T<sampT; T++) {
	meanV1[n][T].resize(osize);
	if (fullCovar) {
	  covrV1[n][T].resize(osize, osize);
	}
      }
      countV1[n].resize(sampT);
      massV1[n].resize(sampT);
    }

    zero_covariance();
  }
}


void Cube::zero_covariance()
{
  for (int T=0; T<sampT; T++) {
    meanV[T].setZero();
    if (fullCovar) {
      covrV[T].setZero();
    }
  }
    
  sampleCounts.setZero();
  sampleMasses.setZero();

  for (int n=0; n<nthrds; n++) {
    for (int T=0; T<sampT; T++) {
      meanV1[n][T].setZero();
      if (fullCovar) {
	covrV1[n][T].setZero();
      }
    }
    workV1[n].setZero();
    countV1[n].setZero();
    massV1[n].setZero();
  }
}


void * Cube::determine_coefficients_thread(void * arg)
{
  int id = *((int*)arg);
  double adb = component->Adiabatic();

  // Number of bodies at this level
  //
  unsigned nbodies = cC->levlist[mlevel].size();

  // Skip empty levels
  //
  if (nbodies) {

    // Partition bodies by thread
    int nbeg = nbodies*(id  )/nthrds;
    int nend = nbodies*(id+1)/nthrds;

    for (int q=nbeg; q<nend; q++) {

      int i = cC->levlist[mlevel][q];

      use[id]++;
      double mass = cC->Mass(i) * adb;

      // Get position
      //
      double x = cC->Pos(i, 0);
      double y = cC->Pos(i, 1);
      double z = cC->Pos(i, 2);

      // Truncate to cube with sides in [0,1]
      //
      if (wrap) {
	auto unitCube = [](double x) {
	  if (x<0.0) x += std::floor(-x) + 1.0;
	  else       x -= std::floor( x);
	  return x;
	};
	
	x = unitCube(x);
	y = unitCube(y);
	z = unitCube(z);
      }

      // Only compute for points inside the unit cube
      //
      if (x<0.0 or x>1.0) continue;
      if (y<0.0 or y>1.0) continue;
      if (z<0.0 or z>1.0) continue;
    
      // Recursion multipliers
      //
      std::complex<double> stepx = std::exp(-kfac*x);
      std::complex<double> stepy = std::exp(-kfac*y);
      std::complex<double> stepz = std::exp(-kfac*z);
      
      // Initial values for recursion
      //
      std::complex<double> startx = std::exp(kfac*(x*nmaxx));
      std::complex<double> starty = std::exp(kfac*(y*nmaxy));
      std::complex<double> startz = std::exp(kfac*(z*nmaxz));
      
      std::complex<double> facx, facy, facz;
      int ix, iy, iz;

      if (requestSubsample) workV1[id].setZero();

      for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
	for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	  for (facz=startz, iz=0; iz<imz; iz++, facz*=stepz) {
	    
	    // Compute wavenumber; the coefficients are stored as:
	    // -nmax,-nmax+1,...,0,...,nmax-1,nmax
	    //
	    int ii = ix - nmaxx;
	    int jj = iy - nmaxy;
	    int kk = iz - nmaxz;
	    
	    if (ii==0 and jj==0 and kk==0) continue;
	    
	    // Normalization
	    double norm = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));
	    
	    expcoef0[id](ix, iy, iz) += - mass * facx * facy * facz * norm;

	    if (requestSubsample) {
	      workV1[id]( ( (ii + nmaxx)*imy + (jj + nmaxy) )*imz + (kk + nmaxz) )
		= facx * facy * facz * norm;
	    }

	    if (deepDebug and ii==1 and jj==0 and kk==0) {
	      auto part = cC->Part(i);
	      if (part->indx < 10) {
		std::complex<double> tst = - mass * facx * facy * facz * norm;
		std::cout << "coef contrib: [" << std::setw(8) << tnow << "] "
			  << std::setw( 6) << part->indx
			  << std::setw(16) << part->pos[0]
			  << std::setw(16) << part->pos[1]
			  << std::setw(16) << part->pos[2]
			  << std::setw(16) << facx.real()
			  << std::setw(16) << facx.imag()
			  << std::setw(16) << facy.real()
			  << std::setw(16) << facy.imag()
			  << std::setw(16) << facz.real()
			  << std::setw(16) << facz.imag()
			  << std::setw(16) << norm * mass
			  << std::setw(16) << tst.real()
			  << std::setw(16) << tst.imag()
			  << std::endl;
	      }
	    }
	    // END: deepDebug
	  }
	  // END: iz loop
	}
	// END: iy loop
      }
      // END: ix loop

      if (requestSubsample) {
	// Which subsample bin?
	//
	int T = q % sampT;

	// Accumulate counts and masses
	//
	countV1[id](T) += 1;
	massV1[id](T)  += mass;

	// Accumulate subsample contributions
	//
	meanV1[id][T] += workV1[id] * mass;
	if (fullCovar)
	  covrV1[id][T] += workV1[id] * workV1[id].adjoint() * mass;
      }

    }
    // END: particle loop
  }
  // END: bodies at this level
  
  return (NULL);
}

void Cube::determine_coefficients(void)
{
  exeTimer timer(this, "Coefficient evaluation");

  //  Coefficients in the Eigen::Tensor are packed as
  //  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax in a single array for each
  //  dimension with z dimension changing most rapidly

  // Zero the coefficients
  //
  for (auto & v : expcoef0) v.setZero();

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
    
  use1 = 0;
  if (multistep==0) used = 0;
  if (mlevel==0) {
    for (int n=0; n<nthrds; n++) use[n] = 0.0;
  }
  
  // Determine whether or not to compute a subsample
  if (mstep==0 or mstep==std::numeric_limits<int>::max()) {
    if (nint>0 && this_step % nint == 0) {
      if (tnow > last) {
	requestSubsample = true;
	last = tnow;
	zero_covariance();
      }
    }
  } else {
    subsampleComputed = false;
  }
  
#if HAVE_LIBCUDA==1
  (*barrier)("Cube::entering cuda coefficients", __FILE__, __LINE__);
  if (component->cudaDevice>=0 and use_cuda) {
    if (cudaAccumOverride) {
      component->CudaToParticles();
      exp_thread_fork(true);
    } else {
      timer.Start1();
      determine_coefficients_cuda();
      DtoH_coefs(mlevel);
      timer.Stop1();
    }
  } else {
    exp_thread_fork(true);
  }
  (*barrier)("Cube::exiting cuda coefficients", __FILE__, __LINE__);
#else
  exp_thread_fork(true);
#endif

  for (int i=0; i<nthrds; i++) use1 += use[i];
  
  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  used = use0;

  for (int i=1; i<nthrds; i++) expcoef0[0] += expcoef0[i];
  
  if (multistep) {

    MPI_Allreduce( expcoef0[0].data(), expcoefN[mlevel]->data(),
		   expcoef0[0].size(),
		   MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  } else {
    
    MPI_Allreduce( expcoef0[0].data(), expcoef.data(), expcoef0[0].size(),
		   MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  }

  // Last level?
  //
  if (multistep and mlevel==multistep) {
      compute_multistep_coefficients();
  }

  // Accumulate mean and covariance subsample contributions
  //
  if (requestSubsample) {

    // Only finalize at the last multistep level
    //
    if ( (multistep and mlevel==multistep) or multistep==0 ) {
      
      // Sum over threads
      //
      for (int n=1; n<nthrds; n++) {
	for (int T=0; T<sampT; T++) {
	  meanV1[0][T] += meanV1[n][T];
	  if (fullCovar) {
	    covrV1[0][T] += covrV1[n][T];
	  }
	  countV1[0](T) += countV1[n](T);
	  massV1[0](T)  += massV1[n](T);
	}
      }

      // Sum over MPI ranks
      //
      for (int T=0; T<sampT; T++) {
	MPI_Allreduce( meanV1[0][T].data(), meanV[T].data(), meanV[T].size(),
		       MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
	if (fullCovar)
	  MPI_Allreduce( covrV1[0][T].data(), covrV[T].data(), covrV[T].size(),
			 MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
      }

      MPI_Allreduce( countV1[0].data(), sampleCounts.data(), sampleCounts.size(),
		     MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
      MPI_Allreduce( massV1[0].data(), sampleMasses.data(), sampleMasses.size(),
		     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      requestSubsample  = false;
      subsampleComputed = true;
    }
  }
    
  // Deep debug for checking a single wave number from cubeics
  //
  if (coefDebug and myid==0) {

    // Create a wavenumber tuple from a flattened index
    auto indices = [&](int indx)
    {
      int NX = 2*this->nmaxx+1, NY = 2*this->nmaxy+1;
      int k  = indx/(NX*NY);
      int j  = (indx - k*NX*NY)/NX;
      int i  = indx - (j + k*NY)*NX;
      
      return std::tuple<int, int, int>{i, j, k};
    };

    if (multistep==0) {

      std::string ofile = "cube_test." + runtag + ".dat";
      std::ofstream out(ofile, ios::app | ios::out);

      if (out) {
	std::multimap<double, int> biggest;

	for (int n=0; n<expcoef.size(); n++) 
	  biggest.insert({std::abs(expcoef.data()[n]), n});

	out << std::string(3*4+3*20, '-') << std::endl
	    << "---- Cube, T=" << tnow    << std::endl
	    << std::string(3*4+3*20, '-') << std::endl
	    << std::setprecision(10);
	
	out << std::setw(4)  << "i"
	    << std::setw(4)  << "j"
	    << std::setw(4)  << "k"
	    << std::setw(20) << "Real"
	    << std::setw(20) << "Imag"
	    << std::setw(20) << "Abs"
	    << std::endl;
	
	int cnt = 0;
	for (auto it = biggest.rbegin(); it!=biggest.rend() and cnt<20; it++, cnt++) {
	  auto [i, j, k] = indices(it->second);
	  auto a = expcoef(i, j, k);
	  out << std::setw(4)  << i-nmaxx
	      << std::setw(4)  << j-nmaxy
	      << std::setw(4)  << k-nmaxz
	      << std::setw(20) << std::real(a)
	      << std::setw(20) << std::imag(a)
	      << std::setw(20) << std::abs(a)
	      << std::endl;
	}
	out << std::string(3*4+4*20, '-') << std::endl;
      } else {
	std::cout << "Error opening <" << ofile << ">" << std::endl;
      }

    } else {

      std::string ofile = "cube_multi." + runtag + ".dat";
      std::ofstream out(ofile, ios::app | ios::out);

      if (out) {
	// Rank coefficients by absolute value
	std::multimap<double, int> biggest;

	// Make the DB of absolute values
	for (int n=0; n<expcoef.size(); n++) 
	  biggest.insert({std::abs(expcoef.data()[n]), n});

	out << std::string(3*4+3*(2+multistep)*20, '-') << std::endl
	    << "---- Cube, T=" << tnow    << std::endl
	    << std::string(3*4+3*20, '-') << std::endl
	    << std::setprecision(10);
	
	out << std::setw(4)  << "i"
	    << std::setw(4)  << "j"
	    << std::setw(4)  << "k"
	    << std::setw(20) << "Real"
	    << std::setw(20) << "Imag"
	    << std::setw(20) << "Abs";

	for (int M=0; M<=multistep; M++) {
	  std::ostringstream ss;  ss << "Real[" << M << "]";
	  out << std::setw(20) << ss.str();
	  ss.str(""); ss << "Imag[" << M << "]";
	  out << std::setw(20) << ss.str();
	  ss.str(""); ss << "Abs["  << M << "]";
	  out << std::setw(20) << ss.str();
	}
	out << std::endl;
	
	int cnt = 0;
	for (auto it = biggest.rbegin(); it!=biggest.rend() and cnt<20; it++, cnt++) {
	  auto [i, j, k] = indices(it->second);
	  auto a = expcoef(i, j, k);
	  out << std::setw(4)  << i-nmaxx
	      << std::setw(4)  << j-nmaxy
	      << std::setw(4)  << k-nmaxz
	      << std::setw(20) << std::real(a)
	      << std::setw(20) << std::imag(a)
	      << std::setw(20) << std::abs(a);
	  for (int M=0; M<=multistep; M++) {
	    auto b = (*expcoefN[M])(i, j, k);
	    out << std::setw(20) << std::real(b)
		<< std::setw(20) << std::imag(b)
		<< std::setw(20) << std::abs(b);
	  }
	  out << std::endl;
	}
	out << std::string(3*4+3*(2+multistep)*20, '-') << std::endl;

      } else {
	std::cout << "Error opening <" << ofile << ">" << std::endl;
      }

    }
  }
  // END: deep debug
}

void * Cube::determine_acceleration_and_potential_thread(void * arg)
{
  int id = *((int*)arg);

  // If we are multistepping, compute accel only at or above <mlevel>
  //
  for (int lev=mlevel; lev<=multistep; lev++) {

    unsigned nbodies = cC->levlist[lev].size();

    if (nbodies==0) continue;

    int nbeg = nbodies*id/nthrds;
    int nend = nbodies*(id+1)/nthrds;

    for (int q=nbeg; q<nend; q++) {
    
      int i = cC->levlist[lev][q];

      // Local accumulators 
      //
      std::complex<double> accx(0), accy(0), accz(0), dens(0), potl(0);
    
      // Get positions
      //
      double x = cC->Pos(i, 0);
      double y = cC->Pos(i, 1);
      double z = cC->Pos(i, 2);

      // Recursion multipliers
      //
      auto stepx = std::exp(kfac*x);
      auto stepy = std::exp(kfac*y);
      auto stepz = std::exp(kfac*z);
    
      // Initial values (note sign change from coefficient accumulation)
      //
      auto startx = std::exp(-kfac*(x*nmaxx));
      auto starty = std::exp(-kfac*(y*nmaxy));
      auto startz = std::exp(-kfac*(z*nmaxz));
    
      std::complex<double> facx, facy, facz;
      int ix, iy, iz;
      
      for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
	for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	  for (facz=startz, iz=0; iz<imz; iz++, facz*=stepz) {
	  
	    std::complex<double> fac = facx*facy*facz*expcoef(ix, iy, iz);
	  
	    // Compute wavenumber; recall that the coefficients are
	    // stored as follows: -nmax,-nmax+1,...,0,...,nmax-1,nmax
	    //
	    int ii = ix-nmaxx;
	    int jj = iy-nmaxy;
	    int kk = iz-nmaxz;
	  
	    // No contribution to acceleration and potential ("swindle")
	    // for zero wavenumber
	    if (ii==0 && jj==0 && kk==0) continue;
	    
	    // Limit to minimum wave number
	    if (abs(ii)<nminx || abs(jj)<nminy || abs(kk)<nminz) continue;
	  
	    // Normalization
	    double norm = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));

	    potl += fac*norm;
	    // dens += fac/norm;
	  
	    accx -= std::complex<double>(0.0, dfac*ii)*fac*norm;
	    accy -= std::complex<double>(0.0, dfac*jj)*fac*norm;
	    accz -= std::complex<double>(0.0, dfac*kk)*fac*norm;
	    
	  }
	}
      }
      
      cC->AddAcc(i, 0, accx.real());
      cC->AddAcc(i, 1, accy.real());
      cC->AddAcc(i, 2, accz.real());
      
      cC->AddPot(i, potl.real());

      // Deep debugging of acceleration
      if (deepDebug) {
	auto part = cC->Part(i);
	if (part->indx < 10) {
	  std::cout << "accel: [" << std::setw(8) << tnow << "] "
		    << std::setw(6) << part->indx;
	  for (int k=0; k<3; k++)
	    std::cout << std::setw(16) << part->pos[k];
	  for (int k=0; k<3; k++)
	    std::cout << std::setw(16) << part->vel[k];
	  for (int k=0; k<3; k++)
	    std::cout << std::setw(16) << part->acc[k];
	  std::cout
	    << std::setw(16) << expcoef(nmaxx+1, nmaxy, nmaxz).real()
	    << std::setw(16) << expcoef(nmaxx+1, nmaxy, nmaxz).imag()
	    << std::endl;
	}
      }
      // END: deep debugging
    }
    // END: particle loop
  }
  // END: multistep levels
    
  return (NULL);
}

void Cube::get_acceleration_and_potential(Component* C)
{
  cC = C;			// "Register" component
  nbodies = cC->Number();	// And compute number of bodies

  // Call the parallel acceleration and potential
  //
  MPL_start_timer();
  determine_acceleration_and_potential();
  MPL_stop_timer();

  // Clear external potential flag
  use_external = false;
}

void Cube::determine_acceleration_and_potential(void)
{
  exeTimer timer(this, "Force evaluation");

  if (play_back) {
    swap_coefs(&expcoefP, &expcoef);
  }

  if (use_external == false) {

    if (multistep && initializing) {
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
      timer.Start1();

      // Copy coefficients from this component to device
      //
      HtoD_coefs();
      //
      // Do the force computation
      //
      determine_acceleration_cuda();
      
      timer.Stop1();
    }
  } else {

    exp_thread_fork(false);

  }
#else

  exp_thread_fork(false);

#endif

  if (play_back) {
    swap_coefs(&expcoef, &expcoefP);
  }
}


void Cube::dump_coefs_h5(const std::string& file)
{
  // Add the current coefficients
  auto cur = std::make_shared<CoefClasses::CubeStruct>();

  cur->time     = tnow;
  cur->geom     = geoname[geometry];
  cur->id       = id;
  cur->time     = tnow;
  cur->nmaxx    = nmaxx;
  cur->nmaxy    = nmaxy;
  cur->nmaxz    = nmaxz;

  cur->allocate();		// Set the storage and copy the
				// coefficients through the map
  *cur->coefs   = expcoef;

  // Check if file exists
  //
  if (std::filesystem::exists(file)) {
    cubeCoefs.clear();
    cubeCoefs.add(cur);
    cubeCoefs.ExtendH5Coefs(file);
  }
  // Otherwise, extend the existing HDF5 file
  //
  else {
    // Copy the YAML config.  We only need this on the first call.
    std::ostringstream sout; sout << conf;
    cur->buf = sout.str();	// Copy to CoefStruct buffer

    // Add the name attribute.  We only need this on the first call.
    cubeCoefs.setName(component->name);

    // Add the default units
    cubeCoefs.setUnits({{"length", "none", 1.0},
			{"mass",   "none", 1.0},
			{"time",   "none", 1.0}});

    // And the new coefficients and write the new HDF5
    cubeCoefs.clear();
    cubeCoefs.add(cur);
    cubeCoefs.WriteH5Coefs(file);
  }
}


void Cube::multistep_update_begin()
{
  if (play_back and not play_cnew) return;

				// Clear the update matricies
  for (int n=0; n<nthrds; n++) {
    for (int M=mfirst[mdrft]; M<=multistep; M++) differ1[n][M].setZero();
  }
}

void Cube::multistep_update_finish()
{
  if (play_back and not play_cnew) return;

  // Combine the update matricies from all nodes
  //
  unsigned sz = (multistep - mfirst[mdrft] + 1)*osize;

  // Zero the buffer space
  //
  std::fill(pack.begin(), pack.begin() + sz, 0.0);

  // Pack the difference matrices
  //
  for (int M=mfirst[mdrft]; M<=multistep; M++) {

    unsigned offset = (M - mfirst[mdrft])*osize;

    for (int i=0; i<osize; i++)
      for (int n=0; n<nthrds; n++) {
	pack[offset + i] += differ1[n][M].data()[i];
      }
  }

  MPI_Allreduce (pack.data(), unpack.data(), sz, 
		 MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  
  // Update the local coefficients
  //
  for (int M=mfirst[mdrft]; M<=multistep; M++) {

    unsigned offset = (M - mfirst[mdrft])*osize;

    for (int i=0; i<osize; i++)
      expcoefN[M]->data()[i] += unpack[offset+i];
  }
}

void Cube::multistep_update(int from, int to, Component *c, int i, int id)
{
  if (play_back and not play_cnew) return;
  if (c->freeze(i)) return;

  double mass = c->Mass(i) * component->Adiabatic();

  double x = c->Pos(i, 0);
  double y = c->Pos(i, 1);
  double z = c->Pos(i, 2);
  
  // Truncate to cube with sides in [0,1]
  //
  if (wrap) {
    auto unitCube = [](double x) {
      if (x<0.0) x += std::floor(-x) + 1.0;
      else       x -= std::floor( x);
      return x;
    };
    
    x = unitCube(x);
    y = unitCube(y);
    z = unitCube(z);
  }

  // Only compute for points inside the unit cube
  //
  if (x<0.0 or x>1.0) return;
  if (y<0.0 or y>1.0) return;
  if (z<0.0 or z>1.0) return;
    
  // Recursion multipliers
  //
  std::complex<double> stepx = std::exp(-kfac*x);
  std::complex<double> stepy = std::exp(-kfac*y);
  std::complex<double> stepz = std::exp(-kfac*z);
    
  // Initial values for recursion
  //
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
	int ii = ix - nmaxx;
	int jj = iy - nmaxy;
	int kk = iz - nmaxz;
	
	if (ii==0 and jj==0 and kk==0) continue;

	// Normalization
	double norm = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));
	
	std::complex<double> val = -mass*facx*facy*facz*norm;

	differ1[id][from](ix, iy, iz) -= val;
	differ1[id][  to](ix, iy, iz) += val;
      }
    }
  }
}

void Cube::compute_multistep_coefficients()
{
  if (play_back and not play_cnew) return;

  // Clean coefficient matrix
  // 
  expcoef.setZero();
    
  // Interpolate to get coefficients above
  // 
  for (int M=0; M<mfirst[mdrft]; M++) {
    
    double numer = static_cast<double>(mdrft            - dstepL[M][mdrft]);
    double denom = static_cast<double>(dstepN[M][mdrft] - dstepL[M][mdrft]);

    double b = numer/denom;	// Interpolation weights
    double a = 1.0 - b;

    for (int i=0; i<osize; i++) {
      expcoef.data()[i] +=
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

    for (int i=0; i<osize; i++) {
      expcoef.data()[i] += expcoefN[M]->data()[i];
    }
  }
}

void Cube::writeCovarH5Params(HighFive::File& file)
{
  file.createAttribute<int>("nminx", HighFive::DataSpace::From(nminx)).write(nminx);
  file.createAttribute<int>("nminy", HighFive::DataSpace::From(nminy)).write(nminy);
  file.createAttribute<int>("nminz", HighFive::DataSpace::From(nminz)).write(nminz);
  file.createAttribute<int>("nmaxx", HighFive::DataSpace::From(nmaxx)).write(nmaxx);
  file.createAttribute<int>("nmaxy", HighFive::DataSpace::From(nmaxy)).write(nmaxy);
  file.createAttribute<int>("nmaxz", HighFive::DataSpace::From(nmaxz)).write(nmaxz);
}


PotAccel::CovarData Cube::getSubsample()
{
  CovarData elem;

  std::get<0>(elem) = sampleCounts;
  std::get<1>(elem) = sampleMasses;
  std::get<2>(elem) = Eigen::Tensor<std::complex<double>, 3>(sampT, 1, osize);
  std::get<3>(elem) = Eigen::Tensor<std::complex<double>, 4>(sampT, 1, osize, osize);

  // Fill the covariance structure with subsamples
  for (int T=0; T<sampT; T++) {
    for (int n1=0; n1<osize; n1++) {
      std::get<2>(elem)(T, 0, n1) = meanV[T](n1);
      for (int n2=0; n2<osize; n2++)
	if (fullCovar)
	  std::get<3>(elem)(T, 0, n1, n2) = covrV[T](n1, n2);
	else
	  std::get<3>(elem)(T, 0, n1, n2) = 0.0;
    } 
  }
    
  return elem;
}

