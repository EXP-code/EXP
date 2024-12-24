#include <filesystem>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <limits>

#include "expand.H"
#include <gaussQ.H>
#include <CylEXP.H>
#include <Cylinder.H>
#include <MixtureBasis.H>
#include <DiskDensityFunc.H>
#include <Timer.H>
#include <exputils.H>
#include <NVTX.H>

//@{
//! These are for testing exclusively (should be set false for production)
static bool cudaAccumOverride = false;
static bool cudaAccelOverride = false;
//@}


const std::set<std::string>
Cylinder::valid_keys = {
  "tk_type",
  "rcylmin",
  "rcylmax",
  "acyl",
  "hcyl",
  "hexp",
  "snr",
  "evcut",
  "nmaxfid",
  "lmaxfid",
  "mmax",
  "mlim",
  "ncylnx",
  "ncylny",
  "ncylr",
  "nmax",
  "ncylodd",
  "ncylrecomp",
  "npca",
  "npca0",
  "nvtk",
  "cachename",
  "eof_file",
  "override",
  "samplesz",
  "rnum",
  "pnum",
  "tnum",
  "ashift",
  "expcond",
  "precond",
  "logr",
  "pcavar",
  "pcaeof",
  "pcavtk",
  "pcadiag",
  "subsamp",
  "try_cache",
  "density",
  "EVEN_M",
  "cmap",
  "cmapr",
  "cmapz",
  "vflag",
  "self_consistent",
  "playback",
  "coefCompute",
  "coefMaster",
  "pyname"
};

Cylinder::Cylinder(Component* c0, const YAML::Node& conf, MixtureBasis *m) :
  Basis(c0, conf)
{
#if HAVE_LIBCUDA==1
  if (m) {
    throw GenericError("Error in Cylinder: MixtureBasis logic is not yet implemented in CUDA", __FILE__, __LINE__, 1016, false);
  }

  // Initialize the circular storage container 
  cuda_initialize();
  initialize_cuda_cyl = true;
#endif

  id              = "Cylinder";
  geometry        = cylinder;
  mix             = m;
  dof             = 3;
				// Default values

  rcylmin         = 0.001;	// Should only change these two in
  rcylmax         = 20.0;	// extreme circumstances

  ncylnx          = 256;	// These defaults should do fine in
  ncylny          = 128;	// most cases, as well
  ncylr           = 2000;

  acyl            = 0.01;
  lmaxfid         = 128;
  nmaxfid         = 64;
  mmax            = 6;
  mlim            = -1;
  hcyl            = 0.002;
  nmax            = 18;
  ncylodd         = 9;
  ncylrecomp      = -1;

  rnum            = 200;
  pnum            = 1;
  tnum            = 80;
  ashift          = 0.0;

  vflag           = 0;
  eof             = 1;
  npca            = 50;
  npca0           = 0;
  defSampT        = 1;
  hexp            = 1.0;
  snr             = 1.0;
  rem             = -1.0;
  self_consistent = true;
  firstime        = true;
  precond         = true;
  cmapR           = 1;
  cmapZ           = 1;
  logarithmic     = false;
  pcavar          = false;
  pcavtk          = false;
  pcadiag         = false;
  pcaeof          = false;
  subsamp         = false;
  nvtk            = 1;
  pcainit         = true;
  coef_dump       = true;
  try_cache       = true;
  dump_basis      = false;
  compute         = false;
  firstime_coef   = true;
  coefMaster      = true;
  lastPlayTime    = -std::numeric_limits<double>::max();
  EVEN_M          = false;
  cachename       = "";
#if HAVE_LIBCUDA==1
  cuda_aware      = true;
#endif

  initialize();

  // Enforce sane values for EOF integration
  //
  rnum = std::max<int>(10, rnum);
  pnum = std::max<int>(1,  pnum);
  tnum = std::max<int>(10, tnum);

  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = ncylnx;
  EmpCylSL::NUMY        = ncylny;
  EmpCylSL::NUMR        = ncylr;
  EmpCylSL::CMAPR       = cmapR;
  EmpCylSL::CMAPZ       = cmapZ;
  EmpCylSL::logarithmic = logarithmic;
  EmpCylSL::VFLAG       = vflag;

  if (cachename.size()==0)
    throw std::runtime_error("EmpCylSL: you must specify a cachename");

  // Make the empirical orthogonal basis instance
  //
  ortho = std::make_shared<CylEXP>
    (nmaxfid, lmaxfid, mmax, nmax, acyl, hcyl, ncylodd, cachename);
  
  // Set azimuthal harmonic order restriction?
  //
  if (mlim>=0)  ortho->set_mlim(mlim);
  if (EVEN_M)   ortho->setEven(EVEN_M);
  ortho->setSampT(defSampT);

  try {
    if (conf["tk_type"]) ortho->setTK(conf["tk_type"].as<std::string>());
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Cylinder: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("Cylinder: error in parsing tk_type");
  }

  
  // Cache file reading or generation
  //
  if (precond) {

    // Set parameters for external dcond function
    //
    eof      = 0;

    bool cache_ok = false;
    int  nOK      = 0;

    // Attempt to read EOF file from cache
    //
    if (std::filesystem::exists(cachename)) {

      cache_ok = ortho->read_cache();
      
      // If new cache is requested, backup existing basis file
      //
      if (!cache_ok) {
	  
	// Backup name for existing file
	string backupfile = cachename + ".bak";
	
	// Only root process renames
	if (myid==0) {
	  try {
	    std::filesystem::rename(cachename, backupfile);
	    std::cout << "---- Cylinder: parameter mismatch.  Renaming <"
		      << cachename << "> to <" << backupfile << "> and "
		      << "recreating the basis" << std::endl;
	  } catch (const std::filesystem::filesystem_error& e) {
	    std::cerr << "Cylinder: new cache requested but "
		      << "error creating backup file <" << backupfile
		      << "> from <" << cachename << ">, message: "
		      << e.code().message()
		      << std::endl;
	    // Set termination flag
	    nOK = 1;
	  }
	}
	
	// Broadcast termination request to all processes
	//
	MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	if (nOK)  {
	  throw std::runtime_error("Cylinder: error initializing from parameters");
	} 
      }
    }
    else {
      if (myid==0)
	std::cout << "---- Cylinder: can't find the EOF basis file <"
		  << cachename << ">" << std::endl;
    }

    // Genererate eof if needed
    //
    if (!cache_ok) {

      if (myid==0)
	std::cout << "---- Cylinder: generating the EOF basis file; "
		  << "this step will take many minutes."
		  << std::endl;
      
      std::shared_ptr<DiskDensityFunc> dfunc;
      if (pyname.size()) dfunc = std::make_shared<DiskDensityFunc>(pyname);

      auto DiskDens = [&](double R, double z, double phi)
      {
	if (dfunc) return (*dfunc)(R, z, phi);
	double f = cosh(z/hcyl);
	return exp(-R/acyl)/(4.0*M_PI*acyl*acyl*hcyl*f*f);
      };

      auto dcond = [&](double R, double z, double phi, int M)
      {
	//
	// No shift for M==0
	//
	if (M==0) return DiskDens(R, z, phi);
	
	//
	// Fold into [-PI/M, PI/M] for M>=1
	//
	double dmult = M_PI/M, phiS;
	if (phi>M_PI)
	  phiS = phi + dmult*(int)((2.0*M_PI - phi)/dmult);
	else
	  phiS = phi - dmult*(int)(phi/dmult);
  
	//
	// Apply a shift along the x-axis
	//
	double x = R*cos(phiS) - ashift*acyl;
	double y = R*sin(phiS);
	return DiskDens(sqrt(x*x + y*y), z, atan2(y, x));
      };

      ortho->generate_eof(rnum, pnum, tnum, dcond);
    }

    firstime = false;
  }

  // Make sure that all structures are initialized to start (e.g. for
  // multi-stepping but this should be done on 1st call to determine
  // coefs by default
  //
  ortho->setup_accumulation();

  
  std::string sep("----    ");

#ifdef DEBUG
  for (int n=0; n<numprocs; n++) {
    if (myid==n) {
      std::cout << std::endl << "Process " << myid << ": Cylinder parameters: "
		<< std::endl << sep << "nmaxfid="     << nmaxfid
		<< std::endl << sep << "lmaxfid="     << lmaxfid
		<< std::endl << sep << "mmax="        << mmax
		<< std::endl << sep << "mlim="        << mlim
		<< std::endl << sep << "nmax="        << nmax
		<< std::endl << sep << "ncylodd="     << ncylodd
		<< std::endl << sep << "rcylmin="     << rcylmin
		<< std::endl << sep << "rcylmax="     << rcylmax
		<< std::endl << sep << "acyl="        << acyl
		<< std::endl << sep << "hcyl="        << hcyl
		<< std::endl << sep << "precond="     << precond
		<< std::endl << sep << "pcavar="      << pcavar
		<< std::endl << sep << "pcaeof="      << pcaeof
		<< std::endl << sep << "nvtk="        << nvtk
		<< std::endl << sep << "npca="        << npca
		<< std::endl << sep << "npca0="       << npca0
		<< std::endl << sep << "pcadiag="     << pcadiag
		<< std::endl << sep << "cachename="   << cachename
		<< std::endl << sep << "selfgrav="    << std::boolalpha << self_consistent
		<< std::endl << sep << "logarithmic=" << logarithmic
		<< std::endl << sep << "vflag="       << vflag
		<< std::endl << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  if (myid==0) {
    std::cout << "---- Cylinder parameters: "
	      << std::endl << sep << "nmaxfid="     << nmaxfid
	      << std::endl << sep << "lmaxfid="     << lmaxfid
	      << std::endl << sep << "mmax="        << mmax
	      << std::endl << sep << "mlim="        << mlim
	      << std::endl << sep << "nmax="        << nmax
	      << std::endl << sep << "ncylodd="     << ncylodd
	      << std::endl << sep << "rcylmin="     << rcylmin
	      << std::endl << sep << "rcylmax="     << rcylmax
	      << std::endl << sep << "acyl="        << acyl
	      << std::endl << sep << "hcyl="        << hcyl
	      << std::endl << sep << "precond="     << precond
	      << std::endl << sep << "pcavar="      << pcavar
	      << std::endl << sep << "pcaeof="      << pcaeof
	      << std::endl << sep << "nvtk="        << nvtk
	      << std::endl << sep << "npca="        << npca
	      << std::endl << sep << "npca0="       << npca0
	      << std::endl << sep << "pcadiag="     << pcadiag
	      << std::endl << sep << "cachename="   << cachename
	      << std::endl << sep << "selfgrav="    << std::boolalpha << self_consistent
	      << std::endl << sep << "logarithmic=" << logarithmic
	      << std::endl << sep << "vflag="       << vflag
	      << std::endl;
  }
#endif
      
  // Test for basis consistency (will generate an exception if maximum
  // error is out of tolerance)
  //
  std::cout << "---- ";
  orthoTest(ortho->orthoCheck(), "Cylinder", "m");

  // Initialize internal variables
  //
  ncompcyl = 0;

  pos.resize(nthrds);
  frc.resize(nthrds);

#ifdef DEBUG
  offgrid.resize(nthrds);
#endif

}

Cylinder::~Cylinder()
{
  // Nothing
}

void Cylinder::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    // These first two should not be user settable . . . but need them for now
    //
    if (conf["rcylmin"   ])    rcylmin  = conf["rcylmin"   ].as<double>();
    if (conf["rcylmax"   ])    rcylmax  = conf["rcylmax"   ].as<double>();

    if (conf["acyl"      ])       acyl  = conf["acyl"      ].as<double>();
    if (conf["hcyl"      ])       hcyl  = conf["hcyl"      ].as<double>();
    if (conf["hexp"      ])       hexp  = conf["hexp"      ].as<double>();
    if (conf["snr"       ])        snr  = conf["snr"       ].as<double>();
    if (conf["evcut"     ])        rem  = conf["evcut"     ].as<double>();
    if (conf["nmaxfid"   ])    nmaxfid  = conf["nmaxfid"   ].as<int>();
    if (conf["lmaxfid"   ])    lmaxfid  = conf["lmaxfid"   ].as<int>();
    if (conf["mmax"      ])       mmax  = conf["mmax"      ].as<int>();
    if (conf["mlim"      ])       mlim  = conf["mlim"      ].as<int>();
    if (conf["ncylnx"    ])     ncylnx  = conf["ncylnx"    ].as<int>();
    if (conf["ncylny"    ])     ncylny  = conf["ncylny"    ].as<int>();
    if (conf["ncylr"     ])      ncylr  = conf["ncylr"     ].as<int>();
    if (conf["nmax"      ])       nmax  = conf["nmax"      ].as<int>();
    if (conf["ncylodd"   ])    ncylodd  = conf["ncylodd"   ].as<int>();
    if (conf["ncylrecomp"]) ncylrecomp  = conf["ncylrecomp"].as<int>();
    if (conf["npca"      ])       npca  = conf["npca"      ].as<int>();
    if (conf["npca0"     ])      npca0  = conf["npca0"     ].as<int>();
    if (conf["nvtk"      ])       nvtk  = conf["nvtk"      ].as<int>();
    if (conf["cachename" ])  cachename  = conf["cachename" ].as<std::string>();
    if (conf["eof_file"  ])  cachename  = conf["eof_file"  ].as<std::string>();
    if (conf["samplesz"  ])   defSampT  = conf["samplesz"  ].as<int>();
    
    if (conf["rnum"      ])       rnum  = conf["rnum"      ].as<int>();
    if (conf["pnum"      ])       pnum  = conf["pnum"      ].as<int>();
    if (conf["tnum"      ])       tnum  = conf["tnum"      ].as<int>();
    if (conf["ashift"    ])     ashift  = conf["ashift"    ].as<double>();
    if (conf["expcond"   ])    precond  = conf["expcond"   ].as<bool>();
    if (conf["precond"   ])    precond  = conf["precond"   ].as<bool>();
    if (conf["logr"      ]) logarithmic = conf["logr"      ].as<bool>();
    if (conf["pcavar"    ])     pcavar  = conf["pcavar"    ].as<bool>();
    if (conf["pcaeof"    ])     pcaeof  = conf["pcaeof"    ].as<bool>();
    if (conf["pcavtk"    ])     pcavtk  = conf["pcavtk"    ].as<bool>();
    if (conf["pcadiag"   ])    pcadiag  = conf["pcadiag"   ].as<bool>();
    if (conf["subsamp"   ])    subsamp  = conf["subsamp"   ].as<bool>();
    if (conf["try_cache" ])  try_cache  = conf["try_cache" ].as<bool>();
    if (conf["EVEN_M"    ])     EVEN_M  = conf["EVEN_M"    ].as<bool>();
    if (conf["cmap"      ])      cmapR  = conf["cmap"      ].as<int>();
    if (conf["cmapr"     ])      cmapR  = conf["cmapr"     ].as<int>();
    if (conf["cmapz"     ])      cmapZ  = conf["cmapz"     ].as<int>();
    if (conf["vflag"     ])      vflag  = conf["vflag"     ].as<int>();
    
    // Deprecation warning
    if (conf["expcond"]) {
      if (myid==0)
	std::cout << "---- Cylinder: parameter 'expcond' is deprecated. "
		  << "It has been renamed to 'precond'. The old parameter "
		  << "will be removed in a later release."
		  << std::endl;
    }

    // Deprecation warning
    if (conf["density"]) {
      if (myid==0)
	std::cout << "---- Cylinder: parameter 'density' is deprecated. "
		  << "The density field will be computed regardless."
		  << std::endl;
    }

    // Deprecation warning
    if (conf["override"]) {
      if (myid==0)
	std::cout << "---- Cylinder: parameter 'override' is deprecated. "
		  << "Basis will be recomputed if needed automatically."
		  << std::endl;
    }

    // Deprecation warning
    if (conf["eof_file"]) {
      if (myid==0)
	std::cout << "---- Cylinder: parameter 'eof_file' is deprecated. "
		  << "and will be removed in a future release. Please "
		  << "use 'cachename' instead."
		  << std::endl;
    }

    if (not conf["ncylodd"]) {
      ncylodd = nmax/4;
    }

    if (conf["self_consistent"])
      self_consistent = conf["self_consistent"].as<bool>();

    if (conf["playback"]) {
      std::string file = conf["playback"].as<std::string>();

      // Check that the file exists
      {
	std::ifstream test(file);
	if (not test) {
	  std::cerr << "Cylinder: process " << myid << " cannot open <"
		    << file << "> for reading" << std::endl;
	  throw std::runtime_error("Cylinder: error reading playback file");
	}
      }

      playback = std::dynamic_pointer_cast<CoefClasses::CylCoefs>(CoefClasses::Coefs::factory(file));

      if (not playback) {
	throw GenericError("Cylinder: failure in downcasting",
			   __FILE__, __LINE__, 1017, false);
      }
      
      // Set tolerance to 2 master time steps
      playback->setDeltaT(dtime*2);

      if (playback->nmax() != nmax) {
	std::ostringstream sout;
	sout << "Cylinder: nmax for playback [" << playback->nmax()
	     << "] does not match specification [" << nmax << "]";
	throw GenericError(sout.str(), __FILE__, __LINE__, 1018, false);
      }

      if (playback->mmax() != mmax) {
	if (myid==0) {
	  std::cerr << "Cylinder: mmax for playback [" << playback->mmax()
		    << "] does not match specification [" << mmax << "]"
		    << std::endl;
	}
	throw std::runtime_error("Cylinder: parameter mismatch");
      }

      P.resize(mmax+1, nmax);

      if (conf["coefCompute"]) play_cnew = conf["coefCompute"].as<bool>();

      if (play_cnew) P1.resize(mmax+1, nmax);

      play_back = true;

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

  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing Cylinder parameters: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("Cylinder::initialize: error parsing YAML");
  }
}

void Cylinder::get_acceleration_and_potential(Component* C)
{
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = std::make_shared<nvTracer>("Cylinder::get_acceleration");

  std::chrono::high_resolution_clock::time_point start0, start1, finish0, finish1;

  start0 = std::chrono::high_resolution_clock::now();

#ifdef DEBUG
  cout << "Process " << myid 
       << ": in Cylinder::get_acceleration_and_potential" << endl;
#endif
				
  cC = C;

  //====================================================
  // Accel & pot using previously computed coefficients 
  //====================================================

  if (use_external) {
    nvTracerPtr tPtr1;
    if (cuda_prof) {
      tPtr1 = std::make_shared<nvTracer>("Cylinder: in external");
    }

    MPL_start_timer();
    determine_acceleration_and_potential();
    MPL_stop_timer();

    use_external = false;

    return;
  }


  //======================================
  // Determine potential and acceleration 
  //======================================

  MPL_start_timer();

  determine_acceleration_and_potential();

  MPL_stop_timer();

  //=======================
  // Recompute PCA analysis
  //=======================

				// No recomputation ever if the
  if (!precond) {		// basis has been precondtioned

				// Only do this check only once per
				// multistep; might as well be at 
				// the end of the multistep sequence
    if ((multistep==0 || mstep==0) && !initializing) {
      ncompcyl++;
      if (ncompcyl == ncylrecomp) {
	ncompcyl = 0;
	eof = 1;
	determine_coefficients();
      }
    }

  }


  //=================
  // Debugging output
  //=================
  if (VERBOSE>3 && myid==1 && component->EJ) {
    std::string toutfile = outdir + "test.orientation." + runtag;
    std::ofstream debugf(toutfile.c_str(), ios::app);
    auto axis = component->orient->currentAxis();
    debugf << tnow << " "
	   << component->orient->currentAxis()[0] << " " 
	   << component->orient->currentAxis()[1] << " " 
	   << component->orient->currentAxis()[2] << " " 
	   << component->orient->currentAxisVar() << " "
	   << component->orient->currentCenter()[0] << " " 
	   << component->orient->currentCenter()[1] << " " 
	   << component->orient->currentCenter()[2] << " " 
	   << component->orient->currentCenterVar() << " "
	   << component->orient->currentCenterVarZ() << " "
	   << component->orient->currentE() << " "
	   << component->orient->currentUsed()
	   << endl;
  }

}

std::ostream& operator<< (std::ostream& os, const Eigen::Vector3d& p)
{
  std::streamsize sp = os.precision();
  os.precision(6);
  for (int i=0; i<3; i++) os << std::setw(16) << p[i];
  os.precision(sp);
  return os;
}


void * Cylinder::determine_coefficients_thread(void * arg)
{
  double r, r2, phi, R2;
  double xx, yy, zz, mas;
  double Rmax2 = rcylmax*rcylmax*acyl*acyl;

  int id = *((int*)arg);
  int indx, nbeg, nend, nbodies;

  use[id] = 0;
  cylmass0[id] = 0.0;

  thread_timing_beg(id);

  vector<double> ctr;
  if (mix) mix->getCenter(ctr);

  if (eof) {

    // Will use all of the bodies independent of level
    //
    nbodies = cC->Number();
    
    if (nbodies==0) {
      thread_timing_end(id);
      return (NULL);
    }
    nbeg = nbodies*id/nthrds;
    nend = nbodies*(id+1)/nthrds;

    unsigned indx;
    PartMapItr n = cC->Particles().begin();
    for (int i=0; i<nbeg; i++) n++; // Move to beginning iterator

    for (int i=nbeg; i<nend; i++, n++) {

      indx = n->first;

      // Frozen particles don't contribute to field
      //
      if (cC->freeze(indx)) continue;
    
      if (mix) {
	for (int j=0; j<3; j++) 
	  pos[id][j] = cC->Pos(indx, j, Component::Local) - ctr[j];
      } else {
	for (int j=0; j<3; j++) 
	  pos[id][j] = cC->Pos(indx, j, 
			       Component::Local | Component::Centered);
      }
      
      if ( (cC->EJ & Orient::AXIS) && !cC->EJdryrun) 
	pos[id] = cC->orient->transformBody() * pos[id];

      xx = pos[id][0];
      yy = pos[id][1];
      zz = pos[id][2];

      r2 = xx*xx + yy*yy;
      r = sqrt(r2);
      R2 = r2 + zz*zz;
    
      if ( R2 < Rmax2) {

	mas = cC->Mass(indx);
	phi = atan2(yy, xx);

	ortho->accumulate_eof(r, zz, phi, mas, id, cC->Part(indx)->level);
	
	use[id]++;
	cylmass0[id] += mas;

      } 
    }

  } else {

    nbodies = cC->levlist[mlevel].size();
    
    if (nbodies==0) {
      thread_timing_end(id);
      return (NULL);
    }
    nbeg = nbodies*id/nthrds;
    nend = nbodies*(id+1)/nthrds;

    double adb = component->Adiabatic();

    for (int i=nbeg; i<nend; i++) {

      indx = cC->levlist[mlevel][i];

      // Frozen particles don't contribute to field
      //
      if (cC->freeze(indx)) continue;
    
      for (int j=0; j<3; j++) 
	pos[id][j] = cC->Pos(indx, j, Component::Local | Component::Centered);

      if ( (cC->EJ & Orient::AXIS) && !cC->EJdryrun) 
	pos[id] = cC->orient->transformBody() * pos[id];

      xx = pos[id][0];
      yy = pos[id][1];
      zz = pos[id][2];

      r2 = xx*xx + yy*yy;
      r = sqrt(r2);
      R2 = r2 + zz*zz;
    
      if ( R2 < Rmax2) {

	mas = cC->Mass(indx) * adb;
	phi = atan2(yy, xx);

	ortho->accumulate(r, zz, phi, mas, indx, id, mlevel, compute);

	use[id]++;
	cylmass0[id] += mas;
	
      } else {

	if (VERBOSE>6) {
	  cout << "Process " << myid 
	       << ": r^2=" << R2
	       << " max r^2=" << Rmax2 
	       << " r2=" << r2 
	       << " z2=" << zz*zz 
	       << " m=" << cylmass0[id] 
	       << " eof=" << eof
	       << endl;

	  if (std::isnan(R2)) {
	    cout << endl
		 << cC->orient->transformBody() << endl
		 << cC->orient->currentAxis()   << endl;
	    throw GenericError("Squared radius is NaN",
			       __FILE__, __LINE__, -1, true);
	  }
	}
      }
      
    }
  }

  thread_timing_end(id);

  return (NULL);
}

void Cylinder::determine_coefficients(void)
{
  if (play_back) {
    determine_coefficients_playback();
    if (play_cnew) determine_coefficients_particles();
  } else {
    determine_coefficients_particles();
  }
}

void Cylinder::determine_coefficients_playback(void)
{
  compute_grid_mass();	// Only performed once to start

  // Do we need new coefficients?
  if (tnow <= lastPlayTime) return;
  lastPlayTime = tnow;

  if (coefMaster) {
    
    if (myid==0) {
      auto ret = playback->interpolate(tnow);

      P = std::get<0>(ret);
      if (not std::get<1>(ret)) stop_signal = 1;

      int nrows = P.rows(), ncols = P.cols();
      MPI_Bcast(&nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(P.data(), P.size(), MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    } else {
      int nrows, ncols;
      MPI_Bcast(&nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
      P.resize(nrows, ncols);
      MPI_Bcast(P.data(), P.size(), MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    }
    
  } else {
    auto ret = playback->interpolate(tnow);

    P = std::get<0>(ret);
    if (not std::get<1>(ret)) stop_signal = 1;
  }

}

void Cylinder::determine_coefficients_particles(void)
{
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = std::make_shared<nvTracer>("Cylinder::determine_coefficients");

  std::chrono::high_resolution_clock::time_point start0, start1, finish0, finish1;

  start0 = std::chrono::high_resolution_clock::now();


  static char routine[] = "determine_coefficients_Cylinder";

  if (!self_consistent && !firstime_coef && !initializing) return;

  if (!precond && firstime) {
				// Try to read cache
    bool cache_ok = false;
    if (try_cache || restart) {
      cache_ok = ortho->read_cache();
				// For a restart, cache must be read
				// otherwise, abort
      if (restart && !cache_ok) {
	if (myid==0) 
	  std::cerr << "Cylinder: can not read cache file on restart" << endl;
	throw std::runtime_error("Cylinder: cache read error");
      }
    }

    eof = cache_ok ? 0 : 1;
    firstime = false;
				// If we can't read the cache, or the cache
				// does not match the requested parameters,
				// remake the emperical orthogonal basis 
    if (eof) {
      determine_coefficients_eof();
    }
  }

  if ( (pcavar or pcaeof) and pcainit) {
    EmpCylSL::PCAVAR = pcavar;
    EmpCylSL::PCAEOF = pcaeof;
    EmpCylSL::PCAVTK = pcavtk;
    EmpCylSL::VTKFRQ = nvtk;
    EmpCylSL::HEXP   = hexp;
    std::ostringstream sout;
    if (pcadiag) 
      sout << runtag << ".pcadiag." << cC->id << "." << cC->name;
    ortho->setHall(sout.str(), component->NewTotal());
    if (myid==0) {
      std::cout << "Cylinder: PCA initialized";
      if (pcadiag) 
	std::cout << ", writing diagnostic output to <"
		  << sout.str() << ">";
      std::cout << std::endl;
    }
    pcainit = false;
  }

  if (pcavar or pcaeof) {
    if (this_step >= npca0)
      compute = (mstep == 0) && !( (this_step-npca0) % npca);
    else
      compute = false;
  }

  ortho->setup_accumulation(mlevel);

  cylmass0.resize(nthrds);

#ifdef LEVCHECK
  for (int n=0; n<numprocs; n++) {
    if (n==myid) {
      if (myid==0) cout << "------------------------" << endl
			<< "Level check in Cylinder:" << endl 
			<< "------------------------" << endl;
      cout << setw(4) << myid << setw(4) << mlevel << setw(4) << eof;
      if (cC->levlist[mlevel].size())
	cout << setw(12) << cC->levlist[mlevel].size()
	     << setw(12) << cC->levlist[mlevel].front()
	     << setw(12) << cC->levlist[mlevel].back() << endl;
      else
	cout << setw(12) << cC->levlist[mlevel].size()
	     << setw(12) << (int)(-1)
	     << setw(12) << (int)(-1) << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0) cout << endl;
#endif
    
#if HAVE_LIBCUDA==1
  if (component->cudaDevice>=0 and use_cuda) {
    if (cudaAccumOverride) {
      component->CudaToParticles();
      exp_thread_fork(true);
    } else {
      start1 = std::chrono::high_resolution_clock::now();
      if (mstep==0) {
	std::fill(use.begin(), use.end(), 0.0);
	std::fill(cylmass0.begin(), cylmass0.end(), 0.0);
      }
      determine_coefficients_cuda(compute);
      DtoH_coefs(mlevel);
      finish1 = std::chrono::high_resolution_clock::now();
    }
  } else {    
    exp_thread_fork(true);
  }
#else
				// Threaded coefficient accumulation loop
  exp_thread_fork(true);
#endif
				// Accumulate counts and mass used to
				// determine coefficients
  int use1=0, use0=0;
  double cylmassT1=0.0, cylmassT0=0.0;
  
  for (int i=0; i<nthrds; i++) {
    use1      += use[i];
    cylmassT1 += cylmass0[i];
  }
				// Turn off timer so as not bias by 
				// communication barrier
  MPL_stop_timer();

  if (not play_back and tnow==resetT) {

    MPI_Allreduce ( &use1, &use0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce ( &cylmassT1, &cylmassT0, 1, MPI_DOUBLE, MPI_SUM, 
		    MPI_COMM_WORLD );

    used    += use0;
    cylmass += cylmassT0;
  }

  MPL_start_timer();

  //=========================
  // Make the coefficients
  // for this level
  //=========================

  if (multistep==0 || !self_consistent) {
    ortho->make_coefficients(compute);
  } else {
    ortho->make_coefficients(mfirst[mstep], compute);
    compute_multistep_coefficients(); // I don't think this is necessary . . .
  }

  //=========================
  // Compute Hall smoothing
  //=========================

  if ((pcavar or pcaeof) and mlevel==0) ortho->pca_hall(compute, subsamp);

  //=========================
  // Apply Hall smoothing
  //=========================

  if (pcavar and used) ortho->set_trimmed(snr, rem);

  //=========================
  // Dump basis on first call
  //=========================

  if ( dump_basis and (this_step==0 || (precond and ncompcyl==0) )
       && ortho->coefs_made_all() && !initializing) {

    if (myid == 0 and multistep==0 || mstep==0) {
      
      nvTracerPtr tPtr2;
      if (cuda_prof) {
	tPtr2 = std::make_shared<nvTracer>("Cylinder::dump basis");
      }

      ortho->dump_basis(runtag.c_str(), this_step);
      
      ostringstream dumpname;
      dumpname << "images" << "." << runtag << "." << this_step;
      ortho->dump_images(dumpname.str(), 5.0*acyl, 5.0*hcyl, 64, 64, true);
      //
      // This next call is ONLY for deep debug
      //
      // dump_mzero(runtag.c_str(), this_step);
    }
  }

  print_timings("Cylinder: coefficient timings");

  finish0 = std::chrono::high_resolution_clock::now();
  
#if HAVE_LIBCUDA==1
  if (component->timers) {
    std::chrono::duration<double> duration0 = finish0 - start0;
    std::chrono::duration<double> duration1 = finish1 - start1;
    
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "== Coefficient evaluation [Cylinder] level="
	      << mlevel << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Time in CPU: " << duration0.count()-duration1.count() << std::endl;
    if (cC->cudaDevice>=0) {
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
	      << "-- Cylinder T=" << std::setw(16) << tnow << std::endl
	      << std::string(60, '-') << std::endl;
    ortho->dump_coefs(std::cout);
    std::cout << std::string(60, '-') << std::endl;
  }

  firstime_coef = false;
}


void Cylinder::determine_coefficients_eof(void)
{
  if (eof==0) return;

  static char routine[] = "determine_coefficients_eof_Cylinder";
  
  ortho->setup_eof();
  ortho->setup_accumulation();

  cylmass = 0.0;
  if (myid==0) cerr << "Cylinder: setup for eof\n";

  cylmass0.resize(nthrds);

				// Threaded coefficient accumulation loop
  exp_thread_fork(true);

				// Accumulate counts and mass used to
				// determine coefficients
  int use0=0, use1=0;
  double cylmassT1=0.0, cylmassT0=0.0;

  for (int i=0; i<nthrds; i++) {
    use1 += use[i];
    cylmassT1 += cylmass0[i];
  }

				// Turn off timer so as not bias by 
				// communication barrier
  MPL_stop_timer();

  MPI_Allreduce ( &use1, &use0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce ( &cylmassT1, &cylmassT0, 1, MPI_DOUBLE, MPI_SUM, 
		  MPI_COMM_WORLD );

  MPL_start_timer();

  if (myid==0) cerr << "Cylinder: eof grid mass=" << cylmassT0 
		    << ", number=" << use0 << "\n";

  ortho->make_eof();
  if (myid==0) cerr << "Cylinder: eof computed\n";

  ortho->make_coefficients();
  if (myid==0) cerr << "Cylinder: coefs computed\n";

  eof = 0;
}


void check_force_values(double phi, double p, double fr, double fz, double fp)
{
  if (
      std::isinf(phi) || std::isnan(phi) ||
      std::isinf(p  ) || std::isnan(p  ) ||
      std::isinf(fr ) || std::isnan(fr ) ||
      std::isinf(fz ) || std::isnan(fz ) ||
      std::isinf(fp ) || std::isnan(fp ) ) 
    {
      cerr << "check_force_values: Illegal value\n";
    }
}


void * Cylinder::determine_acceleration_and_potential_thread(void * arg)
{
  double r, r2, r3, phi;
  double xx, yy, zz;
  double p, p0, fr, fz, fp, pa;

  constexpr double ratmin = 0.75;
  constexpr double maxerf = 3.0;
  constexpr double midpt  = ratmin + 0.5*(1.0 - ratmin);
  constexpr double rsmth  = 0.5*(1.0 - ratmin)/maxerf;

  double ratio, frac, cfrac, mfactor = 1.0;

  // Get the grid actual radius from EmpCylSL
  //
  double R2 = ortho->get_ascale()*ortho->get_rtable();
  R2 = R2*R2;			// Compute the square

  vector<double> ctr;
  if (mix) mix->getCenter(ctr);

  int id = *((int*)arg);

#ifdef DEBUG
  static bool firstime = true;
  std::ofstream out;
  int flg;
  if (firstime && myid==0 && id==0) out.open("debug.tst");
#endif

  thread_timing_beg(id);

  // If we are multistepping, compute accel only at or below <mlevel>
  //
  for (unsigned lev=mlevel; lev<=multistep; lev++) {

    unsigned nbodies = cC->levlist[lev].size();

    if (nbodies==0) continue;

    int nbeg = nbodies*id/nthrds;
    int nend = nbodies*(id+1)/nthrds;
    
#ifdef DEBUG
    cout << "Process " << myid << " id=" << id 
	 << ": nbodies=" << nbodies
	 << " lev=" << lev
	 << " nbeg=" << nbeg
	 << " nend=" << nend << endl;
#endif

    for (int q=nbeg; q<nend; q++) {

      unsigned indx = cC->levlist[lev][q];

      // Deep debug
      /*
      if (indx<1 || indx>cC->CurTotal()) {
	cout << "Process " << myid << " id=" << id 
	     << ": index error in Cylinder q=" << q
	     << " indx=" << indx << endl;
      }
      */

      if (mix) {

	if (use_external) {
	  cC->Pos(pos[id].data(), indx, Component::Inertial);
	  component->ConvertPos(pos[id].data(), Component::Local);
	} else
	  cC->Pos(pos[id].data(), indx, Component::Local);

	// Only apply this fraction of the force
	mfactor = mix->Mixture(pos[id].data());
	for (int k=0; k<3; k++) pos[id][k] -= ctr[k];

      } else {

	if (use_external) {
	  cC->Pos(pos[id].data(), indx, Component::Inertial);
	  component->ConvertPos(pos[id].data(), Component::Local | Component::Centered);
	} else
	  cC->Pos(pos[id].data(), indx, Component::Local | Component::Centered);

      }

      if ( (component->EJ & Orient::AXIS) && !component->EJdryrun) 
	pos[id] = component->orient->transformBody() * pos[id];

      xx    = pos[id][0];
      yy    = pos[id][1];
      zz    = pos[id][2];
      
      r2    = xx*xx + yy*yy;
      r     = sqrt(r2) + DSMALL;
      phi   = atan2(yy, xx);
      pa    = 0.0;

      ratio = sqrt( (r2 + zz*zz)/R2 );

      if (ratio >= 1.0) {
	frac       = 0.0;
	cfrac      = 1.0;
	frc[id][0] = 0.0;
	frc[id][1] = 0.0;
	frc[id][2] = 0.0;
      } else if (ratio > ratmin) {
	frac  = 0.5*(1.0 - erf( (ratio - midpt)/rsmth ));
	cfrac = 1.0 - frac;
      } else {
	cfrac = 0.0;
	frac  = 1.0;
      }
	
      cfrac *= mfactor;
      frac  *= mfactor;

      if (ratio < 1.0) {

	ortho->accumulated_eval(r, zz, phi, p0, p, fr, fz, fp);
#ifdef DEBUG
	check_force_values(phi, p, fr, fz, fp);
#endif
	frc[id][0] = ( fr*xx/r - fp*yy/r2 ) * frac;
	frc[id][1] = ( fr*yy/r + fp*xx/r2 ) * frac;
	frc[id][2] = fz * frac;
	pa         = p  * frac;
	
#ifdef DEBUG
	flg = 1;
#endif
      }

      if (ratio > ratmin) {

	r3 = r2 + zz*zz;
	p = -cylmass/sqrt(r3);	// -M/r
	fr = p/r3;		// -M/r^3

	frc[id][0] += xx*fr * cfrac;
	frc[id][1] += yy*fr * cfrac;
	frc[id][2] += zz*fr * cfrac;
	pa         += p     * cfrac;

#ifdef DEBUG
	offgrid[id]++;
	flg = 2;
#endif
      }
    
      cC->AddPot(indx, pa);

      if ( (component->EJ & Orient::AXIS) && !component->EJdryrun) 
	frc[id] = component->orient->transformOrig() * frc[id];

      for (int j=0; j<3; j++) cC->AddAcc(indx, j, frc[id][j]);

#ifdef DEBUG
      if (firstime && myid==0 && id==0 && q < 5) {
	out << setw(9)  << q          << endl
	    << setw(9)  << indx       << endl
	    << setw(9)  << flg        << endl
	    << setw(18) << xx         << endl
	    << setw(18) << yy         << endl
	    << setw(18) << zz         << endl
	    << setw(18) << frc[0][0]  << endl
	    << setw(18) << frc[0][1]  << endl
	    << setw(18) << frc[0][2]  << endl;
      }
#endif
    }
  }

#ifdef DEBUG
  firstime = false;		// DEBUG
#endif

  thread_timing_end(id);

  return (NULL);
}

static int ocf = 0;

void Cylinder::determine_acceleration_and_potential(void)
{
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = std::make_shared<nvTracer>("Cylinder::determine_acceleration");

  std::chrono::high_resolution_clock::time_point start0, start1, finish0, finish1;

  start0 = std::chrono::high_resolution_clock::now();

  static char routine[] = "determine_acceleration_and_potential_Cyl";
  
  if (play_back) {
    if (play_cnew) getCoefs(P1);
    setCoefs(P);
  }

  if (use_external == false) {

    if (multistep && (self_consistent || initializing)) {
      compute_multistep_coefficients();
    }

  }

#ifdef DEBUG
  for (int i=0; i<nthrds; i++) offgrid[i] = 0;
  cout << "Process " << myid << ": about to fork" << endl;
#endif

#if HAVE_LIBCUDA==1
  if (use_cuda and cC->cudaDevice>=0 and cC->force->cudaAware()) {
    if (cudaAccelOverride) {
      cC->CudaToParticles();
      exp_thread_fork(false);
      cC->ParticlesToCuda();
    } else {
      start1 = std::chrono::high_resolution_clock::now();
      //
      // Copy coeficients from this component to device
      //
      HtoD_coefs();
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
  cout << "Cylinder: process " << myid << " returned from fork" << endl;
  int offtot=0;
  for (int i=1; i<nthrds; i++) offgrid[0] += offgrid[i];
  MPI_Reduce(&offgrid[0], &offtot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid==0) {
    if (use_external)
      cout << endl << "T=" << tnow << "  external offgrid=" << offtot << endl;
    else
      cout << endl << "T=" << tnow << "  self offgrid=" << offtot << endl;
  }    

  unsigned long imin = std::numeric_limits<unsigned long>::max();
  unsigned long imax = 0, kmin = imin, kmax = 0;

  for (auto p : cC->Particles()) {
    imin = std::min<unsigned long>(imin, p.first);
    imax = std::max<unsigned long>(imax, p.first);
    kmin = std::min<unsigned long>(kmin, p.second->indx);
    kmax = std::max<unsigned long>(kmax, p.second->indx);
  }

  cout << "Cylinder: process " << myid << " name=<" << cC->name << "> bodies ["
       << kmin << ", " << kmax << "], ["
       << kmin << ", " << kmax << "]"
       << " #=" << cC->Particles().size() << endl;
#endif

  if (play_back) {
    getCoefs(P);
    if (play_cnew) setCoefs(P1);
  }

  print_timings("Cylinder: acceleration timings");


# if HAVE_LIBCUDA
  if (component->timers) {
    auto finish0 = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> duration0 = finish0 - start0;
    std::chrono::duration<double> duration1 = finish1 - start1;
    std::chrono::duration<double> duration2 = start1  - start0;

    std::cout << std::string(60, '=') << std::endl;
    std::cout << "== Force evaluation [Cylinder::" << cC->name
	      << "] level=" << mlevel << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Time in CPU: " << duration0.count()-duration1.count() << std::endl;
    if (cC->cudaDevice>=0) {
      std::cout << "Time in GPU: " << duration1.count() << std::endl;
      std::cout << "Time before: " << duration2.count() << std::endl;
    }
    std::cout << std::string(60, '=') << std::endl;
  }
#endif
}

void Cylinder::determine_fields_at_point
(double x, double y, double z, 
 double *tdens0, double *tpotl0, 
 double *tdens, double *tpotl, 
 double *tpotX, double *tpotY, double *tpotZ)
{
  double R   = sqrt(x*x + y*y);
  double phi = atan2(y, x);
  double cph = cos(phi), sph = sin(phi);

  double tpotR, tpotP;

  determine_fields_at_point_cyl(R, z, phi, tdens0, tpotl0, tdens, tpotl, 
				&tpotR, tpotZ, &tpotP);
  
  *tpotX = tpotR*cph - tpotP*sph ;
  *tpotY = tpotR*sph + tpotP*cph ;
}


void Cylinder::
determine_fields_at_point_sph(double r, double theta, double phi,
			      double *tdens0, double *tpotl0, 
			      double *tdens, double *tpotl, 
			      double *tpotr, double *tpott, 
			      double *tpotp)

{
  double R = r*sin(theta);
  double z = r*cos(theta);
  double tpotR, tpotZ;

  determine_fields_at_point_cyl(R, z, phi, tdens0, tpotl0, tdens, tpotl, 
				&tpotR, &tpotZ, tpotp);
  
  *tpotr =   tpotR*sin(theta) + tpotZ*cos(theta) ;
  *tpott = (-tpotZ*sin(theta) + tpotR*cos(theta) )/(r+1.0e-10);
}



void Cylinder::determine_fields_at_point_cyl(double r, double z, double phi,
					     double *tdens0, double *tpotl0, 
					     double *tdens, double *tpotl, 
					     double *tpotr, double *tpotz, double *tpotp)
{
  ortho->accumulated_eval(r, z, phi, *tpotl0, *tpotl, *tpotr, *tpotz, *tpotp);
  // Accumulated eval returns forces not potential gradients
  *tpotr *= -1.0;
  *tpotz *= -1.0;
  *tpotp *= -1.0;
  *tdens = ortho->accumulated_dens_eval(r, z, phi, *tdens0);
}

				// Dump coefficients to a file
void Cylinder::dump_coefs(ostream& out)
{
  ortho->dump_coefs_binary(out, tnow);
}

// Dump coefficients to an HDF5 file

void Cylinder::dump_coefs_h5(const std::string& file)
{
  // Add the current coefficients
  auto cur = std::make_shared<CoefClasses::CylStruct>();

  cur->time = tnow;
  cur->geom = geoname[geometry];
  cur->id   = id;
  cur->mmax = mmax;
  cur->nmax = nmax;

  cur->allocate();

  auto & cof = *cur->coefs;	// Reference to the coefficient map

  Eigen::VectorXd cos1(nmax), sin1(nmax);
  
  for (int m=0; m<=mmax; m++) {

    ortho->get_coefs(m, cos1, sin1);

    for (int ir=0; ir<nmax; ir++) {
      if (m==0) {
	cof(m, ir) = {cos1(ir), 0.0};
      } else {
	cof(m, ir) = {cos1(ir), sin1(ir)};
      }
    }
  }

  // Add center
  //
  cur->ctr = component->getCenter(Component::Local | Component::Centered);

  // Check if file exists
  //
  if (std::filesystem::exists(file)) {
    cylCoefs.clear();
    cylCoefs.add(cur);
    cylCoefs.ExtendH5Coefs(file);
  } else {
    // Copy the YAML config.  We only need this on the first call.
    std::ostringstream sout; sout << conf;
    cur->buf = sout.str();	// Copy to CoefStruct buffer

    // Add the name attribute.  We only need this on the first call.
    cylCoefs.setName(component->name);

    // Add the new coefficients and write the new HDF5
    cylCoefs.clear();
    cylCoefs.add(cur);
    cylCoefs.WriteH5Coefs(file);
  }
}

				// Density debug
#include <fstream>

void Cylinder::dump_mzero(const string& name, int step)
{
  const double RMAX = 5.0*acyl;
  const double ZMAX = 5.0*hcyl;
  double r, dr = RMAX/(ncylnx-1);
  double z, dz = 2.0*ZMAX/(ncylny-1);

  float zz;
  string label[] = {".dens0.", ".pot0.", ".fr0.", ".fz0."};
  std::ofstream out[4];

  for (int i=0; i<4; i++) {
    ostringstream ins;
    ins << name << label[i] << step;
    out[i].open(ins.str());

    out[i].write((char *)&ncylnx, sizeof(int));
    out[i].write((char *)&ncylny, sizeof(int));
    out[i].write((char *)&(zz=  0.0), sizeof(float));
    out[i].write((char *)&(zz= RMAX), sizeof(float));
    out[i].write((char *)&(zz=-ZMAX), sizeof(float));
    out[i].write((char *)&(zz= ZMAX), sizeof(float));
  }


				// Ok, write data
  double p, p0, d0, fr, fz, fp;

  for (int k=0; k<ncylny; k++) {

    z = -ZMAX + dz*k;
	
    for (int j=0; j<ncylnx; j++) {
	  
      r = dr*j;

      zz = ortho->accumulated_dens_eval(r, z, 0.0, d0);
      out[0].write((char *)&zz, sizeof(float));

      ortho->accumulated_eval(r, z, 0.0, p0, p, fr, fz, fp);
      out[1].write((char *)&(zz=p ), sizeof(float));
      out[2].write((char *)&(zz=fr), sizeof(float));
      out[3].write((char *)&(zz=fz), sizeof(float));
    }
  }

				// Close and delete streams
  for (int i=0; i<4; i++) {
    try {
      out[i].close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "Cylinder: exception closing file <"
		<< name << label[i] << step << ">: "
		<< e.what() << std::endl;
    }

  }

}

void Cylinder::multistep_update(int from, int to, Component* c, int i, int id)
{
  if (play_back)        return;
  if (!self_consistent) return;
  if (c->freeze(i))     return;

  double mass = c->Mass(i) * component->Adiabatic();

  double xx = c->Pos(i, 0, Component::Local | Component::Centered);
  double yy = c->Pos(i, 1, Component::Local | Component::Centered);
  double zz = c->Pos(i, 2, Component::Local | Component::Centered);

  double r2 = (xx*xx + yy*yy);
  double  r = sqrt(r2);
  double phi = atan2(yy, xx);

  ortho->multistep_update(from, to, r, zz, phi, mass, id);

#ifdef CYL_UPDATE_TABLE
  occt[from][to]++;
#endif
}

void Cylinder::compute_multistep_coefficients() 
{
  if (play_back and not play_cnew) return;
  
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = nvTracerPtr(new nvTracer("Cylinder::compute_multistep_coefficients"));
  ortho->compute_multistep_coefficients(mfirst[mstep]);
}


void Cylinder::multistep_reset() 
{ 
  if (play_back) return;
  
  used    = 0; 
  cylmass = 0.0;
  resetT  = tnow;
  ortho->reset_mass();
  ortho->multistep_reset();
}


static int idbg = 0;
void Cylinder::multistep_debug() 
{
  if (myid==0) {
    cout << endl;
    cout << setw(70) << setfill('-') << '-' << endl;
    ostringstream sout;
    sout << "--- multistep_debug: " << idbg << endl;
    cout << setw(70) << left << sout.str() << endl << right;
    cout << setw(70) << '-' << setfill(' ') << endl;

    ostringstream sout2;
    sout2 << "cylinder.coefs." << runtag << "." << ocf++;
    ofstream out(sout2.str().c_str());
    ortho->dump_coefs(out);
  }

  ortho->multistep_debug();

  if (myid==1) ortho->dump_basis(runtag.c_str(), idbg);

  ostringstream dumpname;
  dumpname << "images" << "." << runtag << "." << idbg;
  ortho->dump_images(dumpname.str(), 5.0*acyl, 5.0*hcyl, 64, 64, true);
  dump_mzero(runtag.c_str(), idbg);
  
  idbg++;
}


void Cylinder::compute_grid_mass()
{
  static bool done = false;
  
  if (done) return;
  done = true;

  // Compute used and cylmass for playback (this means that cylmass
  // will not be the same as the original simulation but it should be
  // close unless the original grid was inappropriate.
  //
  cylmass = 0.0;
  used    = 0;
    
  double Rmax2 = rcylmax*rcylmax*acyl*acyl;
  
  auto p = cC->Particles();
  
  for (auto it=p.begin(); it!=p.end(); ++it) {
    auto n = it->first;

    double R2 = 0.0;
    for (int j=0; j<3; j++)  {
      double pos = cC->Pos(n, j, Component::Local | Component::Centered);
      R2 += pos*pos;
    }
      
    if ( R2 < Rmax2) {
      cylmass += cC->Mass(n);
      used    += 1;
    } 
  }
  // END: particle loop

  MPI_Allreduce(MPI_IN_PLACE, &cylmass, 1, MPI_DOUBLE, MPI_SUM,
		MPI_COMM_WORLD);

  MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);

  ortho->set_mass(cylmass);
}

void Cylinder::occt_output()
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
		<< "Cylinder update matrix at T=" << tnow
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
