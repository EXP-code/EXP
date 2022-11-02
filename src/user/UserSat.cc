#include <UserSat.H>

int UserSat::instances = 0;

const std::set<std::string>
UserSat::valid_keys = {
  "comname",
  "config",
  "core",
  "mass",
  "ton",
  "toff",
  "delta",
  "toffset",
  "orbit",
  "shadow",
  "verbose",
  "r0",
  "phase",
  "omega",
  "trajtype"
};

UserSat::UserSat(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "UserSat";		// ID string

  core = 0.5;			// Satellite core size
  mass = 0.3;			// Satellite mass
  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration
  toffset = 0.0;		// Time offset for orbit

  orbit = false;		// Print out orbit for debugging
  shadow = false;		// Simulate an inverse symmetric satellite
  verbose = false;		// Print messages on zero particles
  r0 = 1.0;			// Radius
  phase = 0.0;			// Initial position angle
  omega = 1.0;			// Angular frequency

  pinning  = false;	        // Pin reference frame to a component
  com_name = "";		// Default component for com
  traj_type = circ;		// Trajectory type (default is circ)

  initialize();

  if (pinning) {		// Look for the fiducial component
    bool found = false;
    for (auto c : comp->components) {
      if ( !com_name.compare(c->name) ) {
	c0 = c;
	found = true;
	break;
      }
    }
    
    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << com_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }
  }

  switch (traj_type) {
  case circ:
    traj = 0;
    break;
  case bound:
    traj = std::make_shared<SatelliteOrbit>(config);
    break;
  case unbound:
    traj = std::make_shared<UnboundOrbit>(config);
    break;
  case linear:
    traj = std::make_shared<LinearOrbit>(config);
    break;
  default:
    if (myid==0) {
      cerr << "UserSat: no such trjectory type="
	   << traj_type << endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 36);
  }

  if (orbit && myid==0) {
    std::ostringstream sout;
    sout << outdir << "UserSat." << runtag << "." << ++instances;
    orbfile = sout.str();
    std::ofstream out (orbfile);
    if (out) {
      out << left << setfill('-')
	  << setw(15) << "#"
	  << setw(15) << "+"
	  << setw(15) << "+"
	  << setw(15) << "+"
	  << setw(15) << "+"
	  << endl << setfill(' ')
	  << setw(15) << "# Time"
	  << setw(15) << "+ Mass"
	  << setw(15) << "+ X-pos"
	  << setw(15) << "+ Y-pos"
	  << setw(15) << "+ Z-pos"
	  << endl << setfill('-')
	  << setw(15) << "#"
	  << setw(15) << "+"
	  << setw(15) << "+"
	  << setw(15) << "+"
	  << setw(15) << "+"
	  << endl << setfill(' ');
    } else {
      std::cerr << "UserSat: could not open orbit diagnostic file <"
		<< orbfile << ">" << std::endl;
    }
      
    tlast = tnow;
  }

  zbflag = true;

  userinfo();
}

UserSat::~UserSat()
{
  // Nothing
}

void UserSat::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine SATELLITE IN FIXED POTENTIAL initialized using ";
  if (traj_type==circ)
    cout << "fixed circular orbit with mass=" << mass 
	 << ", core=" << core
	 << ", r=" << r0 
	 << ", p(0)=" << phase 
	 << ", Omega=" << omega;
  else {
    cout << "specified trjectory of type ";
    switch (traj_type) {
    case bound:
      cout << "<BOUND>";
      break;
    case unbound:
      cout << "<UNBOUND>";
      break;
    case linear:
      cout << "<LINEAR>";
      break;
    default:
      cout << "<UNKNOWN>";
      break;
    }
    cout << " with mass=" << mass 
	 << ", core=" << core
	 << ", config=" << config;
    if (pinning) cout << ", centered on Component <" << c0->name << ">";
  }
  if (shadow)  cout << ", shadowing is on";
  if (verbose) cout << ", verbose messages are on";
  if (orbit)   cout << ", with trajectory logging";

  cout << endl;

  print_divider();
}

void UserSat::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserSat", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["comname"]) {
      com_name = conf["comname"].as<std::string>();
      pinning = true;
    }
    
    if (conf["config"])         config             = conf["config"];
    if (conf["core"])           core               = conf["core"].as<double>();
    if (conf["mass"])           mass               = conf["mass"].as<double>();
    if (conf["ton"])            ton                = conf["ton"].as<double>();
    if (conf["toff"])           toff               = conf["toff"].as<double>();
    if (conf["delta"])          delta              = conf["delta"].as<double>();
    if (conf["toffset"])        toffset            = conf["toffset"].as<double>();
    if (conf["orbit"])          orbit              = conf["orbit"].as<bool>();
    if (conf["shadow"])         shadow             = conf["shadow"].as<bool>();
    if (conf["verbose"])        verbose            = conf["verbose"].as<bool>();
    if (conf["r0"])             r0                 = conf["r0"].as<double>();
    if (conf["phase"])          phase              = conf["phase"].as<double>();
    if (conf["omega"])          omega              = conf["omega"].as<double>();
    
    // Set trajectory type
    if (conf["trajtype"]) {
      std::string val = conf["trajtype"].as<std::string>();
      switch (atoi(val.c_str())) {
      case circ:
	traj_type = circ;
	break;
      case bound:
	traj_type = bound;
	break;
      case unbound:
	traj_type = unbound;
	break;
      case linear:
	traj_type = linear;
	break;
      default:
	if (myid==0) {
	  cerr << "UserSat: no such trjectory type="
	       << val << endl;
	}
	MPI_Abort(MPI_COMM_WORLD, 36);
      }
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserSat: "
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


void UserSat::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1
  if (use_cuda) {
    determine_acceration_and_potential_cuda();
  } else {
    exp_thread_fork(false);
  }
#else
  exp_thread_fork(false);
#endif
}

void * UserSat::determine_acceleration_and_potential_thread(void * arg) 
{
  double pos[3], rs[3], fac, ffac, phi;
  double satmass;
  // Sanity check
  int nbodies = cC->Number();
  if (nbodies != static_cast<int>(cC->Particles().size())) {
    cerr << "UserSat: ooops! number=" << nbodies
	 << " but particle size=" << cC->Particles().size() << endl;
    nbodies = static_cast<int>(cC->Particles().size());
  }
  
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  if (nbodies==0) {		// Return if there are no particles
    if (verbose) {
      if (id==0 && zbflag) {	// Only print message on state change
	cout << "Process " << myid << ": in UserSat, nbodies=0" 
	     << " for Component <" << cC->name << "> at T=" << tnow
	     << endl;
	zbflag = false;
      }
    }
    thread_timing_end(id);
    return (NULL);
  }

  if (id==0) zbflag = true;

  if (traj_type==circ) {
    phi = phase + omega*tnow;
    rs[0] = r0*cos(phi);
    rs[1] = r0*sin(phi);
    rs[2] = 0.0;
  }
  else
    traj->get_satellite_orbit(tnow - toffset, rs);

  satmass = mass * 
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    
  if (shadow) satmass *= 0.5;

  if (orbit && myid==0 && id==0 && mlevel==0 && tnow>tlast) {
    ofstream out (orbfile.c_str(), ios::app);
    if (out) {
      out << setw(15) << tnow << setw(15) << satmass;
      for (int k=0; k<3; k++) out << setw(15) << rs[k];
      out << endl;
      tlast = tnow;
    } else {
      cout << "Error opening trajectory file: " << orbfile << endl;
    }
  }


  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    fac = core*core;

    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(i, k, Component::Inertial);
      if (pinning) pos[k] -= c0->com[k];
      fac += (pos[k] - rs[k])*(pos[k] - rs[k]);
    }
    fac = pow(fac, -0.5);
    
    ffac = -satmass*fac*fac*fac;

    // Add acceration
    for (int k=0; k<3; k++) cC->AddAcc(i, k, ffac*(pos[k]-rs[k]) );
    
    // Add external potential
    cC->AddPotExt(i, -satmass*fac );

    // Add the shadow satellite
    if (shadow) {

      fac = core*core;
      for (int k=0; k<3; k++)
	fac += (pos[k] + rs[k])*(pos[k] + rs[k]);
      fac = pow(fac, -0.5);
    
      ffac = -satmass*fac*fac*fac;

      // Add acceration
      for (int k=0; k<3; k++) cC->AddAcc(i, k, ffac*(pos[k]+rs[k]) );
    
      // Add external potential
      cC->AddPotExt(i, -satmass*fac );
    }

  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerSat(const YAML::Node& conf)
  {
    return new UserSat(conf);
  }
}

class proxysat { 
public:
  proxysat()
  {
    factory["usersat"] = makerSat;
  }
};

static proxysat p;
