#include <cmath>
#include <sstream>

#include "expand.H"
#include <localmpi.H>

#include <UserAgnNoise.H>

UserAgnNoise::UserAgnNoise(const YAML::Node &conf) : ExternalForce(conf)
{
  id = "AGNnoise";

  comp_name = "";		// Component for AGN emulation: mandatory

  tau = 0.1;			// Default half life
  R0  = 0.003;			// Default radius of event region
  eps = 0.1;			// Default fraction of mass lost in region
  loc = 0;			// Default location of interaction
				// time in particle attribute array

  initialize();

  if (comp_name.size()>0) {
				// Look for the fiducial component
    bool found = false;
    for (auto c : comp->components) {
      if ( !comp_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      std::ostringstream sout;
      sout << "Can't find desired component <" << comp_name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }

  } else {
    std::string msg = "UserAgnNoise: desired component name must be specified";
    throw GenericError(msg, __FILE__, __LINE__, 36, false);
  }

  userinfo();

  // Assign uniform generator for Poisson events
  //
  number_01 = std::uniform_real_distribution<double>(0.0, 1.0);

  // Compute time for first event
  //
  tev = tnow - tau*log(number_01(random_gen));

  // Assign initial counter values
  //
  if (not restart) {
    int nbodies = c0->Number();
    PartMapItr it = c0->Particles().begin();
    if (myid==0) {
      // Sanity check: attribute dimension
      if (it->second->dattrib.size() < loc+1) {
	std::ostringstream sout;
	sout << "Number float attribute in Particle must be >= " << loc+2;
	throw GenericError(sout.str(), __FILE__, __LINE__, 37, false);
      }
    }
    // Assign initial attribute values
    for (int q=0; q<nbodies; q++) {
      PartPtr p = it++->second;
      p->dattrib[loc+0] = p->mass;   // Initial mass
      p->dattrib[loc+1] = -32.0*tau; // Large negative value
    }
  }
}

UserAgnNoise::~UserAgnNoise()
{
}

void UserAgnNoise::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine AGN noise initialized "
       << "using component <" << comp_name << "> with" << std::endl
       << "    tau=" << tau
       << " R0="  << R0
       << " eps=" << eps
       << std::endl;

  print_divider();
}

void UserAgnNoise::initialize()
{
  try {
    if (conf["compname"]) comp_name = conf["compname"].as<std::string>();
    if (conf["tau"])      tau       = conf["tau"].as<double>();
    if (conf["R0"])       R0        = conf["R0"].as<double>();
    if (conf["esp"])      eps       = conf["eps"].as<double>();
    if (conf["loc"])      loc       = conf["loc"].as<int>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserAgnNoise: "
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


void UserAgnNoise::determine_acceleration_and_potential(void)
{
  if (cC != c0) return;		// Check that this component is the target

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  // Only compute for top level
  //
  if (multistep && mlevel>0) return;

  // Update masses
  //
  exp_thread_fork(false);

  // Compute next event
  //
  if (tnow > tev) {
    tev = tev - tau*log(number_01(random_gen));
  }

}


void * UserAgnNoise::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = c0->Number();
  int id      = *((int*)arg);
  int nbeg    = nbodies*id/nthrds;
  int nend    = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  // Move iterator to current thread batch
  //
  PartMapItr it = cC->Particles().begin();
  for (int q=0   ; q<nbeg; q++) it++;

  // Update masses and times for this batch
  //
  for (int q=nbeg; q<nend; q++) {
    PartPtr p = it++->second;	// Particle pointer for convenience

    if (tnow>tev) {		// Did an event occur?
      double R = 0.0;
      for (int k=0; k<3; k++) R += p->pos[k]*p->pos[k];
      if (R*R < R0*R0) {
	p->dattrib[loc+1] = tnow;
      }
    }
				// The updated particle mass
    p->mass = p->dattrib[loc]*(1.0 - eps*exp(-(tnow-p->dattrib[loc+1])/tau));
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerAgnNoise(const YAML::Node& conf)
  {
    return new UserAgnNoise(conf);
  }
}

class proxyagnnoise { 
public:
  proxyagnnoise()
  {
    factory["useragnnoise"] = makerAgnNoise;
  }
};

static proxyagnnoise p;
