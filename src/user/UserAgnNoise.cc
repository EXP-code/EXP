#include <cmath>
#include <sstream>

#include "expand.H"
#include <localmpi.H>

#include <UserAgnNoise.H>

UserAgnNoise::UserAgnNoise(const YAML::Node &conf) : ExternalForce(conf)
{
  id = "AGNnoise";

  comp_name = "";		// Component for AGN emulation: mandatory

  tau1 = 0.1;			// Default half life
  tau2 = 0.05;			// Default recovery time
  R0   = 0.003;			// Default radius of event region
  eps  = 0.1;			// Default fraction of mass lost in region
  loc  = 0;			// Default location of interaction
				// time in particle attribute array

  initialize();			// Assign parameters

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
  if (myid==0) {
    tev = tnow - tau1*log(number_01(random_gen));
    if (info) std::cout << "**** AGN noise: next event will be T="
			<< tev << std::endl;
  }
  MPI_Bcast(&tev, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
      p->dattrib[loc+0] = p->mass;    // Initial mass
      p->dattrib[loc+1] = -32.0*tau1; // Large negative value
    }
  }
}

UserAgnNoise::~UserAgnNoise()
{
}

void UserAgnNoise::userinfo()
{
  if (myid) return;		// Only root node prints to stdout

  print_divider();

  cout << "** User routine AGN noise initialized "
       << "using component <" << comp_name << "> with" << std::endl
       << "   tau1=" << tau1
       << " tau2="   << tau2
       << " R0="     << R0
       << " eps="    << eps
       << " info="   << std::boolalpha << info
       << std::endl;

  print_divider();
}

void UserAgnNoise::initialize()
{
  try {
    if (conf["compname"]) comp_name = conf["compname"].as<std::string>();
    if (conf["tau1"])     tau1      = conf["tau1"].as<double>();
    if (conf["tau2"])     tau2      = conf["tau2"].as<double>();
    if (conf["R0"])       R0        = conf["R0"].as<double>();
    if (conf["eps"])      eps       = conf["eps"].as<double>();
    if (conf["loc"])      loc       = conf["loc"].as<int>();
    if (conf["info"])     info      = conf["info"].as<bool>();
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
    if (myid==0) {
      tev = tev - tau1*log(number_01(random_gen));
      if (info) std::cout << "**** AGN noise: next event will be T="
			  << tev << std::endl;
    }
    MPI_Bcast(&tev, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
    p->mass = p->dattrib[loc]*(1.0 - eps*exp(-(tnow-p->dattrib[loc+1])/tau2));
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
