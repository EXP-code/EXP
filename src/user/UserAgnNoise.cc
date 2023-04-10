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

  cuda_aware = true;		// Cuda is implemented

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

  // Check attribute count
  //
  if (myid==0) {
    // Sanity check: ensure attribute dimension is large enough
    //
    PartMapItr it = c0->Particles().begin();
    if (it->second->dattrib.size() < loc+1) {
      std::ostringstream sout;
      sout << "Number of Particle float attributes for component ["
	   << c0->name << "] must be >= " << loc+2;
      throw GenericError(sout.str(), __FILE__, __LINE__, 37, false);
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

void UserAgnNoise::setup_decay(void)
{
  int nbodies = cC->Number();
  PartMapItr it = cC->Particles().begin();
  for (int q=0; q<nbodies; q++) {
    PartPtr p = it++->second;
    p->dattrib[loc+0] = p->mass;    // Initial mass
    p->dattrib[loc+1] = -32.0*tau1; // Large negative value
  }
}

void UserAgnNoise::determine_acceleration_and_potential(void)
{
  if (cC != c0) return;		// Check that this component is the target

  // Initialize on first call
  //
  static bool first_call = true;
  if (first_call) {
    
#if HAVE_LIBCUDA==1
    if (use_cuda) {
      setup_decay_cuda();
    } else {
      setup_decay();
    }
#else
    setup_decay();
#endif
    first_call = false;
  }

  // Only compute for top level
  //
  if (multistep && mlevel>0) return;

  // Update masses
  //
#if HAVE_LIBCUDA==1
  if (use_cuda) {
    determine_acceleration_and_potential_cuda();
  } else {
    exp_thread_fork(false);
  }
#else
  exp_thread_fork(false);
#endif

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
      double R2 = 0.0;
      for (int k=0; k<3; k++) R2 += p->pos[k]*p->pos[k];
      if (R2 < R0*R0) {
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