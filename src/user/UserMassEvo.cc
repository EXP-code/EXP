#include <math.h>
#include <sstream>

#include "expand.H"
#include <localmpi.H>
#include <gaussQ.H>

#include <UserMassEvo.H>

UserMassEvo::UserMassEvo(const YAML::Node &conf) : ExternalForce(conf)
{
  id = "MassEvo";

  comp_name = "";		// Default component for com
  t0 = 0.0;			// Default time offset
  a0 = 1.0;			// Default constant term
  a1 = 0.0;			// Default linear term
  a2 = 0.0;			// Default quadratic term
  a3 = 0.0;			// Default cubic term

  initialize();

  last = fct(tnow);		// Cache value for initial time

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
      cerr << "Process " << myid << ": can't find desired component <"
	   << comp_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  } else {
    if (myid==0) {
      std:: cerr << "UserMassEvo: desired component name must be specified"
	<< std::endl;
	   MPI_Abort(MPI_COMM_WORLD, 36);
    }
  }

  userinfo();
}

UserMassEvo::~UserMassEvo()
{
}

void UserMassEvo::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine MASS EVOLUTION initialized "
       << "using component <" << comp_name << "> with"
       << " t0=" << t0
       << " a0=" << a0
       << " a1=" << a1
       << " a2=" << a2
       << " a3=" << a3
       << std::endl;

  print_divider();
}

void UserMassEvo::initialize()
{
  try {
    if (conf["compname"]) comp_name = conf["compname"].as<std::string>();
    if (conf["t0"])       t0        = conf["t0"].as<double>();
    if (conf["a0"])       a0        = conf["a0"].as<double>();
    if (conf["a1"])       a1        = conf["a1"].as<double>();
    if (conf["a2"])       a2        = conf["a2"].as<double>();
    if (conf["a3"])       a3        = conf["a3"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserMassEvo: "
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


void UserMassEvo::determine_acceleration_and_potential(void)
{
  if (cC != c0) return;		// Check that this component is the target

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif
				// Only compute for top level
  if (multistep && mlevel>0) return;

  curr = fct(tnow);

  exp_thread_fork(false);

  last = curr;
}


void * UserMassEvo::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = c0->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  PartMapItr it = cC->Particles().begin();

  for (int q=0   ; q<nbeg; q++) it++;

  for (int q=nbeg; q<nend; q++) {
    unsigned long i = (it++)->first;
    c0->Part(i)->mass *= curr/last;
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerMassEvo(const YAML::Node& conf)
  {
    return new UserMassEvo(conf);
  }
}

class proxymassevo { 
public:
  proxymassevo()
  {
    factory["usermassevo"] = makerMassEvo;
  }
};

static proxymassevo p;
