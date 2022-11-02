#include <math.h>
#include <sstream>

#include "expand.H"

#include <UserLogPot.H>

const std::set<std::string>
UserLogPot::valid_keys = {
  "R", 
  "b",
  "c",
  "v2"
};

UserLogPot::UserLogPot(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "LogarithmicPotential";

  initialize();

  userinfo();
}

UserLogPot::~UserLogPot()
{
}

void UserLogPot::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine LOGARITHMIC POTENTIAL initialized, " ;
  cout << "Phi = v2/2*log(R^2 + x^2 + y^2/b^2 + z^2/c^2) with R=" 
       << R << ", b=" << b << ", c=" << c << ", v2=" << v2 << endl; 
  
  print_divider();
}

void UserLogPot::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserPeriodic", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["R"])      R     = conf["R"].as<double>();
    if (conf["b"])      b     = conf["b"].as<double>();
    if (conf["c"])      c     = conf["c"].as<double>();
    if (conf["v2"])     v2    = conf["v2"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserLogPot: "
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


void UserLogPot::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);
}


void * UserLogPot::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double xx, yy, zz, rr;

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    xx = cC->Pos(i, 0);
    yy = cC->Pos(i, 1);
    zz = cC->Pos(i, 2);
    rr = R*R + xx*xx + yy*yy/(b*b) + zz*zz/(c*c);

    cC->AddAcc(i, 0,-v2*xx/rr );
    
    cC->AddAcc(i, 1, -v2*yy/(rr*b*b) );

    cC->AddAcc(i, 2, -v2*zz/(rr*c*c) );
    
    cC->AddPotExt(i, 0.5*v2*log(rr) );
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerLogPot(const YAML::Node& conf)
  {
    return new UserLogPot(conf);
  }
}

class proxylogpot { 
public:
  proxylogpot()
  {
    factory["userlogp"] = makerLogPot;
  }
};

static proxylogpot p;
