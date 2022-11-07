#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.H"

#include <UserSlabHalo.H>

const std::set<std::string>
UserSlabHalo::valid_keys = {
  "ctrname",
  "h0",
  "z0",
  "rho0",
  "v0"
};

UserSlabHalo::UserSlabHalo(const YAML::Node& conf) : ExternalForce(conf)
{

  id   = "SlabHalo";		// Halo model file

  rho0 = 1.0;			// Central density
  h0   = 1.0;			// Scale height
  U0   = 4.0*M_PI*rho0*h0*h0;	// Potential coefficient

  z0   = 0.5;			// Position of the midplane

  v0   = sqrt(0.5*U0);	        // Isothermal velocity dispersion

  ctr_name = "";		// Default component for center (none)

  initialize();

  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    for (auto c : comp->components) {
      if ( !ctr_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      std::ostringstream sout;
      sout << "libslabhalo, process " << myid 
	   << ": can't find desired component <"
	   << ctr_name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }

  }
  else
    c0 = NULL;


  userinfo();
}

UserSlabHalo::~UserSlabHalo()
{
}

void UserSlabHalo::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine SLAB HALO initialized, v0=" << v0 
       << ", rho0=" << rho0 << ", h0=" << h0 << ", z0=" << z0
       << ", U0=" << U0;
  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";
  cout << endl;

  print_divider();
}

void UserSlabHalo::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserSlabHalo", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["ctrname"])        ctr_name           = conf["ctrname"].as<string>();
    if (conf["h0"])             h0                 = conf["h0"].as<double>();
    if (conf["z0"])             z0                 = conf["z0"].as<double>();
    
    if (conf["rho0"]) {
      rho0 = conf["rho0"].as<double>();
      U0  = 4.0*M_PI*rho0*h0*h0;
      v0 = sqrt(0.5*U0);
    }
    
    if (conf["v0"]) {
      v0 = conf["v0"].as<double>();
      U0 = 2.0*v0*v0;
      rho0 = U0/(4.0*M_PI*h0*h0);
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserSlabHalo: "
			   << error.what()         << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}


void UserSlabHalo::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);
  
  print_timings("UserSlabHalo: accleration timings");
}


void * UserSlabHalo::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  
  thread_timing_beg(id);
  
  double pos[3];
  
  PartMapItr it = cC->Particles().begin();
  
  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    unsigned long i = (it++)->first;
    // If we are multistepping, compute accel 
    // only at or below this level
    if (multistep && (cC->Part(i)->level < mlevel)) continue;
    
    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(i, k);	// Inertial by default
      if (c0) pos[k] -= c0->center[k];
    }
    
    // Add external acceleration
    cC->AddAcc(i, 2, -U0/h0*tanh((pos[2]-z0)/h0));
    
    // Add external potential
    cC->AddPotExt(i, U0*log(cosh((pos[2]-z0)/h0)) );
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerSlabHalo(const YAML::Node& conf)
  {
    return new UserSlabHalo(conf);
  }
}

class proxyslabhalo { 
public:
  proxyslabhalo()
  {
    factory["userslabhalo"] = makerSlabHalo;
  }
};

static proxyslabhalo p;
