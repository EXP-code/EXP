#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.H"

#include "PeriodicBC.H"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

const std::set<std::string>
PeriodicBC::valid_keys = {
  "compname",
  "sx",
  "sy",
  "sz",
  "cx",
  "cy",
  "cz",
  "btype"
};

PeriodicBC::PeriodicBC(const YAML::Node& conf) : ExternalForce(conf)
{
  (*barrier)("PeriodicBC: BEGIN construction", __FILE__, __LINE__);

  id = "PeriodicBC";		// Periodic boundary condition ID

				// Sizes in each dimension
  L = vector<double>(3, 1.0);	
				// Center offset in each dimension
  offset = vector<double>(3, 0.0);

  bc = "ppp";			// Periodic BC in all dimensions

  comp_name = "";		// Default component (must be specified)

  initialize();			// Assign parameters

  cuda_aware = true;		// Cuda routine is implemented

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
    sout << "PeriodicBC: "
	 << "can't find fiducial component <" << comp_name << ">";
    throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
  }
  
  userinfo();

  (*barrier)("Periodic: END construction", __FILE__, __LINE__);
}

PeriodicBC::~PeriodicBC()
{
}

void PeriodicBC::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** External routine PERIODIC BOUNDARY CONDITION initialized"
       << " using component <" << comp_name << ">" << endl;

  cout << "   Cube sides (x , y , z) = (" 
       << L[0] << " , " 
       << L[1] << " , " 
       << L[2] << " ) " << endl; 

  cout << "Center offset (x , y , z) = (" 
       << offset[0] << " , " 
       << offset[1] << " , " 
       << offset[2] << " ) " << endl; 

  cout << "Boundary type (x , y , z) = (" 
       << bc[0] << " , " 
       << bc[1] << " , " 
       << bc[2] << " ) " << endl;

  print_divider();
}

void PeriodicBC::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("PeriodicBC", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["compname"])       comp_name          = conf["compname"].as<string>();
    
    if (conf["sx"])             L[0]               = conf["sx"].as<double>();
    if (conf["sy"])             L[1]               = conf["sy"].as<double>();
    if (conf["sz"])             L[2]               = conf["sz"].as<double>();
    
    if (conf["cx"])             offset[0]          = conf["cx"].as<double>();
    if (conf["cy"])             offset[1]          = conf["cy"].as<double>();
    if (conf["cz"])             offset[2]          = conf["cz"].as<double>();
    
    if (conf["btype"]) {
      std::string val = conf["btype"].as<std::string>();
      if (strlen(val.c_str()) >= 3) {
	for (int k=0; k<3; k++) {
	  switch (val.c_str()[k]) {
	  case 'p':
	    bc[k] = 'p';		// Periodic
	    break;
	  case 'r':
	    bc[k] = 'r';		// Reflection
	    break;
	  default:
	    bc[k] = 'v';		// Vacuum
	    break;
	  }
	}
      }
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in PeriodicBC: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("PeriodicBC::initialize: error parsing YAML");
  }
    

  // Cuda initialization
#if HAVE_LIBCUDA==1
  if (use_cuda) {
    cuda_initialize();
  }
#endif
}


void PeriodicBC::determine_acceleration_and_potential(void)
{
  if (cC != c0) return;
#if HAVE_LIBCUDA==1		// Cuda compatibility
  if (use_cuda) {
    determine_acceleration_and_potential_cuda();
  } else
  getParticlesCuda(cC);
  exp_thread_fork(false);
#endif
  exp_thread_fork(false);
  
  print_timings("PeriodicBC: thread timings");
}

void * PeriodicBC::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  
  thread_timing_beg(id);

  double pos, delta;
  PartMapItr it = cC->Particles().begin();
  
  std::advance(it, nbeg);

  for (int q=nbeg; q<nend; q++) {
    
				// Index for the current particle
    unsigned long i = (it++)->first;
    
    Particle *p = cC->Part(i);
    double   mi = 0.0;

    for (int k=0; k<3; k++) {

      // Ignore vacuum boundary dimensions
      //
      if (bc[k] == 'v') continue;

      // Increment so that the positions range
      // between 0 and L[k]
      //
      pos = p->pos[k] + offset[k];

      //
      // Reflection BC
      //
      if (bc[k] == 'r') {
	if (pos < 0.0) {
	  delta = -pos - L[k]*floor(-pos/L[k]);
	  p->pos[k] = delta - offset[k];
	  p->vel[k] *= -1.0;
	} 
	if (pos >= L[k]) {
	  delta = pos - L[k]*floor(pos/L[k]);
	  p->pos[k] =  L[k] - delta - offset[k];
	  p->vel[k] *= -1.0;
	}
      }

      //
      // Periodic BC
      //
      if (bc[k] == 'p') {
	if (pos < 0.0) {
	  p->pos[k] += L[k]*floor(1.0+fabs(pos/L[k]));
	  
	}
	if (pos >= L[k]) {
	  p->pos[k] += - L[k]*floor(fabs(pos/L[k]));
	}
      }

      //
      // Sanity check
      //
      if (p->pos[k] < -offset[k] || p->pos[k] >= L[k]-offset[k]) {
	cout << "Process " << myid << " id=" << id 
	     << ": Error in pos[" << k << "]=" << p->pos[k] << endl;
      }
    }
  }
  
  thread_timing_end(id);

  return (NULL);
}

