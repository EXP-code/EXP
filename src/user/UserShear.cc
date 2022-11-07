#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.H"

#include <UserShear.H>


UserShear::UserShear(const YAML::Node& conf) : ExternalForce(conf)
{

  id       = "ShearingSheet";	// Shearing sheet
  r0       = 1;			// R_o in radial units
  s0       = 1;			// Sigma in velocity units
  xoff     = 0.5;	        // Offset in x-direction for shear center
  ctr_name = "";		// Default component for center

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
      sout << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }

  }
  else
    c0 = NULL;

  userinfo();
}

UserShear::~UserShear()
{
}

void UserShear::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine SHEARING SHEET initialized, " ;
  cout << "r_0=" << r0 << "  sigma=" << s0 << "  v_0=" << sqrt(2.0)*s0
       << "  offset=" << xoff;
  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";
  cout << endl;

  print_divider();
}

const std::set<std::string>
UserShear::valid_keys = {
  "radius",
  "velocity",
  "offset",
  "ctrname"
};

void UserShear::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserShear", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["radius"])    r0           = conf["radius"].as<double>();
    if (conf["velocity"])  s0           = conf["velocity"].as<double>();
    if (conf["offset"])    xoff         = conf["offset"].as<double>();
    if (conf["ctrname"])   ctr_name     = conf["ctrname"].as<string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserShear: "
			   << error.what()         << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }


  omega = sqrt(2.0)*s0/r0;
  kappa = 2.0*s0/r0;
  A     =  s0/(sqrt(2.0)*r0);
  B     = -s0/(sqrt(2.0)*r0);
}


void UserShear::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);

  print_timings("UserShear: accleration timings");
}


void * UserShear::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  double x, xd, yd, ax, ay, pot;

  PartMapItr it = cC->Particles().begin();

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    unsigned long i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    x  = cC->Pos(i, 0) - xoff;
    xd = cC->Vel(i, 0);
    yd = cC->Vel(i, 1);

    ax =  2.0*omega*yd + 4.0*A*omega*x;
    ay = -2.0*omega*xd;

    cC->AddAcc(i, 0, ax);
    cC->AddAcc(i, 1, ay);

    // Add external potential
    pot = -2.0*A*omega*x*x;
    cC->AddPotExt(i, pot);
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerShear(const YAML::Node& conf)
  {
    return new UserShear(conf);
  }
}

class proxyhalo { 
public:
  proxyhalo()
  {
    factory["usershear"] = makerShear;
  }
};

static proxyhalo p;
