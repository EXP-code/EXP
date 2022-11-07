#include <math.h>
#include <sstream>

#include "expand.H"

#include <UserMNdisk.H>

const std::set<std::string>
UserMNdisk::valid_keys = {
  "ctrname",
  "a",
  "b",
  "mass",
  "Ton",
  "Toff",
  "DeltaT"
};

UserMNdisk::UserMNdisk(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "MiyamotoNagaiDiskPotential";

  a    = 1.0;			// Disk scale length
  b    = 0.1;			// Disk scale height
  mass = 1.0;			// Total disk mass

  Ton = -20.0;			// Turn on start time
  Toff = 200.0;			// Turn off start time
  DeltaT = 1.0;			// Turn on duration

  ctr_name = "";		// Default component for com
  
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
      sout << "Can't find desired component <" << ctr_name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }

  }
  else
    c0 = NULL;


  userinfo();
}

UserMNdisk::~UserMNdisk()
{
}

void UserMNdisk::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine: Miyamoto-Nagai disk with a=" << a 
       << ", b=" << b
       << ", mass=" << mass;

  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";

  cout << endl;

  print_divider();
}

void UserMNdisk::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserMNdisk", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["ctrname"])   ctr_name    = conf["ctrname"].as<string>();
    if (conf["a"])         a           = conf["a"].as<double>();
    if (conf["b"])         b           = conf["b"].as<double>();
    if (conf["mass"])      mass        = conf["mass"].as<double>();
    if (conf["Ton"])       Ton         = conf["Ton"].as<double>();
    if (conf["Toff"])      Toff        = conf["Toff"].as<double>();
    if (conf["DeltaT"])    DeltaT      = conf["DeltaT"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserMNdisk: "
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

void UserMNdisk::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);
}


void * UserMNdisk::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  vector<double> pos(3);

  double amp = 
      0.5*(1.0 + erf( (tnow - Ton )/DeltaT ))
    * 0.5*(1.0 - erf( (tnow - Toff)/DeltaT )) ;


  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;
    

				// Set center if component is
				// defined, otherwise use origin
    if (c0)
      for (int k=0; k<3; k++) 
	pos[k] = cC->Pos(i, k) - c0->center[k];
    else
      for (int k=0; k<3; k++) 
	pos[k] = cC->Pos(i, k);
    
    double xx = pos[0];
    double yy = pos[1];
    double zz = pos[2];

    double rr = sqrt( xx*xx + yy*yy );
    double zb = sqrt( zz*zz + b * b );
    double ab = a + zb;
    double dn = sqrt( rr*rr + ab*ab );
    
    double pot = -mass/dn;
    double fr  = -mass*rr/(dn*dn*dn);
    double fz  = -mass*zz*ab/(zb*dn*dn*dn);

				// Add acceleration by disk
    cC->AddAcc(i, 0, amp * fr*xx/(rr+1.0e-10) );
    cC->AddAcc(i, 1, amp * fr*yy/(rr+1.0e-10) );
    cC->AddAcc(i, 2, amp * fz );

				// Add external potential
    cC->AddPotExt(i, pot);

  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerMNdisk(const YAML::Node& conf)
  {
    return new UserMNdisk(conf);
  }
}

class proxymndisk { 
public:
  proxymndisk()
  {
    factory["usermndisk"] = makerMNdisk;
  }
};

static proxymndisk p;
