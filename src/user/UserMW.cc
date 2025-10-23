#include <math.h>
#include <sstream>

#include "expand.H"

#include "UserMW.H"

UserMW::UserMW(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "MilkyWayPotential";

  // !!! add all parameter defaults here (must be defined in UserMW.H) !!! Done

  G = 4.3e-6; // The gravitational constant

  M_halo = 5.4e11; // Total mass of the Navarro-Frank-White (NFW) halo component
  rs_halo = 15.62; // Scale radius of the Navarro-Frank-White (NFW) halo component

  M_disk = 6.8e10; // Total mass of the Miyamoto-Nagai (MN) disk component
  a_disk = 3.0; // Scale length of the Miyamoto-Nagai (MN) disk component
  b_disk = 0.28; // Scale height of the Miyamoto-Nagai (MN) disk component

  M_nucl = 1.71e9; // Total mass of the Hernquist (HN) nucleus component
  c_nucl = 0.07; // Concentration of the Hernquist (HN) nucleus component

  M_bulge = 5e9; // Total mass of the Hernquist (HN) bulge component
  c_bulge = 1.0; // Concentration of the Hernquist (HN) bulge component

  // these are BOILERPLATE
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
      cerr << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;


  userinfo();
}

UserMW::~UserMW()
{
}

void UserMW::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  // !!! Fill out this part with the different parameters.
  //   this gets printed into the log file, so we can use as a diagnostic down the road !!! Done

  cout << "** User routine: Milky Way potential with "
       << "M_halo="      << M_halo
       << ", rs_halo="    << rs_halo
       << ", M_disk=" << M_disk
       << ", a_disk=" << a_disk
       << ", b_disk=" << b_disk
       << ", M_nucl=" << M_nucl
       << ", c_nucl=" << c_nucl
       << ", M_bulge=" << M_bulge
       << ", c_bulge=" << c_bulge;

  // BOILERPLATE
  if (c0)
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";

  cout << endl;

  print_divider();
}

void UserMW::initialize()
{
  // Fill out with all variable names
  try {
    if (conf["Mhalo"])         M_halo   = conf["Mhalo"].as<double>(); // Total mass of the NFW halo component component
    if (conf["rshalo"])        rs_halo  = conf["rshalo"].as<double>(); // Scale radius of the Navarro-Frank-White (NFW) halo component

    if (conf["Mdisk"])         M_disk   = conf["Mdisk"].as<double>(); // Total mass of the Miyamoto-Nagai (MN) disk component
    if (conf["adisk"])         a_disk   = conf["adisk"].as<double>(); // Scale length of the Miyamoto-Nagai (MN) disk component
    if (conf["bdisk"])         b_disk   = conf["bdisk"].as<double>(); // Scale height of the Miyamoto-Nagai (MN) disk component

    if (conf["Mnucl"])         M_nucl   = conf["Mnucl"].as<double>(); // Total mass of the Hernquist (HN) nucleus component
    if (conf["cnucl"])         c_nucl   = conf["cnucl"].as<double>(); // Concentration of the Hernquist (HN) nucleus component

    if (conf["Mbulge"])         M_bulge   = conf["Mbulge"].as<double>(); // Total mass of the Hernquist (HN) bulge component
    if (conf["cbulge"])         c_bulge   = conf["cbulge"].as<double>(); // Concentration of the Hernquist (HN) bulge component

    if (conf["ctrname"])   ctr_name    = conf["ctrname"].as<string>();// BOILERPLATE
    if (conf["Ton"])       Ton         = conf["Ton"].as<double>();    // BOILERPLATE
    if (conf["Toff"])      Toff        = conf["Toff"].as<double>();   // BOILERPLATE
    if (conf["DeltaT"])    DeltaT      = conf["DeltaT"].as<double>(); // BOILERPLATE
  }
  catch (YAML::Exception & error) { // All BOILERPLATE
    if (myid==0) std::cout << "Error parsing parameters in UserMW: "
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

// BOILERPLATE
void UserMW::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);
}


void * UserMW::determine_acceleration_and_potential_thread(void * arg)
{
  // BOILERPLATE
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  // Allocate particle workspace
  vector<double> pos(3);

  // This is relevant if we are turning on or off the MW.
  double amp =
      0.5*(1.0 + erf( (tnow - Ton )/DeltaT ))
    * 0.5*(1.0 - erf( (tnow - Toff)/DeltaT )) ;

  // Make particle map to loop through
  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  // Loop through particles
  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;

    // BOILERPLATE
    // If we are multistepping, compute accel
		// only at or below this level
    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    // BOILERPLATE
		// Set center if component is
		// defined, otherwise use origin
    if (c0)
      for (int k=0; k<3; k++)
	      pos[k] = cC->Pos(i, k) - c0->center[k];
    else
      for (int k=0; k<3; k++)
	      pos[k] = cC->Pos(i, k);

    // Extract particle positions
    double xx = pos[0];
    double yy = pos[1];
    double zz = pos[2];

    // Make 2d and 3d radii
    double r2 = sqrt( xx*xx + yy*yy ); // R
    double r3 = sqrt( xx*xx + yy*yy + zz*zz ); // r

    // !!! Fill in with TOTAL potential and force (accelerations) !!! Done
    // !!! definitions for forces can be added in UserMW.H for cleanliness !!!
    // !!! can call above, or call in-line, your choice !!!

    double ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, ax4, ay4, az4; // Is there a way to make this neater?

    NFW_dphi_dr(xx, yy, zz, &ax1, &ay1, &az1);
    HN_nucl_dphi_dr(xx, yy, zz, &ax2, &ay2, &az2);
    HN_bulge_dphi_dr(xx, yy, zz, &ax3, &ay3, &az3);
    MN_dphi_dR_dz(xx, yy, zz, &ax4, &ay4, &az4);

    double pot = NFW_pot(r3) + HN_nucl_pot(r3) + HN_bulge_pot(r3) + MN_pot(zz, r2);
    double xacc = ax1 + ax2 + ax3 + ax4;
    double yacc = ay1 + ay2 + ay3 + ay4;
    double zacc = az1 + az2 + az3 + az4;

		// Add acceleration to particle
    cC->AddAccExt(i, 0, xacc );
    cC->AddAccExt(i, 1, yacc );
    cC->AddAccExt(i, 2, zacc );

		// Add external potential
    cC->AddPotExt(i, pot);

  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerMW(const YAML::Node& conf)
  {
    return new UserMW(conf);
  }
}

class proxymw {
public:
  proxymw()
  {
    factory["usermw"] = makerMW;
  }
};

static proxymw p;
