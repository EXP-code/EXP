#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.H"

#include <YamlCheck.H>
#include <EXPException.H>
#include <UserBarTrigger.H>

const std::set<std::string>
UserBarTrigger::valid_keys = {
  "impact",
  "theta",
  "smass",
  "svel",
  "stime",
  "lmax",
  "ctrname"
};

UserBarTrigger::UserBarTrigger(const YAML::Node& conf) : ExternalForce(conf)
{
				// Name for this force (only used for
				// diagnostic output
  id = "BarTrigger";

  // These are all the parameters, clearly, you can add more as
  // necessary . . .

  smass  = 0.1;	     		// Satellite mass
  impact = 0.1;			// Impact parameter
  theta  = 0.0;			// Inclination
  svel   = 1.0;	     		// Satellite velocity
  stime  = -10.0;		// Time origin offset
  lmax   = 2;			// Maximum harmonic order
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
      cerr << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;


  userinfo();

  /*======================================================================
    Put your initialization of the bar trigger force here

    Given that you have a C-like routine, you have two choices:

    1) You can keep your current structure and simply add it
       to this class.  In this case you would split up your
       code into two routines.  The first, call it
       
       my_initialization_routine(...);

       takes all of the input parameters as arguments and does all
       of the computation and caching of tables that you need to
       compute the force.  The force routine itself takes the form:

       void my_force_routine(double tnow, double* pos, double& pot, 
                             double* acc);

       That is, you pass in a pointer to a vector of x, y, z 
       positions ("pos") and the routine returns the potential
       value ("pot") and an acceleration vector ("acc").
       
    1) Alternatively, you can keep your code separate from this
       code by wrapping your existing code in a new class, call it
       
       to this class.  In this case you would split up your
       code into two routines.  The first, call it my_force.
       You would then create an instance of your class as follows:
       
       my_force *model = new my_force(...);

       whose constructor takes all of the input parameters as arguments.
       You would also create a member function 
       
       void
       model->get_pot_acc(double tnow, double* pos, double& pot, double* acc);
       
       which returns the potential value ("pot") and accelration vector
       ("acc) when called with the position vector ("pos") and current time.

    ======================================================================*/

  /*
    my_initialization_routine(lmax, impact, theta, smass, svel, stime);

    // or wrap your code in a class structure to do this, e.g.
  
    model = new my_force(lmax, impact, theta, smass, svel, stime);
  */
}

UserBarTrigger::~UserBarTrigger()
{
  // If you made any new classes, e.g. such as my_new_force, delete them here
  // delete model;
}

void UserBarTrigger::userinfo()
{
  if (myid) return;		// Return if node master node

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  print_divider();

  cout << "** User routine BAR TRIGGER initialized, " ;
  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";

  cout << ", Impact parameter=" << impact << ", Inclination=" << theta
       << ", Satellite mass=" << smass << ", Satellite vel=" << svel 
       << ", Time offset=" << stime << ", Lmax=" << lmax
       << endl;

  print_divider();
}

void UserBarTrigger::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserBarTrigger", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign parameters from YAML config
  //
  try {
    if (conf["impact"])   impact       = conf["impact"].as<double>();
    if (conf["theta"])    theta        = conf["theta"].as<double>();
    if (conf["smass"])    smass        = conf["smass"].as<double>();
    if (conf["svel"])     svel         = conf["svel"].as<double>();
    if (conf["stime"])    stime        = conf["stime"].as<double>();
    if (conf["lmax"])     lmax         = conf["lmax"].as<int>();
    if (conf["ctrname"])  ctr_name     = conf["ctrname"].as<string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserBarTrigger: "
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


void UserBarTrigger::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);

  print_timings("UserBarTrigger: accleration timings");
}


void * UserBarTrigger::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  std::vector<double> pos(3, 0.0), acc(3, 0.0);
  double pot=0.0;

				// This is an iterator for the bodies
  PartMapItr it = cC->Particles().begin();

				// Iterate through the list to get the
				// particles to be computed by this
				// process
  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    unsigned long i = (it++)->first;

				// If we are multistepping, compute accel 
				// only at or below this level
    if (multistep && (cC->Part(i)->level < mlevel)) continue;

				// Computes the position about the center
    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(i, k);	// Inertial by default
      if (c0) pos[k] -= c0->center[k];
    }

    //=======================================================
    // Put your code to compute the bar trigger force here
    //=======================================================

    /*
      my_force_routine(pos, pot, acc);

      // or your class's member function
      
      model->get_pot_acc(pos, pot, acc);
    */
      

    // Add external accerlation
    for (int k=0; k<3; k++)
      cC->AddAcc(i, k, acc[k]);

    // Add external potential
    cC->AddPotExt(i, pot);

  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerBarTrigger(const YAML::Node& conf)
  {
    return new UserBarTrigger(conf);
  }
}

class proxytrigger { 
public:
  proxytrigger()
  {
    factory["usertrigger"] = makerBarTrigger;
  }
};

static proxytrigger p;
