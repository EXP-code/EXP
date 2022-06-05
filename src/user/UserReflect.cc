#include <sys/timeb.h>
#include <math.h>
#include <sstream>
#include <random>

#include "expand.H"

#include <UserReflect.H>


UserReflect::UserReflect(const YAML::Node& conf) : ExternalForce(conf)
{

  id = "ReflectBC";		// Reflect boundary condition ID

  comp_name = "";		// Default component (must be specified)

  radius = 1.0;			// Radius for the reflective boundary

  debug = false;		// Debug info

  initialize();

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
    cerr << "UserReflect: process " << myid 
	 << " can't find fiducial component <" << comp_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }
  
  if (debug) {
    wrong_dir = new unsigned [nthrds];
    too_big   = new unsigned [nthrds];
  } else {
    wrong_dir = 0;
    too_big   = 0;
  }

  // Set random number seeds for procs and threads
  for (int n=0; n<nthrds; n++) gen[n].seed(11 + nthrds*myid + n);

  userinfo();
}

UserReflect::~UserReflect()
{
  delete [] wrong_dir;
  delete [] too_big;
}

void UserReflect::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine REFLECTIVE SPHERICAL BOUNDARY CONDITION initialized"
       << " using component <" << comp_name << "> with radius=" << radius;
  if (debug) cout << ", with debug output ON";
  cout << endl;

  print_divider();
}

void UserReflect::initialize()
{
  try {
    if (conf["compname"])  comp_name   = conf["compname"].as<string>();
    if (conf["radius"])    radius      = conf["radius"].as<double>();
    if (conf["debug"])     debug       = conf["debug"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserReflect: "
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


void UserReflect::determine_acceleration_and_potential(void)
{
  if (cC != c0) return;

  if (debug) {
    for (int n=0; n<nthrds; n++) wrong_dir[n] = too_big[n] = 0;
  }

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);

  if (debug) {
    for (int n=1; n<nthrds; n++) {
      wrong_dir[0] += wrong_dir[n];
      too_big[0]   += too_big[n];
    }
    unsigned total_dir=0, total_big=0;
    MPI_Reduce(wrong_dir, &total_dir, 1, MPI_UNSIGNED, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(too_big,    &total_big, 1, MPI_UNSIGNED, MPI_SUM, 0,
	       MPI_COMM_WORLD);

    if (myid==0) {
      if (total_dir)
	cout << "UserReflect [" << cC->name <<"]: found "
	     << total_dir << " in opposite direction" << endl;
      if (total_big)
	cout << "UserReflect [" << cC->name <<"]: found "
	     << total_big << " particles *way* beyond boundary" << endl;
    }
  }

}


void * UserReflect::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  
  vector<double> pos(3), vel(3), dif(3);
  PartMapItr it = cC->Particles().begin();
  double rr, rv, delr;
  
  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    
				// Index for the current particle
    unsigned long i = (it++)->first;
    
				// If we are multistepping, compute BC
				// only at or above this level
    if (!multistep || (cC->Part(i)->level >= mlevel)) {


				// For access to the particle position
				// and velocity
      Particle *p = cC->Part(i);


				// Shift the coordinate center
      for (int k=0; k<3; k++) {
	pos[k] = p->pos[k] - cC->center[k];
	if (cC->com_system) pos[k] -= cC->com0[k];
      }

				// Compute the radius
      rr = 0.0;
      for (int k=0; k<3; k++) rr += pos[k]*pos[k];
      rr = sqrt(rr);
      
				// If inside the boundary, next particle
      if (rr < radius) continue;
      
				// Compute radial projection of velocity

				// Shift the velocity center
      for (int k=0; k<3; k++) {
	vel[k] = p->vel[k];
	if (cC->com_system) vel[k] -= cC->cov0[k];
      }

      rv = 0.0;
      for (int k=0; k<3; k++) {
	rv += pos[k]/rr * vel[k];
      }
				// rv is the radial velocity
				// Particle should be outgoing
      if (rv > 0.0) {
				// Reverse the component of radial velocity
				//
	for (int k=0; k<3; k++)
	  vel[k] -= 2.0*rv * pos[k]/rr;

				// Undo the center shift and reassign
	for (int k=0; k<3; k++) {
	  if (cC->com_system) vel[k] += cC->cov0[k];
	  p->vel[k] = vel[k];
	}

      } else if (debug) wrong_dir[id]++;

				// Reflect the position from the edge
				//
      delr = rr - radius;
      if (fabs(delr) > 0.2*radius) {
				// Make it sane, probably the result of CBA,
				// by choosing a random radius
	delr = radius * unit(gen[id]);
	for (int k=0; k<3; k++)
	  pos[k] = delr * pos[k]/rr;

      } else {

	for (int k=0; k<3; k++)
	  pos[k] -= 2.0*delr * pos[k]/rr;

      }

				// Undo the center shifts and reassign
      for (int k=0; k<3; k++) {
	pos[k] += cC->center[k];
	if (cC->com_system) pos[k] += cC->com0[k];
	p->pos[k] = pos[k];
      }

    }
    
  }
  
  return (NULL);
}


extern "C" {
  ExternalForce *makerReflect(const YAML::Node& conf)
  {
    return new UserReflect(conf);
  }
}

class proxyrefl { 
public:
  proxyrefl()
  {
    factory["userreflect"] = makerReflect;
  }
};

static proxyrefl p;
