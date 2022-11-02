#include <sys/timeb.h>
#include <math.h>
#include <sstream>
#include <random>

#include <expand.H>
#include <global.H>

#include <UserSNheat.H>

//
// Physical units
//

static double pc = 3.086e18;		// cm
// static double a0 = 2.0*0.054e-7;	// cm (2xBohr radius)
// static double boltz = 1.381e-16;	// cgs
static double year = 365.25*24*3600;	// seconds per year
// static double mp = 1.67e-24;		// g
static double msun = 1.989e33;		// g

double UserSNheat::Lunit = 3.0e5*pc;
double UserSNheat::Munit = 1.0e12*msun;
double UserSNheat::Tunit = sqrt(Lunit*Lunit*Lunit/(Munit*6.673e-08));
double UserSNheat::Vunit = Lunit/Tunit;
double UserSNheat::Eunit = Munit*Vunit*Vunit;

const std::set<std::string>
UserSNheat::valid_keys = {
  "compname",    
  "X",
  "Y",
  "Z",
  "dT",
  "dE",
  "radius",
  "delay",
  "number",
  "Lunit",
  "Tunit",
  "Munit"
};

UserSNheat::UserSNheat(const YAML::Node& conf) : ExternalForce(conf)
{

  id = "SupernovaHeating";	// ID

				// Location heat source
  origin = vector<double>(3, 0.0);	
				// Radius of spherical region for heating
				// 0.0001 is 30 pc
  radius = 0.0001;
				// Delay in system units
  delay = 0.0;
				// Spacing in years
  dT = 1.0e+04;
				// Energy per SN in erg
  dE = 1.0e+51;
				// Number of SN
  N = 100;
				// Default component (must be specified)
  comp_name = "";
				// Debugging output switch
  verbose = false;

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
    cerr << "UserSNheat: process " << myid 
	 << " can't find fiducial component <" << comp_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }
  
  userinfo();

  Vunit = Lunit/Tunit;
  Eunit = Munit*Vunit*Vunit;

  dT *= year/Tunit;
  dE /= Eunit;

  if (myid==0 && verbose) {
    cout << "UserSNheat: dT=" << dT << " dE=" << dE << endl;
  }
  //

  firstime = true;
  ncount   = 0;
  mm0      = vector<double>(nthrds);
  ke0      = vector<double>(nthrds);
  ke1      = vector<double>(nthrds);
  mom      = vector<double>(3);
  pp0      = vector< vector<double> >(nthrds);

  for (int k=0; k<nthrds; k++) pp0[k] = vector<double>(3, 0.0);

  threading_init();

  // Recover from old log file
  //
  if (myid==0) {

    // Name of original and backup log files
    //
    filename = outdir + "SNheat." + runtag;
    string backupfile = filename + ".bak";

    // Attempt to open old log file, does it exist?
    // If so, backup by renaming
    //
    ifstream chk(filename.c_str());
    if (chk) {
      chk.close();

      // Attempt to rename
      //
      if (rename(filename.c_str(), backupfile.c_str())) {
	perror("UserSNheat::Run()");
	std::ostringstream message;
	message << "UserSNheat: could not rename log file <" 
		<< filename << "> to <" << backupfile << ">";
	throw GenericError(message.str(), __FILE__, __LINE__);
      }
      
      // Open (original) backup file for reading
      //
      ifstream in(backupfile.c_str());
      if (!in) {
	cerr << "UserSNheat: error opening original log file <" 
	     << backupfile << "> for reading" << endl
	     << "UserSNheat: assuming no SN so far" << endl;
      } else {
	
	// Open new output stream for writing
	//
	ofstream out(filename.c_str());
	if (!out) {
	  ostringstream message;
	  message << "UserSNheat: error opening new log file <" 
		  << filename << "> for writing";
	  throw GenericError(message.str(), __FILE__, __LINE__);
	}
	
	const size_t cbufsiz = 4096;
	char cbuffer[cbufsiz];
	double T, Tphys;
	int n;
      
	// Now read the backup file and write the new log file
	//
	while (in) {
	  in.getline(cbuffer, cbufsiz);
	  if (!in) break;
	  
	  istringstream sin(cbuffer);
	  
	  sin >> T;
	  if (T>tnow) break;
	  
	  sin >> Tphys;
	  sin >> n;
	  sin >> ncount;
	  
	  out.write((const char *)cbuffer, cbufsiz);
	}
      }
    }
  }

  // Share the current SN count (either initial for from restart)
  //
  MPI_Bcast(&ncount, 1, MPI_INT, 0, MPI_COMM_WORLD);

}

UserSNheat::~UserSNheat()
{
}


int UserSNheat::arrivalTime(double dt)
{
  double L = exp(-dt/dT);
  double p = 1.0;
  int    k = 0;
  do {
    k++;
    p *= unit(random_gen);
  } while (p>L);

  return k - 1;
}

void UserSNheat::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine stochastic Supernova heating initialized"
       << " using component <" << comp_name << ">"
       << " with delay time=" << delay << ", time interval dT=" << dT 
       << ", SN energy dE=" << dE << ", number SN=" << N 
       << ", bubble radius=" << radius 
       << endl
       << "   Lunit=" << Lunit << ", Tunit=" << Tunit << ", Munit=" << Munit
       << endl
       << "   Origin (x , y , z) = (" 
       << origin[0] << " , " 
       << origin[1] << " , " 
       << origin[2] << " ) " << endl; 

  print_divider();
}

void UserSNheat::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserSNheat", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["compname"])       comp_name          = conf["compname"].as<string>();
    if (conf["verbose"])        verbose            = conf["verbose"].as<int>();
    
    if (conf["X"])              origin[0]          = conf["X"].as<double>();
    if (conf["Y"])              origin[1]          = conf["Y"].as<double>();
    if (conf["Z"])              origin[2]          = conf["Z"].as<double>();
    
    if (conf["dT"])             dT                 = conf["dT"].as<double>();
    if (conf["dE"])             dE                 = conf["dE"].as<double>();
    if (conf["radius"])         radius             = conf["radius"].as<double>();
    if (conf["delay"])          delay              = conf["delay"].as<double>();
    if (conf["number"])         N                  = conf["number"].as<int>();
    
    if (conf["Lunit"])          Lunit              = conf["Lunit"].as<double>();
    if (conf["Tunit"])          Tunit              = conf["Tunit"].as<double>();
    if (conf["Munit"])          Munit              = conf["Munit"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserSNheat: "
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

void UserSNheat::determine_acceleration_and_potential(void)
{
  if (cC != c0)               return;
  if (multistep and mlevel>0) return;
  if (tnow < delay)           return;
  if (ncount > N)             return;
  
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  if (!firstime) {
    
    if (myid==0) {
      nSN = arrivalTime(tnow - tlast);
      if (nSN) {
	if (verbose)
	  cout << "UserSNheat: T=" << setw(12) << tnow 
	       << " [" << setw(12) << tnow*Tunit/year << " years]"
	       << "     SN=" << setw(4) << nSN 
	       << "     so far=" << setw(4) << ncount << endl;
	ofstream out(filename.c_str(), ios::app | ios::out);
	if (out) {
	  out << setw(18) << tnow
	      << setw(18) << tnow*Tunit/year
	      << setw(8)  << nSN
	      << setw(8)  << ncount
	      << endl;
	}
      }
    }
    MPI_Bcast(&nSN, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (nSN) {
      plist = vector< set<int> >(nthrds);
      exp_thread_fork(false);
      ncount += nSN;
    }
    print_timings("UserSNheat: thread timings");

  }

  tlast    = tnow;
  firstime = false;

}

void * UserSNheat::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  
  thread_timing_beg(id);

  PartMapItr it = cC->Particles().begin();

  mm0[id] = ke0[id] = ke1[id] = 0.0;
  for (int k=0; k<3; k++) pp0[id][k] = 0.0;
  
  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    
				// Index for the current particle
    unsigned long i = (it++)->first;
    
    Particle *p = cC->Part(i);
    
    double dist = 0.0;
    for (int k=0; k<3; k++) {
      double pos = p->pos[k] - origin[k];
      dist += pos*pos;
    }

    if (dist < radius*radius) {
      plist[id].insert(i);
      mm0[id] += p->mass;
      for (int k=0; k<3; k++) pp0[id][k] += p->mass*p->vel[k];
    }
  }
  
  Pbarrier(id, 37);

  if (id==0) {			// Thread 0 add up contributions to
				// mass and momentum from remaining
				// threads
    for (int i=1; i<nthrds; i++) {
      mm0[0] += mm0[i];
      for (int k=0; k<3; k++) pp0[0][k] += pp0[i][k];
    }

    MPI_Allreduce(&mm0[0],    &mass0,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pp0[0][0], &mom[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  Pbarrier(id, 38);

  if (mass0 > 0.0) {
				// Compute the center of mass velocity
    for (int k=0; k<3; k++) mom[k] /= mass0;

    for (auto s : plist[id]) {
      Particle *p = cC->Part(s);
      for (int k=0; k<3; k++) {
	double vel = p->vel[k] - mom[k];
	ke0[id] += 0.5*p->mass * vel*vel;
      }
    }

    Pbarrier(id, 39);

    if (id==0) {
      for (int i=1; i<nthrds; i++) ke0[0] += ke0[i];
      MPI_Allreduce(&ke0[0], &ketot0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    Pbarrier(id, 40);
    
    double disp = sqrt(2.0/3.0*(dE*nSN+ketot0)/mass0);
    for (auto s : plist[id]) {
      Particle *p = cC->Part(s);
      for (int k=0; k<3; k++) {
	double vel = disp*norm(random_gen);
	ke1[id] += 0.5*p->mass * vel*vel;
	p->vel[k] = mom[k] + vel;
      }
    }

    Pbarrier(id, 41);
    
    if (id == 0) {
      for (int j=1; j<nthrds; j++) ke1[0] += ke1[j];
      MPI_Allreduce(&ke1[0], &ketot1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      factor = sqrt((dE*nSN+ketot0)/ketot1);
    }

    Pbarrier(id, 42);
    
    for (auto s : plist[id]) {
      Particle *p = cC->Part(s);
      for (int k=0; k<3; k++) p->vel[k] *= factor;
    }
    
    if (id==0 && myid==0 && verbose)
      cout << "UserSNheat: " << "mass=" << mass0 
	   << ", factor=" << factor << ", snE=" << dE*nSN
	   << ", ke0=" << ketot0 << ", ke1=" << ketot1
	   << endl;

  } else {
    if (id==0 && myid==0)
      cout << "UserSNheat: "
	   << "No points in heating sphere of radius " << radius 
	   << " at time " << tnow << endl;
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerSNheat(const YAML::Node& conf)
  {
    return new UserSNheat(conf);
  }
}

class proxysnheat { 
public:
  proxysnheat()
  {
    factory["usersnheat"] = makerSNheat;
  }
};

static proxysnheat p;
