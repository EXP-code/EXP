#include <mpi.h>
#include <SatFix.H>
#include <YamlCheck.H>
#include <EXPException.H>


SatFix::SatFix(const YAML::Node& conf) : ExternalForce(conf)
{
  verbose = true;
  debug   = false;
  comp_name = "Points";		// Default component for fixing

  initialize();

				// Look for the fiducial component
  bool found = false;
  for (auto cc : comp->components) {
    if ( !comp_name.compare(cc->name) ) {
      c0 = cc;
      found = true;
      break;
    }
  }

  unsigned total = c0->CurTotal();

  // Find out who has particles, make sure that there are an even number
  if (2*(total/2) != total) {
    std::ostringstream sout;
    sout << "SatFix: component <" << comp_name 
	 << "> has an odd number of particles!!! nbodies_tot=" 
	 << total;
    throw GenericError(sout.str(), __FILE__, __LINE__, 36, false);
  }

  if (!found) {
    std::ostringstream sout;
    sout << "Process " << myid << ": SatFix can't find desired component <"
	 << comp_name << ">";
    throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
  }

  owner = std::vector<int>(total);

  userinfo();
}

SatFix::~SatFix()
{
}

void SatFix::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();
  cout << "** Enforces mirror coordinates for adjacent particles on component: " 
       << c0->name << endl;
  print_divider();
}

const std::set<std::string>
SatFix::valid_keys = {
  "compname",
  "verbose",
  "debug"
};

void SatFix::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("SatFix", "parameter", unmatched, __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["compname"])       comp_name          = conf["compname"].as<string>();
    if (conf["verbose"])        verbose            = conf["verbose"].as<bool>();
    if (conf["debug"])          debug              = conf["debug"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in SatFix: "
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

void SatFix::get_acceleration_and_potential(Component* C)
{
  if (C != c0) return;

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(c0);
#endif

  compute_list();

  MPL_start_timer();

  //
  // Begin list
  //
  unsigned total = c0->CurTotal();

  for (unsigned n=0; n<total-1; n+=2) {

				// Who has the odd particle?
    if (myid==owner[n]) {
				// Temporary!  Remove this . . .
#if 0
      print_body(n+1, "send");
#endif

      check_body(n+1);

      MPI_Send(C->Part(n+1)->pos, 3, MPI_DOUBLE, owner[n+1], 133, 
	       MPI_COMM_WORLD);
      MPI_Send(C->Part(n+1)->vel, 3, MPI_DOUBLE, owner[n+1], 134, 
	       MPI_COMM_WORLD);
      MPI_Send(C->Part(n+1)->acc, 3, MPI_DOUBLE, owner[n+1], 135, 
	       MPI_COMM_WORLD);
    }

				// Who has the even particle?
    if (myid==owner[n+1]) {
      MPI_Recv(C->Part(n+2)->pos, 3, MPI_DOUBLE, owner[n], 133, 
	       MPI_COMM_WORLD, &status);
      MPI_Recv(C->Part(n+2)->vel, 3, MPI_DOUBLE, owner[n], 134, 
	       MPI_COMM_WORLD, &status);
      MPI_Recv(C->Part(n+2)->acc, 3, MPI_DOUBLE, owner[n], 135, 
	       MPI_COMM_WORLD, &status);

				// Temporary!  Remove this . . .
#if 0
      print_body(n+2, "recv");
#endif

				// Change sign of phase space
      for (int k=0; k<3; k++) {
	C->Part(n+2)->pos[k] *= -1.0;C->Part(n+2)->vel[k] *= -1.0;
	C->Part(n+2)->acc[k] *= -1.0;
      }
      
      check_body(n+2);
    }

  }

  MPL_stop_timer();
}

void SatFix::compute_list()
{
  unsigned total = c0->CurTotal();

				// Get body list
  for (unsigned n=0; n<total; n++) {
    if (c0->Particles().find(n+1) != c0->Particles().end())
      owner[n] = myid+1;
    else
      owner[n] = 0;
  }

  MPI_Allreduce(MPI_IN_PLACE, &owner[0], total, MPI_UNSIGNED,
		MPI_SUM, MPI_COMM_WORLD);
  
  for (unsigned n=0; n<total; n++) {
    owner[n] -= 1;
    if (owner[n]<0) {
      cout << "SatFix: error in ownership list" << endl;
    }
  }

}

void SatFix::check_body(int n)
{
  if (debug) {
    bool ferror = false;
    for (int k=0; k<3; k++) {
      if (std::isnan(c0->Part(n)->pos[k])) ferror = true;
      if (std::isnan(c0->Part(n)->vel[k])) ferror = true;
      if (std::isnan(c0->Part(n)->acc[k])) ferror = true;
    }
    if (ferror) {
      cout << "Process " << myid << ": error in coordindates, n=" << n << "!" << endl;
    }
  }
}


void SatFix::print_body(int n, string& mesg)
{
  if (debug) {
    cout << endl << "Process " << myid << ": " << mesg << " pos=";
    for (int k=0; k<3; k++) cout << setw(16) << c0->Part(n)->pos[k];
    cout << endl << "Process " << myid << ": " << mesg << " vel=";
    for (int k=0; k<3; k++) cout << setw(16) << c0->Part(n)->vel[k];
    cout << endl << "Process " << myid << ": " << mesg << " acc=";
    for (int k=0; k<3; k++) cout << setw(16) << c0->Part(n)->acc[k];
  }
}

