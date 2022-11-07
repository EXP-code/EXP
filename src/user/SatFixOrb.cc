#include <mpi.h>
#include <cassert>
#include <SatFixOrb.H>
#include <YamlCheck.H>
#include <EXPException.H>

const std::set<std::string>
SatFixOrb::valid_keys = {
  "compname",
  "config",
  "toffset",
  "verbose",
  "debug"
};

SatFixOrb::SatFixOrb(const YAML::Node& conf) : ExternalForce(conf)
{
  verbose = true;
  debug   = false;
  toffset = 0.0;
  tag     = 0;			// Tag to match for enforcing orbit
  comp_nam = "Points";		// Default component for fixing
  config   = "conf.file";	// Configuration file for spherical orbit

  initialize();

				// Look for the fiducial component
  bool found = false;
  for (auto c : comp->components) {
    if ( !comp_nam.compare(c->name) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  if (!found) {
    std::ostringstream sout;
    sout << "Process " << myid << ": SatFixOrb can't find desired component <"
	 << comp_nam << ">";
    throw GenericError(sout.str(), __FILE__, __LINE__ 35, false);
  }

  unsigned total = c0->CurTotal();

  // Find out who has particles, make sure that there are an even number
  //
  if (2*(total/2) != total) {
    std::ostringstream sout;
    sout << "SatFixOrb: component <" << comp_nam 
	 << "> has an odd number of particles!!! nbodies_tot=" 
	 << total;
    throw GenericError(sout.str(), __FILE__, __LINE__ 36, false);
  }

  // Make sure that the particles are tagged
  //
  int in=0, out;
  if (c0->Number()) in = c0->Part(0)->iattrib.size();
  MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (out<1) {
    std::string msg = "SatFixOrb: bodies need an integer tag!";
    throw GenericError(sout.str(), __FILE__, __LINE__ 37, false);
  }

  last = vector<unsigned int>(numprocs, 0);

  orb = std::make_shared<SatelliteOrbit>(conf);

  userinfo();
}

SatFixOrb::~SatFixOrb()
{
  delete orb;
}

void SatFixOrb::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** Enforces mirror coordinates for adjacent particles on component: " 
       << c0->name << endl 
       << "and enforces an analytic orbit from config file <" 
       << config << "> with Toffset=" << toffset << endl;

  print_divider();
}

void SatFixOrb::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("SatFixOrb", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["compname"])       comp_nam           = conf["compname"].as<string>();
    if (conf["config"])         config             = conf["config"].as<string>();
    if (conf["toffset"])        toffset            = conf["toffset"].as<double>();
    if (conf["verbose"])        verbose            = conf["verbose"].as<bool>();
    if (conf["debug"])          debug              = conf["debug"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in SatFixOrb: "
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

void SatFixOrb::get_acceleration_and_potential(Component* C)
{
  if (C != c0) return;

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(c0);
#endif

  compute_list();

  MPL_start_timer();

  enforce();
				// Temporary!  Remove this . . .
				// Check for info from previous process
  if (recv[myid]>=0) {
    
    MPI_Recv(C->Part(0)->pos, 3, MPI_DOUBLE, recv[myid], 133, 
	     MPI_COMM_WORLD, &status);
    MPI_Recv(C->Part(0)->vel, 3, MPI_DOUBLE, recv[myid], 134, 
	     MPI_COMM_WORLD, &status);
    MPI_Recv(C->Part(0)->acc, 3, MPI_DOUBLE, recv[myid], 135, 
	     MPI_COMM_WORLD, &status);
    
				// Temporary!  Remove this . . .
#if 0
    print_recv();
#endif
				// Change sign of phase space
    for (int k=0; k<3; k++) {
      C->Part(0)->pos[k] *= -1.0;
      C->Part(0)->vel[k] *= -1.0;
      C->Part(0)->acc[k] *= -1.0;
    }

    check_recv();
  }

				// Send info to next process
  if (send[myid]>=0) {
    
				// Temporary!  Remove this . . .
#if 0
    print_send();
#endif

    check_send();

    MPI_Send(C->Part(end)->pos, 3, MPI_DOUBLE, send[myid], 133, 
	     MPI_COMM_WORLD);
    MPI_Send(C->Part(end)->vel, 3, MPI_DOUBLE, send[myid], 134, 
	     MPI_COMM_WORLD);
    MPI_Send(C->Part(end)->acc, 3, MPI_DOUBLE, send[myid], 135, 
	     MPI_COMM_WORLD);
  }

  for (int n=begin; n<end; n+=2) {

    check_body(n);

    for (int k=0; k<3; k++) {
      C->Part(n+1)->pos[k] = -C->Part(n)->pos[k];
      C->Part(n+1)->vel[k] = -C->Part(n)->vel[k];
      C->Part(n+1)->acc[k] = -C->Part(n)->acc[k];
    }
  } 

  MPL_stop_timer();
}

void SatFixOrb::enforce()
{
  double rs[3];

  for (int n=0; n<c0->particle_count(myid); n++) {
    if (tag == c0->Part(n)->iattrib[0]) {
      orb->get_satellite_orbit(tnow - toffset, &rs[0]);
      for (int k=0; k<3; k++) c0->Part(n)->pos[k] = rs[k];
#if 1				// For debugging orbit fixing
      ostringstream sout;
      sout << outdir << "SatFixOrb.test." << myid;
      ofstream out(sout.str().c_str(), ios::app);
      out << setw(15) << tnow;
      for (int k=0; k<3; k++) out << setw(15) << c0->Part(n)->pos[k];
      out << endl;
#endif      
      break;
    }
  }

}

void SatFixOrb::compute_list()
{
				// Get body count
  ncount = c0->particle_count();

				// Check for change in body count
  bool recompute = false;
  for (int n=0; n<numprocs; n++) {
    if (last[n] != ncount[n]) recompute = true;
  }
  
  if (recompute) {
				// Deal with end points
    send = vector<int>(numprocs, -1);
    recv = vector<int>(numprocs, -1);
    bool remainder = false;
    int from=0, number;

    for (int n=0; n<numprocs; n++) {
      number = ncount[n];

      if (number>0) {
				// Particle to be sent from previous node
	if (remainder) {
	  send[from] = n;
	  recv[n] = from;
	  number--;
	  remainder = false;
	}
				// Particlde to be sent to next node
	if ( 2*(number/2) != number ) {
	  from = n;
	  remainder = true;
	}
	
      }
    }

    begin = 0;
    end   = ncount[myid];

    if (recv[myid]>=0) begin++;
    if (send[myid]>=0) end--;

    last = ncount;
    recompute = false;

				// Report particle allocation calculation
    if (verbose) {

      if (myid==0) {

	int begin1 = begin;
	int end1   = end;

	cout << "====================================================" << endl
	     << "Component: " << c0->name << endl
	     << "Node list for mirroring" << endl 
	     << "-----------------------" << endl
	     << setw(6) << "Node" 
	     << setw(6) << "Count"
	     << setw(6) << "Send" 
	     << setw(6) << "Recv" 
	     << setw(6) << "Begin" 
	     << setw(6) << "End" 
	     << endl;
	
	for (int n=0; n<numprocs; n++) {

	  if (n) {
	    MPI_Recv(&begin1, 1, MPI_INT, n, 131, MPI_COMM_WORLD, &status);
	    MPI_Recv(&end1,   1, MPI_INT, n, 132, MPI_COMM_WORLD, &status);
	  }

	  cout << setw(6) << n << setw(6) << ncount[n];
	  if (send[n]<0) cout << setw(6) << "*";
	  else cout << setw(6) << send[n];
	  if (recv[n]<0) cout << setw(6) << "*";
	  else cout << setw(6) << recv[n];
	  cout << setw(6) << begin1 << setw(6) << end1 << endl;
	}

	cout << "====================================================" << endl;
	
      } else {
	MPI_Send(&begin, 1, MPI_INT, 0, 131, MPI_COMM_WORLD);
	MPI_Send(&end,   1, MPI_INT, 0, 132, MPI_COMM_WORLD);
      }
				// Should always have an even number of 
				// particles on each node after end points 
				// are removed
      int total = end - begin;
      if (total>0) assert( 2*(total/2) == total );
    }
    
  }

}

void SatFixOrb::check_send()
{
  if (debug) {
    bool ferror = false;
    for (int k=0; k<3; k++) {
      if (std::isnan(c0->Part(end)->pos[k])) ferror = true;
      if (std::isnan(c0->Part(end)->vel[k])) ferror = true;
      if (std::isnan(c0->Part(end)->acc[k])) ferror = true;
    }
    if (ferror) {
      cout << "Process " << myid << ": error in coordindates to be sent!" << endl;
    }
  }
  
}

void SatFixOrb::check_body(int n)
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

void SatFixOrb::check_recv()
{
  if (debug) {
    bool ferror = false;
    for (int k=0; k<3; k++) {
      if (std::isnan(c0->Part(0)->pos[k])) ferror = true;
      if (std::isnan(c0->Part(0)->vel[k])) ferror = true;
      if (std::isnan(c0->Part(0)->acc[k])) ferror = true;
    }
    if (ferror) {
      cout << "Process " << myid << ": error in receiving coordindates!" << endl;
    }
  }
}


void SatFixOrb::print_recv()
{
  if (debug) {
    cout << endl << "Process " << myid << ": received pos=";
    for (int k=0; k<3; k++) cout << setw(16) << c0->Part(0)->pos[k];
    cout << endl << "Process " << myid << ": received vel=";
    for (int k=0; k<3; k++) cout << setw(16) << c0->Part(0)->vel[k];
    cout << endl << "Process " << myid << ": received acc=";
    for (int k=0; k<3; k++) cout << setw(16) << c0->Part(0)->acc[k];
  }
}

void SatFixOrb::print_send()
{
  if (debug) {
    cout << endl << "Process " << myid << ": sent pos=";
    for (int k=0; k<3; k++) cout << setw(16) << c0->Part(end)->pos[k];
    cout << endl << "Process " << myid << ": sent vel=";
    for (int k=0; k<3; k++) cout << setw(16) << c0->Part(end)->vel[k];
    cout << endl << "Process " << myid << ": sent acc=";
    for (int k=0; k<3; k++) cout << setw(16) << c0->Part(end)->acc[k];
  }
}

