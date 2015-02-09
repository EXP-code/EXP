#include <mpi.h>
#include <SatFix.H>

SatFix::SatFix(string &line) : ExternalForce(line)
{
  verbose = true;
  debug   = false;
  comp_name = "Points";		// Default component for fixing

  initialize();

				// Look for the fiducial component
  bool found = false;
  list<Component*>::iterator cc;
  Component *c;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    if ( !comp_name.compare(c->name) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  // Find out who has particles, make sure that there are an even number
  if (2*(c0->nbodies_tot/2) != c0->nbodies_tot) {
    if (myid==0) cerr << "SatFix: component <" << comp_name 
		      << "> has an odd number of particles!!! nbodies_tot=" 
		      << c0->nbodies_tot << "\n";
    MPI_Abort(MPI_COMM_WORLD, 36);
  }

  if (!found) {
    cerr << "Process " << myid << ": SatFix can't find desired component <"
	 << comp_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }

  owner = vector<int>(c0->nbodies_tot);

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

void SatFix::initialize()
{
  string val;

  if (get_value("compname", val))   comp_name = val;
  if (get_value("verbose",  val))   if (atoi(val.c_str())) verbose = true;
  if (get_value("debug",    val))   if (atoi(val.c_str())) debug = true;
}

void SatFix::get_acceleration_and_potential(Component* C)
{
  if (C != c0) return;

  compute_list();

  MPL_start_timer();

  //
  // Begin list
  //
  for (unsigned n=0; n<c0->nbodies_tot-1; n+=2) {

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
				// Get body list
  for (unsigned n=0; n<c0->nbodies_tot; n++) {
    if (c0->Particles().find(n+1) != c0->Particles().end())
      owner[n] = myid+1;
    else
      owner[n] = 0;
  }

  MPI_Allreduce(MPI_IN_PLACE, &owner[0], c0->nbodies_tot, MPI_UNSIGNED,
		MPI_SUM, MPI_COMM_WORLD);
  
  for (unsigned n=0; n<c0->nbodies_tot; n++) {
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

