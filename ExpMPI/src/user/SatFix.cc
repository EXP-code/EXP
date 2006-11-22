#include <mpi.h>
#include <SatFix.H>

SatFix::SatFix(string &line) : ExternalForce(line)
{
  verbose = false;
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
  if (2*(c0->nbodies_tot)/2 != c0->nbodies_tot) {
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

  plocate = vector<int>(numprocs, 0);
  vector<int> plocate1 = vector<int>(numprocs, 0);
  plocate1[myid] = c0->Number();

  MPI_Allreduce(&plocate1[0], &plocate[0], numprocs, 
		MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  userinfo();
}

SatFix::~SatFix()
{
}

void SatFix::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();
  cout << "Enforces mirror coordinates for adjacent particles\n";
  print_divider();
}

void SatFix::initialize()
{
  string val;

  if (get_value("compname", val))   comp_name = val;
  if (get_value("verbose", val))    if (atoi(val.c_str())) verbose = true;
}

void SatFix::get_acceleration_and_potential(Component* C)
{
  cC = C;			// "Register" component
  nbodies = cC->Number();	// And compute number of bodies

  if (verbose) cout << "Process " << myid << ": nbodies=" << nbodies << endl;

  if (!nbodies) return;		// Do we have any particle?

  MPI_Status status;
  
  MPL_start_timer();

  int recv = 0;
  int send = 0;
		       // Check for info from previous process
  if (myid) {
    
    MPI_Recv(&recv, 1, MPI_INT, myid-1, 132, MPI_COMM_WORLD, &status);
    
    if (verbose) cout << "Process " << myid << ": recv=" << recv << " from Process " << myid-1 << endl;

    if (recv) {
      MPI_Recv(C->Part(0)->pos, 3, MPI_DOUBLE, myid-1, 133, MPI_COMM_WORLD, &status);
      MPI_Recv(C->Part(0)->vel, 3, MPI_DOUBLE, myid-1, 134, MPI_COMM_WORLD, &status);
      MPI_Recv(C->Part(0)->acc, 3, MPI_DOUBLE, myid-1, 135, MPI_COMM_WORLD, &status);
      
				// Change sign of phase space
      for (int k=0; k<3; k++) {
	C->Part(0)->pos[k] *= -1.0;
	C->Part(0)->vel[k] *= -1.0;
	C->Part(0)->acc[k] *= -1.0;
      }
    }

  }

  if (!(nbodies-recv)) return;	// No more bodies?  
  int nbodiesE = 2*((nbodies-recv)/2);
  int lastn = 0;

  if (verbose) cout << "Process " << myid << ": nbodiesE=" << nbodiesE << endl;

  for (int n=recv; n<nbodiesE; n+=2) {

    for (int k=0; k<3; k++) {
      C->Part(n+1)->pos[k] = C->Part(n)->pos[k];
      C->Part(n+1)->vel[k] = C->Part(n)->vel[k];
      C->Part(n+1)->acc[k] = C->Part(n)->acc[k];
    }
    if (verbose) lastn = n+1;
  } 

  if (verbose) cout << "Process " << myid << ": last particle used=" << lastn << endl;

				// Send particle to mirror?
  if (myid+1 < numprocs) {

    if (nbodies-recv != nbodiesE) send = 1;

    MPI_Send(&send, 1, MPI_INT, myid+1, 132, MPI_COMM_WORLD);

    if (send) {

      if (verbose) cout << "Process " << myid 
			<< ": sending particle " << nbodies -1 << " to Process " << myid+1 << endl;

      MPI_Send(C->Part(nbodies-1)->pos, 3, MPI_DOUBLE, myid+1, 133, MPI_COMM_WORLD);
      MPI_Send(C->Part(nbodies-1)->vel, 3, MPI_DOUBLE, myid+1, 134, MPI_COMM_WORLD);
      MPI_Send(C->Part(nbodies-1)->acc, 3, MPI_DOUBLE, myid+1, 135, MPI_COMM_WORLD);
    }
  }

  MPL_stop_timer();
}

