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

  last = vector<int>(numprocs, 0);

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
  if (C != c0) return;

				// Get body count
  vector<int> ncount = c0->particle_count();
				// Check for change in body count
  bool recompute = false;
  for (int n=0; n<numprocs; n++) {
    if (last[n] != ncount[n]) recompute = true;
  }
  
  if (verbose) 
    cout << "Process " << myid << ": nbodies=" << ncount[myid] << endl;

  if (!ncount[myid]) return;	// Do we have any particles?

  MPI_Status status;
  
  if (recompute) {
				// Deal with end points
    send = vector<int>(numprocs, -1);
    recv = vector<int>(numprocs, -1);
    bool remainder = false;
    int from, number;

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
				// Debug
    if (verbose && myid==0) {

      cout << "Node list" << endl 
	   << "---------" << endl
	   << setw(6) << "Node" 
	   << setw(6) << "Send" 
	   << setw(6) << "Recv" << endl;
      
      for (int n=0; n<numprocs; n++) {
	cout << setw(6) << n;
      if (send[n]<0) cout << setw(6) << "*";
      else cout << setw(6) << send[n];
      if (recv[n]<0) cout << setw(6) << "*";
      else cout << setw(6) << recv[n];
      cout << endl;
      }
    }

    begin = 0;
    end   = ncount[myid];

    if (recv[myid]>=0) begin++;
    if (send[myid]>=0) end--;

    last = ncount;
    recompute = false;

    if (verbose) {
				// Should always have an even number of 
				// particles on each node after end points 
				// are removed
      int total = end - begin;
      if (total>0) assert( 2*(total/2) == total );
    }
    
  }

  MPL_start_timer();

				// Check for info from previous process
  if (recv[myid]>=0) {
    
    if (verbose) cout << "Process " << myid << ": receiving from Process " 
		      << recv[myid] << endl;

    MPI_Recv(C->Part(0)->pos, 3, MPI_DOUBLE, recv[myid], 133, 
	     MPI_COMM_WORLD, &status);
    MPI_Recv(C->Part(0)->vel, 3, MPI_DOUBLE, recv[myid], 134, 
	     MPI_COMM_WORLD, &status);
    MPI_Recv(C->Part(0)->acc, 3, MPI_DOUBLE, recv[myid], 135, 
	     MPI_COMM_WORLD, &status);
    
				// Change sign of phase space
    for (int k=0; k<3; k++) {
      C->Part(0)->pos[k] *= -1.0;
      C->Part(0)->vel[k] *= -1.0;
      C->Part(0)->acc[k] *= -1.0;
    }

  }

				// Send info to next process
  if (send[myid]>=0) {
    
    if (verbose) cout << "Process " << myid << ": sending to Process " 
		      << send[myid] << endl;

    MPI_Recv(C->Part(end)->pos, 3, MPI_DOUBLE, send[myid], 133, 
	     MPI_COMM_WORLD, &status);
    MPI_Recv(C->Part(end)->vel, 3, MPI_DOUBLE, send[myid], 134, 
	     MPI_COMM_WORLD, &status);
    MPI_Recv(C->Part(end)->acc, 3, MPI_DOUBLE, send[myid], 135, 
	     MPI_COMM_WORLD, &status);

  }

  for (int n=begin; n<end; n+=2) {
    for (int k=0; k<3; k++) {
      C->Part(n+1)->pos[k] = -C->Part(n)->pos[k];
      C->Part(n+1)->vel[k] = -C->Part(n)->vel[k];
      C->Part(n+1)->acc[k] = -C->Part(n)->acc[k];
    }
  } 

  MPL_stop_timer();
}

