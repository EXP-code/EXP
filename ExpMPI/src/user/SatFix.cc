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

  if (!found) {
    cerr << "Process " << myid << ": SatFix can't find desired component <"
	 << comp_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }

  // Find out who has particles, make sure that there are only two
  if (c0->nbodies_tot != 2) {
    if (myid==0) cerr << "SatFix: component <" << comp_name 
		      << "> has more than two particles!!! nbodies_tot=" 
		      << c0->nbodies_tot << "\n";
    MPI_Abort(MPI_COMM_WORLD, 36);
  }
  
  plocate = vector<int>(numprocs, 0);
  vector<int> plocate1 = vector<int>(numprocs, 0);
  plocate1[myid] = c0->particles.size(); 

  MPI_Allreduce(&plocate1[0], &plocate[0], numprocs, 
		MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (verbose && myid==0) 
    cout << "****************SatFix Debug****************\n";

  myid0 = myid1 = -1;
  for (int i=0; i<numprocs; i++) {
    if (verbose && myid==0) 
      cout << setw(5) << i << setw(10) << plocate[i] << endl;

    if (myid0<0 && myid1<0 && plocate[i] != 0) {
      myid0 = i;
      if (plocate[i]>1) myid1 = i;
    }

    if (myid0>=0 && myid1<0 && plocate[i] != 0) myid1 = i;
  }

  if (verbose && myid==0)  {
    cout << "Id 0=" << myid0 << ", Particle 0\n";
    if (myid0 == myid1)
      cout << "Id 1=" << myid1 << ", Particle 1\n";
    else
      cout << "Id 1=" << myid1 << ", Particle 0\n";
    
    cout << "********************************************\n";
  }    

  userinfo();
}

SatFix::~SatFix()
{
}

void SatFix::userinfo()
{
  if (myid) return;		// Return if node master node
  cout << "** Enforces mirror coordinates for adjacent particles         **\n";
  cout << "****************************************************************\n";
}

void SatFix::initialize()
{
  string val;

  if (get_value("compname", val))   comp_name = val;
  if (get_value("verbose", val))    if (atoi(val.c_str())) verbose = true;
}

void SatFix::get_acceleration_and_potential(vector<Particle>* P)
{
  particles = P;		// "Register" particles
  nbodies = (*particles).size(); // And compute number of bodies

  MPI_Status status;
  
  MPL_start_timer();

  if (myid0 == myid1 && myid0 == myid) {
    for (int k=0; k<3; k++) {
      (*P)[1].pos[k] = -(*P)[0].pos[k];
      (*P)[1].vel[k] = -(*P)[0].vel[k];
      (*P)[1].acc[k] = -(*P)[0].acc[k];
    }
  } else if (myid0 == myid) {	// Send particle to mirror
    MPI_Send(&((*P)[0].pos[0]), 3, MPI_DOUBLE, myid1, 133, MPI_COMM_WORLD);
    MPI_Send(&((*P)[0].vel[0]), 3, MPI_DOUBLE, myid1, 134, MPI_COMM_WORLD);
    MPI_Send(&((*P)[0].acc[0]), 3, MPI_DOUBLE, myid1, 135, MPI_COMM_WORLD);
  } else if (myid1 == myid) {	// Get particle from fiducial
    MPI_Recv(&((*P)[0].pos[0]), 3, MPI_DOUBLE, myid0, 133, MPI_COMM_WORLD, 
	     &status);
    MPI_Recv(&((*P)[0].vel[0]), 3, MPI_DOUBLE, myid0, 134, MPI_COMM_WORLD,
	     &status);
    MPI_Recv(&((*P)[0].acc[0]), 3, MPI_DOUBLE, myid0, 135, MPI_COMM_WORLD,
	     &status);
				// Change sign of phase space
    for (int k=0; k<3; k++) {
      (*P)[0].pos[k] *= -1.0;
      (*P)[0].vel[k] *= -1.0;
      (*P)[0].acc[k] *= -1.0;
    }
  }

  MPL_stop_timer();
}

