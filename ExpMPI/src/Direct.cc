static char rcsid[] = "$Id$";

#include "expand.h"
#include <localmpi.h>

#include <Direct.H>

static const int MSGTAG=103;
static const int NDIM=5;

Direct::Direct(string& line) : PotAccel(line)
{
  soft_indx = 0;

  initialize();

  if (component->ndattrib<soft_indx+1) {
    if (myid==0) cerr << "Direct: particle softening data missing\n";
    MPI_Abort(MPI_COMM_WORLD, 103);
    exit(0);
  }
    
				// Assign the ring connections
  to_proc = (myid+1) % numprocs;
  from_proc = (myid+numprocs-1) % numprocs;

				// Buffer pointers
  tmp_buffer = NULL;
  bod_buffer = NULL;

}

Direct::~Direct()
{
  delete [] tmp_buffer;
  delete [] bod_buffer;
}

void Direct::initialize(void)
{
  string val;

  if (get_value("soft_indx", val)) soft_indx = atoi(val.c_str());
}

void Direct::get_acceleration_and_potential(vector<Particle>* P)
{
  particles = P;
  nbodies = particles->size();

  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  determine_acceleration_and_potential();
}

void Direct::determine_acceleration_and_potential(void)
{
				// Determine size of largest nbody list
  ninteract = component->particles.size();
  max_bodies = ninteract;
  MPI_Reduce(&ninteract, &max_bodies, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

				// Allocate buffers to handle largest list
  delete [] tmp_buffer;
  delete [] bod_buffer;
  int buffer_size = max_bodies*NDIM;
  tmp_buffer = new double [buffer_size];
  bod_buffer = new double [buffer_size];

				// Load body buffer with local interactors
  double *p = bod_buffer;
  for (int i=0; i<ninteract; i++) {
    *(p++) = component->particles[i].mass;
    *(p++) = component->particles[i].pos[0];
    *(p++) = component->particles[i].pos[1];
    *(p++) = component->particles[i].pos[2];
    *(p++) = component->particles[i].dattrib[soft_indx];
  }

				// Do the local interactors
  local = true;
  exp_thread_fork(false);
  local = false;

				// Do the ring . . . 
  for(int n=1; n<numprocs; n++) {

    MPI_Request req1, req2;
    MPI_Status stat;

				// Copy current to temp buffer
    memcpy(tmp_buffer, bod_buffer, buffer_size*sizeof(double));

				// Get NEW buffer from right
    MPI_Irecv(bod_buffer, buffer_size, MPI_DOUBLE, from_proc, MSGTAG, 
	      MPI_COMM_WORLD, &req1);

				// Send OLD buffer to left
    MPI_Isend(tmp_buffer, ninteract*NDIM, MPI_DOUBLE, to_proc, MSGTAG, 
	      MPI_COMM_WORLD, &req2);

    MPI_Wait(&req2, &stat);
    MPI_Wait(&req1, &stat);

				// How many particles did we get?
    MPI_Get_count(&stat, MPI_DOUBLE, &ninteract);
    ninteract /= NDIM;
	
				// Accumulate the interactions
    exp_thread_fork(false);
  }

  // Clear external potential flag
  use_external = false;
}

void * Direct::determine_acceleration_and_potential_thread(void * arg)
{
  double rr, rfac;
  double mass, pos[3], eps;
  double *p;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  use[id] = 0;

  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;
    
    use[id]++;

    p = bod_buffer;
    for (int j=0; j<ninteract; j++) {

				// Get current interaction particle
      mass = *(p++);
      for (int k=0; k<3; k++) pos[k] = *(p++);
      eps = *(p++);

				// Compute distance
      rr = 0.0;
      for (int k=0; k<3; k++)
	rr += 
	  ((*particles)[i].pos[k] - pos[k]) *
	  ((*particles)[i].pos[k] - pos[k]) ;

				// Check for coincident particles (e.g. same)
      if (local && !use_external && rr<1.0e-10) continue;

				// Add softening
      rr += eps*eps;
      rr = sqrt(rr);

				// Acceleration
      rfac = 1.0/(rr*rr*rr);

      for (int k=0; k<3; k++)
	(*particles)[i].acc[k] += -mass *
	  ((*particles)[i].pos[k] - pos[k]) * rfac;

				// Potential
      if (use_external)
	(*particles)[i].potext += -mass/rr;
      else
	(*particles)[i].pot += -mass/rr;
    }

  }
  
}

void Direct::determine_coefficients(void) {}
void * Direct::determine_coefficients_thread(void *arg) {}

