static char rcsid[] = "$Id$";

#include "expand.h"

#include <Direct.H>

Direct::Direct(string& line) : PotAccel(line)
{
  soft_indx = 0;

  initialize();

  if (component->ndattrib<soft_indx+1) {
    if (myid==0) cerr << "Direct: particle softening data missing\n";
    MPI_Abort(MPI_COMM_WORLD, 103);
    exit(0);
  }
    
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
  exp_thread_fork(false);

  // Clear external potential flag
  use_external = false;
}

void * Direct::determine_acceleration_and_potential_thread(void * arg)
{
  double rr, rfac;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  use[id] = 0;

  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;


    for (int j=0; j<i; j++) {

      rr = (*particles)[i].dattrib[soft_indx];
      rr = rr*rr;
      for (int k=0; k<3; k++)
	rr += (*particles)[i].pos[k] - (*particles)[j].pos[k];
      rr = sqrt(rr);

      rfac = 1.0/(rr*rr*rr);

      for (int k=0; k<3; k++)
	(*particles)[i].acc[k] += -(*particles)[j].mass *
	  ((*particles)[i].pos[k] - (*particles)[j].pos[k]) * rfac;

      if (use_external)
	(*particles)[i].potext += -(*particles)[j].mass/rr;
      else
	(*particles)[i].pot += -(*particles)[j].mass/rr;
    }

  }
  
}

void Direct::determine_coefficients(void) {}
void * Direct::determine_coefficients_thread(void *arg) {}

