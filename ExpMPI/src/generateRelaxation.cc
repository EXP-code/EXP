#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include "expand.h"

#include <generateRelaxation.H>

generateRelaxation::generateRelaxation(string& line) : ExternalForce(line)
{
  done = 0;
}

void generateRelaxation::initialize()
{

}

void * generateRelaxation::determine_acceleration_and_potential_thread
(void *arg)
{
  if (done) return (NULL);

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  if (id==0) epos = (*particles)[0].dattrib.size();

  double esave;

  for (int i=nbeg; i<nend; i++) {

    esave = 0.0;
    for (int j=0; j<3; j++)
      esave += (*particles)[i].vel[j] * (*particles)[i].vel[j];
  
    esave *= 0.5*(*particles)[i].mass;
    esave += (*particles)[i].mass*
      ( (*particles)[i].pot + (*particles)[i].potext );
    (*particles)[i].dattrib.push_back(esave);
    
  }
  
  done = 1;
}


