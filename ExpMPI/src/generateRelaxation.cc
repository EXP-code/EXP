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

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  if (id==0) epos = cC->Part(0)->dattrib.size();

  double esave;

  for (int i=nbeg; i<nend; i++) {

    esave = 0.0;
    for (int j=0; j<3; j++)
      esave += cC->Vel(i, j) * cC->Vel(i, j);
  
    esave *= 0.5*cC->Mass(i);
    esave += cC->Mass(i) *
      ( cC->Part(i)->pot + cC->Part(i)->potext );
    cC->Part(i)->dattrib.push_back(esave);
    
  }
  
  done = 1;
}


