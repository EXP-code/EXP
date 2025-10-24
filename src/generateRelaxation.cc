
#include "expand.H"

#include "generateRelaxation.H"

generateRelaxation::generateRelaxation(const YAML::Node& conf) : ExternalForce(conf)
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

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {

    i = (it++)->first;

    esave = 0.0;
    for (int j=0; j<3; j++)
      esave += cC->Vel(i, j) * cC->Vel(i, j);
  
    esave *= 0.5*cC->Mass(i);
    esave += cC->Mass(i) *
      ( cC->Part(i)->pot + cC->Part(i)->potext );
    cC->Part(i)->dattrib.push_back(esave);
    
  }
  
  done = 1;

  return (NULL);
}


