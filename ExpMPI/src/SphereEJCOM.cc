#include <values.h>

#include "expand.h"

#include <SphereEJCOM.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

SphereEJCOM::SphereEJCOM(string&line) : SphereTwoCenter(line)
{
  id = "SphereEJCOM";
				// Defaults
  cfac = 1.0;

  				// Get initialization info
  initialize();
}

void SphereEJCOM::initialize()
{
  string val;

  if (get_value("cfac", val)) cfac = atof(val.c_str());
}


double SphereEJCOM::mixture(Particle& p)
{
  double dej=0.0, dif=0.0;

  for (int k=0; k<3; k++) {
    dej += 
      (p.pos[k] - component->center[k]) *
      (p.pos[k] - component->center[k]) ;

    dif += 
      (component->com[k] - component->center[k]) *
      (component->com[k] - component->center[k]) ;
  }

  return erf(cfac*sqrt(dej/(dif+1.0e-10)));
}
