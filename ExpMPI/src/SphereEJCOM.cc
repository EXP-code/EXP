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
  alpha = 1.0;

  				// Get initialization info
  initialize();
}

void SphereEJCOM::initialize()
{
  string val;

  if (get_value("cfac", val)) cfac = atof(val.c_str());
  if (get_value("alpha", val)) alpha = atof(val.c_str());
}


double SphereEJCOM::mixture(double* pos)
{
  double dej=0.0, dif=0.0;

  for (int k=0; k<3; k++) {
    dej += 
      (pos[k] - component->center[k]) *
      (pos[k] - component->center[k]) ;

    if (component->com_system)
      dif += 
	(component->comI[k] - component->center[k]) *
	(component->comI[k] - component->center[k]) ;
    else
      dif += 
	(component->com[k] - component->center[k]) *
	(component->com[k] - component->center[k]) ;
  }

  double value = erf(cfac*pow(dej/(dif+1.0e-10), 0.5*alpha));

  accum_histo(value);

  return value;
}
