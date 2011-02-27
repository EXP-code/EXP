#include <values.h>

#include "expand.h"

#include <EJcom.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

EJcom::EJcom(string&line) : TwoCenter(line)
{
  id = "EJcom";

  // Defaults
  //
  cfac  = 1.0;
  alpha = 1.0;

  // Get initialization info
  //
  initialize();
}

void EJcom::initialize()
{
  string val;

  if (get_value("cfac",  val)) cfac  = atof(val.c_str());
  if (get_value("alpha", val)) alpha = atof(val.c_str());
}


double EJcom::mixture(double* pos)
{
  double del=0.0, dif=0.0;

  for (int k=0; k<3; k++) {
    del += (pos[k]   - inner[k]) * (pos[k]   - inner[k]);
    dif += (outer[k] - inner[k]) * (outer[k] - inner[k]);
  }

  double value = erf(cfac*pow(del/(dif+1.0e-10), 0.5*alpha));
  
  if (multistep==0 || mstep==0) accum_histo(value);

  return value;
}
