
#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include <values.h>

#include "expand.h"

#include <gaussQ.h>
#include <SphereTwoCenter.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

SphereTwoCenter::SphereTwoCenter(string& line) : PotAccel(line)
{
  id = "SphereTwoCenter SL";
				// Defaults
  rmin = 1.0e-3;
  rmax = 2.0;
  rs = 0.067*rmax;
  numr = 2000;
  cmap = 1;
  Lmax = 4;
  nmax = 10;

  firstime_accel = true;
  self_consistent = true;

				// Get initialization info
  initialize();

  SLGridSph::mpi = 1;		// Turn on MPI

				// Generate Sturm-Liouville grid
  ortho = new SLGridSph(Lmax, nmax, numr, rmin, rmax, cmap, rs);
  
				// Generate two expansion grids
  exp_ej  = new SphericalBasisMixtureSL(line, this, ej);
  exp_com = new SphericalBasisMixtureSL(line, this, com);
  center = new double [3];
}


void SphereTwoCenter::initialize()
{
  string val;

  if (get_value("rmin", val)) rmin = atof(val.c_str());
  if (get_value("rmax", val)) rmax = atof(val.c_str());
  if (get_value("rs", val)) rs = atof(val.c_str());
  if (get_value("numr", val)) numr = atoi(val.c_str());
  if (get_value("cmap", val)) cmap = atoi(val.c_str());
  if (get_value("Lmax", val)) Lmax = atoi(val.c_str());
  if (get_value("nmax", val)) nmax = atoi(val.c_str());
  if (get_value("self_consistent", val)) {
    if (atoi(val.c_str())) self_consistent = true;
    else self_consistent = false;
  }

}

SphereTwoCenter::~SphereTwoCenter(void)
{
  delete ortho;
  delete exp_ej;
  delete exp_com;
  delete [] center;
}

void SphereTwoCenter::get_acceleration_and_potential(Component* curComp)
{
  cC = curComp;			// "Register" component
  nbodies = cC->Number();	// And compute number of bodies

  
  bool use_external1 = use_external;

				// Set center to Component center
  for (int k=0; k<3; k++) center[k] = component->center[k];
  if (use_external) exp_ej -> SetExternal();
  exp_ej -> get_acceleration_and_potential(cC);
  exp_ej -> ClearExternal();

				// Reset external force flag
  use_external  = use_external1;

				// Set center to Component center of mass
  if (component->com_system)
    for (int k=0; k<3; k++) center[k] = component->comI[k];
  else
    for (int k=0; k<3; k++) center[k] = component->com[k];

  if (use_external) exp_com ->  SetExternal();
  exp_com -> get_acceleration_and_potential(cC);
  exp_com -> ClearExternal();


  // Clear external potential flag
  use_external = false;
}

