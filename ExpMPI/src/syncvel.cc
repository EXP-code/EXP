/*
  Synchronize particle coordinates when outputing
  particle data to body data file or when computing diagnostics
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void synchronize_velocity(int sign)
{
  double rsign;

  if (sign)
    rsign = -1.0;
  else
    rsign = 1.0;

  list<Component*>::iterator cc;
  Component *c;
  unsigned ntot;
  
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    
    ntot = c->Number();
    for (int i=0; i<ntot; i++) {

      for (int k=0; k<c->dim; k++) 
	c->AddVel(i, k, rsign*(c->Acc(i, k) - c->acc0[k])*0.5*dtime );
    }
    
    if (c->com_system) {
      for (int k=0; k<c->dim; k++) c->cov0[k] += rsign*c->acc0[k]*0.5*dtime;
    }

  }
  
  // Increment times

  tvel=tvel+rsign*0.5*dtime;
  tnow=tvel;

}
