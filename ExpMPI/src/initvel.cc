/*
  Initialize the velocities of the bodies for the initial timestep
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void init_velocity(void)
{
  list<Component*>::iterator cc;
  vector<Particle>::iterator p, pend;
  Component *c;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    pend = c->particles.end();
    for (p=c->particles.begin(); p != pend; p++) {

      for (int k=0; k<c->dim; k++) p->vel[k] += 0.5*dtime*p->acc[k];
    }

    if (c->com_system) {
      for (int k=0; k<c->dim; k++) c->cov0[k] += 0.5*dtime*c->acc0[k];
    }

  }

  // Increment velocity time, system time

  tvel=tvel+0.5*dtime;
  tnow=tvel;
}
