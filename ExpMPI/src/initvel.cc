/*
  Initialize the velocities of the bodies for the initial timestep
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void init_velocity(void)
{
  list<Component*>::iterator c;
  vector<Particle>::iterator p, pend;

  for (c=comp.components.begin(); c != comp.components.end(); c++) {

    pend = (*c)->particles.end();
    for (p=(*c)->particles.begin(); p != pend; p++) {

      for (int k=0; k<(*c)->dim; k++) p->vel[k] += 0.5*dtime*p->acc[k];
    }
  }

  // Increment velocity time, system time

  tvel=tvel+0.5*dtime;
  tnow=tvel;
}
