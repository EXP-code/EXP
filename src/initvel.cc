/*
  Initialize the velocities of the bodies for the initial timestep
*/

#include "expand.h"

void init_velocity(void)
{
  unsigned ntot;

  if (eqmotion) {

    for (auto c : comp.components) {

      ntot = c->Number();
      for (int i=0; i<ntot; i++) {

	Particle *p = c->Part(i);

	for (int k=0; k<c->dim; k++) 
	  p->vel[k] += 0.5*dtime*p->acc[k];
      }
      
      if (c->com_system) {
	for (int k=0; k<c->dim; k++) c->cov0[k] += 0.5*dtime*c->acc0[k];
      }
      
    }
    
  }

  // Increment velocity time, system time
}
