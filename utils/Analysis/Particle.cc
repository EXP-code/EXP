#include "Particle.h"

#ifdef STANDALONE
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

Particle::Particle()
{
  // Initialize basic fields
  //
  level = 0;
  mass = 0.0;
  for (int k=0; k<3; k++)
    pos[k] = vel[k] = 0.0;
  indx = 0;
}

Particle::Particle(const Particle &p)
{
  level = p.level;
  mass = p.mass;
  for (int k=0; k<3; k++) {
    pos[k] = p.pos[k];
    vel[k] = p.vel[k];
  }
  indx = p.indx;
}
