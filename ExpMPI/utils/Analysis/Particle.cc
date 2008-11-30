#include "Particle.h"

Particle::Particle()
{
  // Initialize basic fields
  //
  level = 0;
  mass = 0.0;
  for (int k=0; k<3; k++)
    pos[k] = vel[k] = 0.0;
}

Particle::Particle(const Particle &p)
{
  level = p.level;
  mass = p.mass;
  for (int k=0; k<3; k++) {
    pos[k] = p.pos[k];
    vel[k] = p.vel[k];
  }
}
