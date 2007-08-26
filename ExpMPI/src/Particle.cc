#include "expand.h"

#include <Particle.H>

Particle::Particle()
{
  // Initialize basic fields

  mass = pot = potext = 0.0;
  for (int k=0; k<3; k++)
    pos[k] = vel[k] = acc[k] = 0.0;
  level = 0;
  indx = 0;
}

Particle::Particle(const Particle &p)
{
  mass = p.mass;
  for (int k=0; k<3; k++) {
    pos[k] = p.pos[k];
    vel[k] = p.vel[k];
    acc[k] = p.acc[k];
  }
  pot = p.pot;
  potext = p.potext;
  iattrib = p.iattrib;
  dattrib = p.dattrib;
  level = p.level;
  indx = p.indx;
}


