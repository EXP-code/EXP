#include "expand.h"

#include <Particle.H>

Particle::Particle()
{
  // Initialize basic fields

  mass = pot = potext = 0.0;
  for (int k=0; k<3; k++)
    pos[k] = vel[k] = acc[k] = 0.0;
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

}

bool Particle::freeze(void)
{
  double r2 = 0.0;

  for (int i=0; i<3; i++) r2 += pos[i]*pos[i];
  if (r2 > rmax_tidal*rmax_tidal) 
    {
      pos[0] = pos[1] = pos[2] = 1.0e10;
      return true;
    }
  else return false;
}

