#include "expand.h"

#include <Particle.H>

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

