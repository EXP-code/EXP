#include "expand.h"

#include <Particle.H>

float Particle::effort_default = 1.0e-12;

Particle::Particle()
{
  // Initialize basic fields

  mass = pot = potext = 0.0;
  for (int k=0; k<3; k++)
    pos[k] = vel[k] = acc[k] = 0.0;
  level   = 0;
  dtreq   = -1;
  scale   = -1;
  effort  = effort_default;
  indx    = 0;
  tree    = 0u;
  key     = 0u;
}

Particle::Particle(unsigned niatr, unsigned ndatr)
{
  // Initialize basic fields

  mass = pot = potext = 0.0;
  for (int k=0; k<3; k++)
    pos[k] = vel[k] = acc[k] = 0.0;
  level   = 0;
  dtreq   = -1;
  scale   = -1;
  effort  = effort_default;
  indx    = 0;
  tree    = 0u;
  key     = 0u;
  iattrib = vector<int>(niatr, 0);
  dattrib = vector<double>(ndatr, 0);
}

Particle::Particle(const Particle &p)
{
  mass = p.mass;
  for (int k=0; k<3; k++) {
    pos[k] = p.pos[k];
    vel[k] = p.vel[k];
    acc[k] = p.acc[k];
  }
  pot     = p.pot;
  potext  = p.potext;
  iattrib = p.iattrib;
  dattrib = p.dattrib;
  level   = p.level;
  dtreq   = p.dtreq;
  scale   = p.scale;
  effort  = p.effort;
  indx    = p.indx;
  tree    = p.tree;
  key     = p.key;
}


