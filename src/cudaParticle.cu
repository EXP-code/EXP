#include "cudaParticle.cuH"


void ParticleHtoD(const Particle & h, cudaParticle & d)
{
  d.mass = h.mass;
  for (int k=0; k<3; k++) {
    d.pos[k] = h.pos[k];
    d.vel[k] = h.vel[k];
    d.acc[k] = h.acc[k];
  }
  d.pot    = h.pot;
  d.potext = h.potext;
  d.level  = h.level;
  d.indx   = h.indx;
}

void ParticleDtoH(const cudaParticle & d, Particle & h)
{
  h.mass = d.mass;
  for (int k=0; k<3; k++) {
    h.pos[k] = d.pos[k];
    h.vel[k] = d.vel[k];
    h.acc[k] = d.acc[k];
  }
  h.pot    = d.pot;
  h.potext = d.potext;
  h.level  = d.level;
  h.indx   = d.indx;
}
