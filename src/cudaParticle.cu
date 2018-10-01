#include <iostream>
#include <iomanip>
#include "cudaParticle.cuH"


void ParticleHtoD(PartPtr h, cudaParticle & d, int beg=0, int end=0)
{
  d.mass = h->mass;
  for (int k=0; k<3; k++) {
    d.pos[k] = h->pos[k];
    d.vel[k] = h->vel[k];
    d.acc[k] = h->acc[k];
  }
  d.pot    = h->pot;
  d.potext = h->potext;
#ifdef DATTRIB_CUDA
  if (end<h->dattrib.size()) {
    if (int n=beg; n<end; n++) d.datr[n-nbeg] = h->dattrib[n];
  } else {
    std::cerr << "Wrong attribute size in ParticleHtoD" << std::endl;
  }
#endif
  d.level  = h->level;
  d.indx   = h->indx;
}

void ParticleDtoH(const cudaParticle & d, PartPtr h, int beg=0, int end=0)
{
  h->mass = d.mass;
  for (int k=0; k<3; k++) {
    h->pos[k] = d.pos[k];
    h->vel[k] = d.vel[k];
    h->acc[k] = d.acc[k];
  }
  h->pot    = d.pot;
  h->potext = d.potext;
#ifdef DATTRIB_CUDA
  if (int n=beg; n<end; n++) h->dattrib[n] = d.datr[n-nbeg];
#endif
  h->level  = d.level;
  h->indx   = d.indx;
}

__host__
std::ostream& operator<< (std::ostream& os, const cudaParticle& p)
{
  std::streamsize sp = os.precision();
  os.precision(6);
  os << std::setw(10) << p.indx
     << std::setw( 4) << p.level
     << std::setw(16) << p.mass;
  for (int k=0; k<3; k++) os << std::setw(16) << p.pos[k];
  // for (int k=0; k<3; k++) os << std::setw(16) << p.vel[k];
  for (int k=0; k<3; k++) os << std::setw(16) << p.acc[k];
  // os << std::setw(16) << p.pot << std::setw(16) << p.potext;
  os.precision(sp);
  return os;
}
