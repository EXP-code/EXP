// -*- C++ -*-

#include <iostream>
#include <iomanip>
#include "cudaParticle.cuH"


int ParticleHtoD(PartPtr h, cudaParticle & d, int beg, int end)
{
  d.mass = h->mass;
  for (int k=0; k<3; k++) {
    d.pos[k] = h->pos[k];
    d.vel[k] = h->vel[k];
    d.acc[k] = h->acc[k];
  }
  d.pot    = h->pot;
  d.potext = h->potext;
  d.scale  = h->scale;
#if DATTRIB_CUDA>0
  // Skip attributes if end is 0
  //
  if (end) {
    if (end - beg < DATTRIB_CUDA) {
      if (end <= h->dattrib.size()) {
	for (int n=beg; n<end; n++) d.datr[n-beg] = h->dattrib[n];
      } else {
	std::cerr << "Wrong attribute size in ParticleHtoD" << std::endl;
	return 1;
      }
    } else {
	std::cerr << "Requested number of attributes ["
		  << end - beg << "] in ParticleHtoD is larger than storage="
		  << DATTRIB_CUDA << std::endl;
	return 1;
    }
  }
#endif
  d.dtreq    = h->dtreq;
  d.lev[0]   = h->level;
  d.lev[1]   = h->level;
  d.indx     = h->indx;

  return 0;
}

void ParticleDtoH(const cudaParticle & d, PartPtr h, int beg, int end)
{
  h->mass = d.mass;
  for (int k=0; k<3; k++) {
    h->pos[k] = d.pos[k];
    h->vel[k] = d.vel[k];
    h->acc[k] = d.acc[k];
  }
  h->pot    = d.pot;
  h->potext = d.potext;
  h->scale  = d.scale;
#if DATTRIB_CUDA>0
  if (end) {
    if (end - beg < DATTRIB_CUDA) {
      if (end <= h->dattrib.size()) {
	for (int n=beg; n<end; n++) h->dattrib[n] = d.datr[n-beg];
      } else {
	std::cerr << "Wrong attribute size in ParticleDtoH" << std::endl;
      }
    } else {
	std::cerr << "Requested number of attributes ["
		  << end - beg << "] in ParticleDtoH is larger than storage="
		  << DATTRIB_CUDA << std::endl;
    }
  }
#endif
  h->dtreq  = d.dtreq;
  h->level  = d.lev[0];
  h->indx   = d.indx;
}

__host__
std::ostream& operator<< (std::ostream& os, const cudaParticle& p)
{
  std::streamsize sp = os.precision();
  os.precision(6);
  // Index, levels, mass
  os << std::setw(10) << p.indx
     << std::setw( 4) << p.lev[0]
     << std::setw( 4) << p.lev[1]
     << std::setw(16) << p.mass;

  // Position
  for (int k=0; k<3; k++) os << std::setw(16) << p.pos[k];

  // Velocity
  for (int k=0; k<3; k++) os << std::setw(16) << p.vel[k];

  // Acceleration
  for (int k=0; k<3; k++) os << std::setw(16) << p.acc[k];

  // Potential
  os << std::setw(16) << p.pot << std::setw(16) << p.potext;

#if DATTRIB_CUDA>0
  // Double attributes
  for (int n=0; n<DATTRIB_CUDA; n++) os << std::setw(16) << p.datr[n];
#endif

  os.precision(sp);
  return os;
}
