#ifndef PARTICLE_CUH
#define PARTICLE_CUH

#include <cuda.h>
#include <cuda_runtime.h>

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>

#include <cudaUtil.cuH>

#include <Particle.H>

struct cudaParticle
{
  cuFP_t mass;
  cuFP_t pos[3];
  cuFP_t vel[3];
  cuFP_t acc[3];
  cuFP_t pot;
  cuFP_t potext;
  unsigned level;
  unsigned indx;
} __attribute__((aligned));

struct LessCudaSeq
{
  __host__ __device__
  bool operator()(const cudaParticle& p1, const cudaParticle& p2) const
  {
    return p1.indx < p2.indx;
  }
};

struct LessCudaLev
{
  __host__ __device__
  bool operator()(const cudaParticle& p1, const cudaParticle& p2) const
  {
    return p1.level < p2.level;
  }
};


void ParticleHtoD(const Particle & h, cudaParticle & d);

void ParticleDtoH(const cudaParticle & d, Particle & h);

extern __host__
std::ostream& operator<< (std::ostream& os, const cudaParticle& p);

#endif
