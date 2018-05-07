#ifndef PARTICLE_CUH
#define PARTICLE_CUH

#include <cuda.h>
#include <cuda_runtime.h>

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>

#include <cudaUtil.H>
#include <cudaTimer.H>

#include <Particle.H>

struct cudaParticle
{
  float mass;
  float pos[3];
  float vel[3];
  float acc[3];
  float pot;
  float potext;
  unsigned level;
  unsigned indx;
} __attribute__((aligned));

template <typename T>
struct LessCudaSeq
{
  __host__ __device__
  bool operator()(const T& p1, const T& p2) const
  {
    return p1.indx < p2.indx;
  }
};

template <typename T>
struct LessCudaLev
{
  __host__ __device__
  bool operator()(const T& p1, const T& p2) const
  {
    return p1.level < p2.level;
  }
};


void ParticleHtoD(const Particle & h, cudaParticle & d);

void ParticleDtoH(const cudaParticle & d, Particle & h);

#endif
