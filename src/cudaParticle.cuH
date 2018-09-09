// -*- C++ -*-

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

//! Simplified particle structure for use in CUDA kernel code
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

//! For thrust sorts on sequence
struct LessCudaSeq
{
  __host__ __device__
  bool operator()(const cudaParticle& p1, const cudaParticle& p2) const
  {
    return p1.indx < p2.indx;
  }
};

//! For thrust sorts on on level
struct LessCudaLev
{
  __host__ __device__
  bool operator()(const cudaParticle& p1, const cudaParticle& p2) const
  {
    return p1.level < p2.level;
  }
};


//! Copy EXP particle to CUDA particle structure
void ParticleHtoD(PartPtr h, cudaParticle & d);

//! Copy CUDA particle structure to EXP particle
void ParticleDtoH(const cudaParticle & d, PartPtr h);

//! Print cudaParticle to output stream
extern __host__
std::ostream& operator<< (std::ostream& os, const cudaParticle& p);

//! Functor for extracting level info from cudaParticle structures
struct cuPartToLevel :
  public thrust::unary_function<cudaParticle, unsigned>
{
  __host__ __device__
  unsigned operator()(const cudaParticle& p) const
  {
    return p.level;
  }
};

#endif
