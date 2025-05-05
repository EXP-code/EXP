// -*- C++ -*-

#ifndef PARTICLE_CUH
#define PARTICLE_CUH

// Sets the fixed default number of real attributes.  This currently
// defined by CMake
//
// #define DATTRIB_CUDA 4

#include <cuda.h>
#include <cuda_runtime.h>

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>

#include <config_exp.h>

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
  cuFP_t  scale;
#if DATTRIB_CUDA>0
  cuFP_t datr[DATTRIB_CUDA];
#endif
  cuFP_t dtreq;
  unsigned indx;
  int lev[2];
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
  int _t;

  __host__ __device__
  LessCudaLev(int t=0) : _t(t) { }

  __host__ __device__
  bool operator()(const cudaParticle& p1, const cudaParticle& p2) const
  {
    return p1.lev[_t] < p2.lev[_t];
  }
};

//! For thrust sorts on on level and update
struct LessCudaLev2
{
  __host__ __device__
  bool operator()(const cudaParticle& p1, const cudaParticle& p2) const
  {
    return (p1.lev[0] < p2.lev[0]) or (p1.lev[0] == p2.lev[0] and p1.lev[1] < p2.lev[1]);
  }
};


//! Copy EXP particle to CUDA particle structure.  Attributes will not
//! be loaded if end==0.
int ParticleHtoD(PartPtr h, cudaParticle & d, int beg=0, int end=0);

//! Copy CUDA particle structure to EXP particle.  Attributes will not
//! be loaded if end==0.
void ParticleDtoH(const cudaParticle & d, PartPtr h, int beg=0, int end=0);

//! Print cudaParticle to output stream
extern __host__
std::ostream& operator<< (std::ostream& os, const cudaParticle& p);

//! Functor for extracting level info from cudaParticle structures
#if CUDART_VERSION < 12040
struct cuPartToLevel :
  public thrust::unary_function<cudaParticle, unsigned int>
#else
struct cuPartToLevel
#endif
{
  int _t;

  __host__ __device__
  cuPartToLevel(int t=0) : _t(t) {}

  __host__ __device__
  unsigned operator()(const cudaParticle& p) const
  {
    return p.lev[_t];
  }
};

//! Functor for extracting change level info from cudaParticle structures
#if CUDART_VERSION < 12040
struct cuPartToChange :
  public thrust::unary_function<cudaParticle, thrust::pair<int, int>>
#else
struct cuPartToChange
#endif
{
  __host__ __device__
  thrust::pair<int, int> operator()(const cudaParticle& p) const
  {
    return {p.lev[0], p.lev[1]};
  }
};

#endif
