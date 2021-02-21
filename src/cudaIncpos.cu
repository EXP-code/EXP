// -*- C++ -*-

#include <Component.H>
#include "expand.h"
#include <cudaUtil.cuH>
#include "cudaParticle.cuH"

#include <boost/make_shared.hpp>

/*

__global__ void coordDrift
(dArray<cudaParticle> in, cuFP_t dt, int dim, int stride, PII lohi)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Particle counter
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {

#ifdef BOUNDS_CHECK
      if (npart>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle p = in._v[npart];
    
      for (int k=0; k<dim; k++) p.pos[k] += p.vel[k]*dt;
    }
  }
}
*/


struct cudaIncPos : public thrust::unary_function<cudaParticle, cudaParticle>
{
  const cuFP_t _dt;
  const int    _dim;

  cudaIncPos(cuFP_t dt, int dim) : _dt(dt), _dim(dim) { }

  __host__ __device__
  cudaParticle operator()(cudaParticle& p)
  {
    for (int k=0; k<_dim; k++) p.pos[k] += p.vel[k]*_dt;
    return p;
  }
};

void incr_position_cuda(cuFP_t dt, int mlevel)
{
  for (auto c : comp->components) {

    auto cr = c->cuStream;

    if (multistep) {

      auto ret = c->CudaSortByLevel(cr, mlevel, multistep);
      
      thrust::transform(thrust::cuda::par.on(cr->stream),
			cr->cuda_particles.begin()+ret.first, cr->cuda_particles.end(),
			cr->cuda_particles.begin()+ret.first, cudaIncPos(dt, c->dim));
    } else {
      thrust::transform(thrust::cuda::par.on(cr->stream),
			cr->cuda_particles.begin(), cr->cuda_particles.end(),
			cr->cuda_particles.begin(), cudaIncPos(dt, c->dim));
    }
  }

  /*
    // Sort particles and get size
    //
    PII lohi = c->CudaSortByLevel(cr, mlevel, multistep);

    // Compute grid
    //
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
    if (N>0) {

      if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

      unsigned int Nthread = gridSize*BLOCK_SIZE;

      // Do the work
      //
      coordDrift<<<gridSize, BLOCK_SIZE>>>
	(toKernel(c->cuStream->cuda_particles), dt, c->dim, stride, lohi);
    }
  }
  */

}
