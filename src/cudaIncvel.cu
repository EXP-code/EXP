// -*- C++ -*-

#include <Component.H>
#include "expand.h"
#include <cudaUtil.cuH>
#include "cudaParticle.cuH"

#include <boost/make_shared.hpp>

/*
__global__ void velocityKick
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
    
      for (int k=0; k<dim; k++) p.vel[k] += p.acc[k]*dt;
    }
  }
}
*/

struct cudaIncVel : public thrust::unary_function<cudaParticle, cudaParticle>
{
  const cuFP_t _dt;
  const int    _dim;

  cudaIncVel(cuFP_t dt, int dim) : _dt(dt), _dim(dim) { }

  __host__ __device__
  cudaParticle operator()(cudaParticle& p)
  {
    for (int k=0; k<_dim; k++) p.vel[k] += p.acc[k]*_dt;
    return p;
  }
};

void incr_velocity_cuda(cuFP_t dt, int mlevel)
{
  for (auto c : comp->components) {

    auto cr = c->cuStream;

    if (multistep) {

      auto ret = c->CudaGetLevelRange(cr, mlevel, multistep);
      
      std::cout << "[" << myid << ", " << mlevel << "]: #="
		<< ret.second - ret.first << std::endl;

      thrust::transform(// thrust::cuda::par.on(cr->stream),
			thrust::cuda::par,
			cr->cuda_particles.begin()+ret.first, cr->cuda_particles.end(),
			cr->cuda_particles.begin()+ret.first, cudaIncVel(dt, c->dim));
    } else {
      thrust::transform(thrust::cuda::par.on(cr->stream),
			cr->cuda_particles.begin(), cr->cuda_particles.end(),
			cr->cuda_particles.begin(), cudaIncVel(dt, c->dim));
    }
  }

  /*

    // Sort particles and get size
    //
    PII lohi = c->CudaGetLevelRange(cr, mlevel, multistep);

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
      velocityKick<<<gridSize, BLOCK_SIZE>>>
	(toKernel(cr->cuda_particles), dt, c->dim, stride, lohi);
    }
  }
  */
}
