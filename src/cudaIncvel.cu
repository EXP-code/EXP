#include <Component.H>
#include <expand.h>
#include <cudaUtil.cuH>
#include <cudaReduce.cuH>
#include <cudaParticle.cuH>

#include <boost/make_shared.hpp>

#define USE_KERNEL

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
      cudaParticle & p = in._v[npart];
    
      for (int k=0; k<dim; k++) p.vel[k] += p.acc[k]*dt;
    }
  }
}


__global__ void velocityDebug
(dArray<cudaParticle> in, int stride, PII lohi)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Particle counter
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second and npart < in._s) {

      cudaParticle p = in._v[npart];
    
      printf("%d vel a=(%13.6e %13.6e %13.6e) p=%13.6e\n", i, p.acc[0], p.acc[1], p.acc[2], p.pot);
    }
  }
}

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

#ifndef USE_KERNEL		// Use Thrust by default

    auto lo = cr->cuda_particles.begin();
    auto hi = cr->cuda_particles.end();

    if (multistep) {
      auto ret = c->CudaGetLevelRange(mlevel, multistep);

      lo += ret.first;
      
      /* DEBUGGING info
      std::cout << "incVel load <" << c->name << "> myid="
		<< myid << " at level=" << mlevel
		<< ": " << ret.second - ret.first << std::endl;
      */
    }

    if (thrust::distance(lo, hi) >= 0) {

      thrust::transform(// thrust::cuda::par.on(cr->stream),
			thrust::cuda::par, lo, hi, lo, cudaIncVel(dt, c->dim));
    }
#else
    PII lohi = {0, cr->cuda_particles.size()};

    if (multistep) {		// Get particle range
      lohi = c->CudaGetLevelRange(mlevel, multistep);
    }

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, c->cudaDevice);

    // Compute grid
    //
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
    if (N>0) {
      
      if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

      // Do the work
      //
      velocityKick<<<gridSize, BLOCK_SIZE>>>
	(toKernel(c->cuStream->cuda_particles), dt, c->dim, stride, lohi);
    }
#endif

    // DEBUGGING output
    if (false) {
      PII lohi(0, std::min<int>(5, cr->cuda_particles.size()));

      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, c->cudaDevice);

      // Compute grid
      //
      unsigned int N         = lohi.second - lohi.first;
      unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
      unsigned int gridSize  = N/BLOCK_SIZE/stride;
      
      if (N>0) {
	
	if (N > gridSize*BLOCK_SIZE*stride) gridSize++;
	
	// Do the work
	//
	velocityDebug<<<gridSize, BLOCK_SIZE>>>
	  (toKernel(cr->cuda_particles), stride, lohi);
      }
    }
    // END: DEBUG
  }
  // END: component loop
}

// -*- C++ -*-

