// -*- C++ -*-

#include <Component.H>
#include <expand.H>
#include <cudaUtil.cuH>
#include <cudaReduce.cuH>
#include <cudaParticle.cuH>

__global__ void coordDrift
(dArray<cudaParticle> P, dArray<int> I, cuFP_t dt, int dim, int stride, PII lohi)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Particle counter
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {

#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
    
      for (int k=0; k<dim; k++) p.pos[k] += p.vel[k]*dt;
    }
  }
}

__global__ void positionDebug
(dArray<cudaParticle> P, dArray<int> I, int stride, PII lohi)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Particle counter
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second and npart < I._s) {

      cudaParticle & p = P._v[I._v[npart]];
      cuFP_t sumP = 0.0, sumV = 0.0, sumA = 0.0;
      for (int k=0; k<3; k++) {
	sumP += p.pos[k]*p.pos[k];
	sumV += p.vel[k]*p.vel[k];
	sumA += p.acc[k]*p.acc[k];
      }
      sumP = sqrt(sumP);
      sumV = sqrt(sumV);
      sumA = sqrt(sumA);
    
      printf("[%d, %d] pos r=%13.6e v=%13.6e a=%13.6e\n", i, p.indx,
	     sumP, sumV, sumA);
    }
  }
}


void incr_position_cuda(cuFP_t dt, int mlevel)
{
  for (auto c : comp->components) {

    auto cr = c->cuStream;

    PII lohi = {0, cr->cuda_particles.size()};

    if (multistep) {		// Get particle range
      lohi = c->CudaGetLevelRange(mlevel, mlevel);
    }

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, c->cudaDevice);
    cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

    // Compute grid
    //
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
    if (N>0) {
      
      if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

      // Do the work
      //
      coordDrift<<<gridSize, BLOCK_SIZE>>>
	(toKernel(cr->cuda_particles),
	 toKernel(cr->indx1), dt, c->dim, stride, lohi);
    }

    // DEBUGGING output
    //
    if (false) {
      PII lohi(0, std::min<int>(3, cr->cuda_particles.size()));

      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, c->cudaDevice);
      cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

      // Compute grid
      //
      unsigned int N         = lohi.second - lohi.first;
      unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
      unsigned int gridSize  = N/BLOCK_SIZE/stride;
      
      if (N>0) {
	
	if (N > gridSize*BLOCK_SIZE*stride) gridSize++;
	
	// Do the work
	//
	positionDebug<<<gridSize, BLOCK_SIZE>>>
	  (toKernel(cr->cuda_particles), toKernel(cr->indx1), stride, lohi);
      }
    }
    // END: DEBUG

    c->CudaToParticles();
  }
  // END: component loop
}
