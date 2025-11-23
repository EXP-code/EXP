// -*- C++ -*-

#include "cudaUtil.cuH"
#include "cudaReduce.cuH"

#include "Component.H"
#include "PeriodicBC.H"

// Global device symbols for CUDA kernel
//
__device__ __constant__
cuFP_t cudaBCoffset[3], cudaBCside[3];

__device__ __constant__
char cudaBC[3];

// Cuda implementation of bar force
//
__global__ void
userPeriodicBCKernel(dArray<cudaParticle> P, dArray<int> I,
		     int stride, PII lohi)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {
      
#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      
      // Loop over coordinates
      //
      cuFP_t delta;

      for (int k=0; k<3; k++) {

	// Ignore vacuum boundary dimensions
	//
	char T = cudaBC[k];
	if (T == 'v') continue;

	cuFP_t offset = cudaBCoffset[k];
	cuFP_t pos    = p.pos[k] + offset;
	cuFP_t L      = cudaBCside[k];

	// Reflection BC
	//
	if (T == 'r') {
	  if (pos < 0.0) {
	    delta = -pos - L*floor(-pos/L);
	    p.pos[k] = delta - offset;
	    p.vel[k] *= -1.0;
	  } 
	  if (pos >= L) {
	    delta = pos - L*floor(pos/L);
	    p.pos[k] =  L - delta - offset;
	    p.vel[k] *= -1.0;
	  }
	}

	// Periodic BC
	//
	if (T == 'p') {
	  if (pos < 0.0) {
	    p.pos[k] += L*floor(1.0+fabs(pos/L));
	    
	  }
	  if (pos >= L) {
	    p.pos[k] += - L*floor(fabs(pos/L));
	  }
	}
	
	// Sanity check
	//
	if (p.pos[k] < -offset || p.pos[k] >= L-offset) {
	  printf("Error in periodic BC on #%d: pos[%d]=%e not in [%e,%e)\n",
		 npart, k, p.pos[k], -offset, L - offset);
	}
      }
      // END coordinate loop
      
    }
    // END: particle index block

  }
  // END: stride loop

}


__global__
void testConstantsPeriodicBC(cuFP_t tnow)
{
  printf("-------------------------\n");
  printf("---PeriodicBC constants--\n");
  printf("-------------------------\n");
  printf(" Time  = %e\n", tnow );
  printf(" Sides = %e, %e, %e\n",
	 cudaBCside[0], cudaBCside[1], cudaBCside[2] );
  printf("Offset = %e, %e, %e\n",
	 cudaBCoffset[0], cudaBCoffset[1], cudaBCoffset[2] );
  printf(" BCs   = %c, %c, %c\n", cudaBC[0], cudaBC[1], cudaBC[2] );
  printf("-------------------------\n");
}

void PeriodicBC::cuda_initialize()
{
  cuFP_t vec[3];
  for (int k=0; k<3; k++) vec[k] = offset[k];

  cuda_safe_call(cudaMemcpyToSymbol(cudaBCoffset, vec, sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cudaBCoffset");
  
  for (int k=0; k<3; k++) vec[k] = L[k];

  cuda_safe_call(cudaMemcpyToSymbol(cudaBCside, vec, sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cudaBCside");
  
  char h_cudaBC[3];
  for (int k=0; k<3; k++) h_cudaBC[k] = bc[k];

  cuda_safe_call(cudaMemcpyToSymbol(cudaBC, h_cudaBC, sizeof(char)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cudaBC");
  
  if (myid==0 and VERBOSE>4) {
    auto cr = cC->cuStream;
    testConstantsPeriodicBC<<<1, 1, 0, cr->stream>>>(tnow);
    cudaDeviceSynchronize();
    cuda_check_last_error_mpi("cudaDeviceSynchronize", __FILE__, __LINE__, myid);
  }
}


void PeriodicBC::determine_acceleration_and_potential_cuda()
{
  // Sanity check
  //
  int nbodies = cC->Number();
  if (nbodies != static_cast<int>(cC->Particles().size())) {
    std::cerr << "PeriodicBC: ooops! number=" << nbodies
	      << " but particle size=" << cC->Particles().size() << endl;
    nbodies = static_cast<int>(cC->Particles().size());
  }
  
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cC->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  // Stream structure iterators
  //
  auto cr = cC->cuStream;

  // Get particle index range for levels [mlevel, multistep]
  //
  PII lohi = cC->CudaGetLevelRange(mlevel, multistep);

  // Compute grid
  //
  unsigned int N         = lohi.second - lohi.first;
  unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
  if (N>0) {

    if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

    // Do the work
    //
    userPeriodicBCKernel<<<gridSize, BLOCK_SIZE, 0, cr->stream>>>
      (toKernel(cr->cuda_particles), toKernel(cr->indx1), stride, lohi);
  }
}
