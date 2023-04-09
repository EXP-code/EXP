// -*- C++ -*-

#include <cudaUtil.cuH>
#include <cudaReduce.cuH>

#include <Component.H>
#include "UserAgnNoise.H"

// Global device symbols for CUDA kernel
//
__device__ __constant__
cuFP_t userAgnR0, userAgnTau2, userAgnEps;

__device__ __constant__
int userAgnLoc;

// Cuda implementation of AGN mass update
//
__global__ void
userAgnNoiseKernel(dArray<cudaParticle> P, int stride,
		   cuFP_t tnow, cuFP_t tev)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int npart = tid*stride + n;	// Index in the stride

    if (npart < P._s) {
      
      cudaParticle & p = P._v[npart];
      
      if (tnow > tev) {
	cuFP_t rr = 0.0;
	for (int k=0; k<3; k++) rr += p.pos[k]*p.pos[k];
	if (rr < userAgnR0*userAgnR0) {
	  p.datr[userAgnLoc+1] = tnow;
	}
      }

      p.mass = p.datr[userAgnLoc]*(1.0 - userAgnEps*exp(-(tnow - p.datr[userAgnLoc+1])/userAgnTau2));
      
    } // Particle index block

  } // END: stride loop

}


__global__
void testConstantsAgnNoise(cuFP_t tnow)
{
  printf("-------------------------\n");
  printf("---UserAgn constants-----\n");
  printf("-------------------------\n");
  printf("   Time   = %e\n", tnow        );
  printf("   R0     = %e\n", userAgnR0   );
  printf("   Tau2   = %e\n", userAgnTau2 );
  printf("   eps    = %e\n", userAgnEps  );
  printf("   loc    = %d\n", userAgnLoc  );
  printf("-------------------------\n");
}


void UserAgnNoise::determine_acceleration_and_potential_cuda()
{
  // Sanity check
  //
  int nbodies = cC->Number();
  if (nbodies != static_cast<int>(cC->Particles().size())) {
    std::cerr << "UserAgnNoise: ooops! number=" << nbodies
	      << " but particle size=" << cC->Particles().size() << endl;
    nbodies = static_cast<int>(cC->Particles().size());
  }
  
  if (nbodies==0) {		// Return if there are no particles
    if (VERBOSE>4) {
      std::cout << "Process " << myid << ": in UserAgnNoise, nbodies=0" 
		<< " for Component <" << cC->name << "> at T=" << tnow
		<< std::endl;
    }
    return;
  }

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cC->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  // Stream structure iterators
  //
  auto cr = cC->cuStream;

  // VERBOSE diagnostic output on first ncheck calls
  //
  static int cntr = 0;
  const  int ncheck = 10;
  
  if (cntr==0) {

    cuda_safe_call(cudaMemcpyToSymbol(userAgnR0, &R0, sizeof(cuFP_t),
				      size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying userAgnR0");
  
    cuda_safe_call(cudaMemcpyToSymbol(userAgnTau2, &tau2, sizeof(cuFP_t),
				      size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying userAgnTau2");
  
    cuda_safe_call(cudaMemcpyToSymbol(userAgnEps, &eps, sizeof(cuFP_t),
				      size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying userSatCen");
    
    cuda_safe_call(cudaMemcpyToSymbol(userAgnLoc, &loc, sizeof(int),
				      size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying userSatPos");
  }
    
  if (cntr < ncheck) {
    if (myid==0) {
      testConstantsAgnNoise<<<1, 1, 0, cr->stream>>>(tnow);
      cudaDeviceSynchronize();
      cuda_check_last_error_mpi("cudaDeviceSynchronize", __FILE__, __LINE__, myid);
    }
    cntr++;
  }

  // Compute grid
  //
  unsigned int N         = cr->cuda_particles.size();
  unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
  if (N>0) {

    if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

    // Do the work
    //
    userAgnNoiseKernel<<<gridSize, BLOCK_SIZE, 0, cr->stream>>>
      (toKernel(cr->cuda_particles), stride, tnow, tev);
  }
}
