// -*- C++ -*-

#include <cudaUtil.cuH>
#include <cudaReduce.cuH>

#include <Component.H>
#include "UserSat.H"

// Global device symbols for CUDA kernel
//
__device__ __constant__
cuFP_t userBarAmp, userBarPosAng, userBarCen[3], userBarB5;

__device__ __constant__
int userBarSoft;


// Cuda implementation of bar force
//
__global__ void
userBarForceKernel(dArray<cudaParticle> P, dArray<int> I,
		   int stride, PII lohi)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Particle index

    // Temporary relative position
    cuFP_t fpos[3];

    // Angular dependence from position angle
    cuFP_t cos2p = cos(2.0*userBarPosAng);
    cuFP_t sin2p = sin(2.0*userBarPosAng);

    if (npart < lohi.second) {
      
#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      
      // Compute position relative to bar
      cuFP_t rr = 0.0, fac, ffac, nn;
      for (int k=0; k<3; k++) {
	fpos[k] = p.pos[k] - userBarCen[k];
	rr += fpos[k] * fpos[k];
      }
      rr = sqrt(rr);

      cuFP_t xx = fpos[0];
      cuFP_t yy = fpos[1];
      cuFP_t zz = fpos[2];
      cuFP_t pp = (xx*xx - yy*yy)*cos2p + 2.0*xx*yy*sin2p;
    
      if (userBarSoft) {
	fac = 1.0 + rr/b5;
	ffac = -userBarAmp / pow(fac, 6.0);
	nn = pp / (b5*rr);
      } else {
	fac = 1.0 + pow(rr/b5, 5.0);
	ffac = -userBarAmp / (fac*fac);
	nn = pp * pow(rr/b5, 3.0)/ (b5*b5);
      }

      // Add acceleration
      p.acc[0] += ffac * ( 2.0*( xx*cos2p + yy*sin2p)*fac - 5.0*nn*xx);
      p.acc[1] += ffac * ( 2.0*(-yy*cos2p + xx*sin2p)*fac - 5.0*nn*yy);
      p.acc[2] += ffac * ( -5.0*nn*zz);

      // Add external potential
      p.potext += -ffac*pp*fac;
      
    } // Particle index block

  } // END: stride loop

}


__global__
void testConstantsUserSat(cuFP_t tnow)
{
  printf("-------------------------\n");
  printf("---UserSat constants-----\n");
  printf("-------------------------\n");
  printf("   Time   = %e\n", tnow          );
  printf("   Amp    = %e\n", userBarAmp    );
  printf("   Amp    = %e\n", userBarNumFac );
  printf("   b5     = %e\n", userBarB5     );
  printf("   Center = %e, %e, %e\n",
	 userBarCen[0], userBarCen[1], userBarCen[2] );
  if (userBarSoft) printf("   Soft   = true\n" );
  else             printf("   Soft   = false\n");
  printf("-------------------------\n");
}

void UserBar::determine_acceration_and_potential_cuda()
{
  // Sanity check
  //
  int nbodies = cC->Number();
  if (nbodies != static_cast<int>(cC->Particles().size())) {
    std::cerr << "UserSat: ooops! number=" << nbodies
	      << " but particle size=" << cC->Particles().size() << endl;
    nbodies = static_cast<int>(cC->Particles().size());
  }
  
  if (nbodies==0) {		// Return if there are no particles
    if (verbose) {
      cout << "Process " << myid << ": in UserBar, nbodies=0" 
	   << " for Component <" << cC->name << "> at T=" << tnow
	   << endl;
    }
    return;
  }

  double amp = afac * numfac * amplitude/fabs(amplitude) 
    * 0.5*(1.0 + erf( (tnow - Ton )/DeltaT ))
    * 0.5*(1.0 - erf( (tnow - Toff)/DeltaT )) ;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cC->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  // Stream structure iterators
  //
  auto cr = cC->cuStream;

  // Assign expansion center
  //
  auto cn = cC->getCenter();
  cuFP_t ctr[3], dtmp;
  for (int k=0; k<3; k++) ctr[k] = cn[k];

  cuFP_t dtmp;

  cuda_safe_call(cudaMemcpyToSymbol(userBarAmp, &(dtmp=amp), sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userBarAmp");
  
  cuda_safe_call(cudaMemcpyToSymbol(userBarPosAng, &(dtmp=posang), sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userBarPosAng");
  
  cuda_safe_call(cudaMemcpyToSymbol(userBarCen, ctr, sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userBarCen");

  cuda_safe_call(cudaMemcpyToSymbol(userBarB5, &(dtmp=b5), sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userBarB5");

  int cuSoft = 0;
  if (soft) cuSoft = 1;
  cuda_safe_call(cudaMemcpyToSymbol(userBarSoft, &cuSoft,  sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userBarSoft");

  // VERBOSE diagnostic output on first ncheck calls
  //
  static int cntr = 0;
  const  int ncheck = 100;
  
  if (myid==0 and VERBOSE>4 and cntr < ncheck) {
    testConstantsUserBar<<<1, 1, 0, cr->stream>>>(tnow);
    cudaDeviceSynchronize();
    cuda_check_last_error_mpi("cudaDeviceSynchronize", __FILE__, __LINE__, myid);
    cntr++;
  }

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
    userBarForceKernel<<<gridSize, BLOCK_SIZE, 0, cr->stream>>>
      (toKernel(cr->cuda_particles), toKernel(cr->indx1), stride, lohi);
  }
}
