// -*- C++ -*-

#include "expand.H"
#include <Component.H>
#include <cudaReduce.cuH>

// Define this to see per operation sub timings
//
// #define VERBOSE_TIMING

// Global symbols for time step selection
//
__device__ __constant__
cuFP_t cuDynfracS, cuDynfracD, cuDynfracV, cuDynfracA, cuDynfracP, cuDtime;

__device__ __constant__
int cuMultistep, cuShiftlev, cuDTold;

void cuda_initialize_multistep_constants()
{
  // Copy constants to device
  //
  cuFP_t z;

  cuda_safe_call(cudaMemcpyToSymbol(cuDynfracS, &(z=dynfracS), sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuDynfracS");

  cuda_safe_call(cudaMemcpyToSymbol(cuDynfracD, &(z=dynfracD), sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuDynfracD");

  cuda_safe_call(cudaMemcpyToSymbol(cuDynfracV, &(z=dynfracV), sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuDynfracV");

  cuda_safe_call(cudaMemcpyToSymbol(cuDynfracA, &(z=dynfracA), sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuDynfracA");

  cuda_safe_call(cudaMemcpyToSymbol(cuDynfracP, &(z=dynfracP), sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuDynfracP");

  cuda_safe_call(cudaMemcpyToSymbol(cuDtime, &(z=dtime), sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuDtime");

  cuda_safe_call(cudaMemcpyToSymbol(cuMultistep, &multistep, sizeof(int), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuMultistep");

  cuda_safe_call(cudaMemcpyToSymbol(cuShiftlev, &shiftlevl, sizeof(int), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuShiftlev");

  int tmp = DTold ? 1 : 0;

  cuda_safe_call(cudaMemcpyToSymbol(cuDTold, &(tmp), sizeof(int), size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying cuDTold");
}


__global__
void testConstantsMultistep()
{
  printf("------------------------\n"     );
  printf("---Multistep constants--\n"   );
  printf("------------------------\n"   );
  printf("   DynS   = %e\n", cuDynfracS );
  printf("   DynD   = %e\n", cuDynfracD );
  printf("   DynV   = %e\n", cuDynfracV );
  printf("   DynA   = %e\n", cuDynfracA );
  printf("   DynP   = %e\n", cuDynfracP );
  printf("   Dtime  = %f\n", cuDtime    );
  printf("   Multi  = %d\n", cuMultistep);
  printf("   Shift  = %d\n", cuShiftlev );
  if (cuDTold)
    printf("   DTold  = true\n"         );
  else
    printf("   DTold  = false\n"        );
  printf("------------------------\n"   );
}

__global__ void
timestepSetKernel(dArray<cudaParticle> P, int stride, int mstep, int mdrft)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int npart = tid*stride + n;
    if (npart < P._s) P._v[npart].dtreq = 0.0;
  }
  // END: stride loop

}

__global__ void
timestepKernel(dArray<cudaParticle> P, dArray<int> I,
	       cuFP_t cx, cuFP_t cy, cuFP_t cz,
	       int minactlev, int dim, int stride, PII lohi,
	       bool switching, bool shift)
{
  const int tid    = blockDim.x * blockIdx.x + threadIdx.x;
  const cuFP_t eps = 1.0e-20;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Index into the sorted array

    if (npart < lohi.second) {
      
#ifdef BOUNDS_CHECK
      if (npart >= I._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      //                      ^    ^
      //                      |    |
      // Particle index ------+    |
      // Sort index ---------------+
      
      cuFP_t xx = p.pos[0] - cx;
      cuFP_t yy = p.pos[1] - cy;
      cuFP_t zz = p.pos[2] - cz;
      
      cuFP_t dtd=1.0/eps, dtv=1.0/eps, dta=1.0/eps, dtA=1.0/eps, dts=1.0/eps;

      if (cuDTold) {

	// dtv = eps* r/v         -- roughly, crossing time
	// dta = eps* v/a         -- force scale
	// dtA = eps* sqrt(r/a)   -- acceleration time
	
	cuFP_t rtot = sqrt(xx*xx + yy*yy + zz*zz);
	cuFP_t vtot = 0.0;
	cuFP_t atot = 0.0;

	for (int k=0; k<dim; k++) {
	  vtot += p.vel[k]*p.vel[k];
	  atot += p.acc[k]*p.acc[k];
	}
	vtot = sqrt(vtot) + 1.0e-18;
	atot = sqrt(atot) + 1.0e-18;
	
	if (p.scale>0.0) dts = cuDynfracS*p.scale/vtot;

	dtv = cuDynfracV*rtot/vtot;
	dta = cuDynfracA*vtot/atot;
	dtA = cuDynfracP*sqrt(rtot/atot);

      } else {

	// dtd = eps* rscale/v_i    -- char. drift time scale
	// dtv = eps* min(v_i/a_i)  -- char. force time scale
	// dta = eps* phi/(v * a)   -- char. work time scale
	// dtA = eps* sqrt(phi/a^2) -- char. "escape" time scale

	cuFP_t dtr  = 0.0;
	cuFP_t vtot = 0.0;
	cuFP_t atot = 0.0;
	
	for (int k=0; k<dim; k++) {
	  dtr  += p.vel[k]*p.acc[k];
	  vtot += p.vel[k]*p.vel[k];
	  atot += p.acc[k]*p.acc[k];
	}

	cuFP_t ptot = fabs(p.pot + p.potext);
	
	if (p.scale>0) dts = cuDynfracS*p.scale/fabs(sqrt(vtot)+eps);
	
	dtd = cuDynfracD * 1.0/sqrt(vtot+eps);
	dtv = cuDynfracV * sqrt(vtot/(atot+eps));
	dta = cuDynfracA * ptot/(fabs(dtr)+eps);
	dtA = cuDynfracP * sqrt(ptot/(atot+eps));
      }

      
      // Smallest time step
      //
      cuFP_t dt = dts;
      if (dt > dtd) dt = dtd;
      if (dt > dtv) dt = dtv;
      if (dt > dta) dt = dta;
      if (dt > dtA) dt = dtA;
      if (dt < eps) dt = eps;
      
      if (p.dtreq<=0) {
	p.dtreq = dt;
      } else {
	if (dt > p.dtreq) dt = p.dtreq;
	else              p.dtreq = dt;
      }

      if (switching or minactlev==0) {

	// Time step wants to be LARGER than the maximum
	//
	p.lev[1] = 0;
	if (dt<cuDtime)
	  p.lev[1] = (int)floor(log(cuDtime/dt)/log(2.0));
    
	// Time step wants to be SMALLER than the maximum
	//
	if (p.lev[1]>cuMultistep) p.lev[1] = cuMultistep;
      
	// Limit new level to minimum active level
	//
	if (p.lev[1]<minactlev)   p.lev[1] = minactlev;
	
	// Enforce n-level shifts at a time
	//
	if (cuShiftlev and shift) {
	  if (p.lev[1] > p.lev[0]) {
	    if (p.lev[1] - p.lev[0] > cuShiftlev)
	      p.lev[1] = p.lev[0] + cuShiftlev;
	  } else if (p.lev[0] > p.lev[1]) {
	    if (p.lev[0] - p.lev[1] > cuShiftlev)
	      p.lev[1] = p.lev[0] - cuShiftlev;
	  }
	}
      }
      // Active level = 0
	
    } // Particle index block
    
  } // END: stride loop

}

// Reset target level to current level
//
__global__ void
timestepFinalizeKernel(dArray<cudaParticle> P, dArray<int> I,
		       int stride, PII lohi)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Index into the sorted array

    if (npart < lohi.second) {
      
#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      
      if (p.lev[0] != p.lev[1]) p.lev[0] = p.lev[1];

    } // Particle index block
    
  } // END: stride loop

}

void cuda_initialize_multistep()
{
  // Constants to device once only
  //
  cuda_initialize_multistep_constants();

  // VERBOSE diagnostic output
  //
  if (myid==0 and VERBOSE>4)
    testConstantsMultistep<<<1, 1>>>();
}

void cuda_compute_levels()
{
  // DEBUGGING
  if (false) {
    for (int n=0; n<numprocs; n++) {
      if (n==myid) testConstantsMultistep<<<1, 1>>>();
      std::cout << std::endl;
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  // END DEBUGGING

  //
  // Begin the update
  //
  for (auto c : comp->components) {
    c->force->multistep_update_begin();
    if (not c->force->cudaAware()) c->ParticlesToCuda();
  }

  cudaDeviceProp deviceProp;

#ifdef VERBOSE_TIMING
  double time1 = 0.0, time2 = 0.0, timeSRT = 0.0, timeADJ = 0.0, timeCOM = 0.0;
  auto start0 = std::chrono::high_resolution_clock::now();
  auto start  = std::chrono::high_resolution_clock::now();
#endif

  for (auto c : comp->components) {
    
    cudaGetDeviceProperties(&deviceProp, c->cudaDevice);
    cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

    // Set maximum time step
    //
    if (multistep and c->NoSwitch() and c->DTreset() and mdrft==1) {

      // Compute grid
      //
      unsigned int N         = c->cuStream->cuda_particles.size();
      unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
      unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
      if (N>0) {
	if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

	timestepSetKernel<<<gridSize, BLOCK_SIZE>>>
	  (toKernel(c->cuStream->cuda_particles), stride, mstep, mdrft);
      }
    }

    PII lohi = {0, c->cuStream->cuda_particles.size()};
    if (multistep) lohi = c->CudaGetLevelRange(mfirst[mdrft], multistep);
      
    // DEBUGGING
    if (false and multistep>0) {
      for (int n=0; n<numprocs; n++) {
	if (n==myid) testConstantsMultistep<<<1, 1>>>();
	std::cout << std::string(60, '-') << std::endl
		  << "[" << myid << ", " << c->name
		  << "]: mlevel=" << mfirst[mdrft]
		  << " mstep=" << mstep << " mdrft=" << mdrft
		  << " (lo, hi) = (" << lohi.first << ", " << lohi.second << ")"
		  << std::endl << std::string(60, '-') << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    // END DEBUGGING

    // Compute grid
    //
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
    if (N>0) {

      if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

      // Get the current expansion center
      //
      auto ctr = c->getCenter(Component::Local | Component::Centered);
      
      // Allow intrastep level switches
      //
      bool switching = true;

      // Allowing switching on first call to set levels
      //
      bool notFirstCall = this_step!=0 or mstep!=0;
      if (c->NoSwitch() and notFirstCall) switching = false;

      // Do the work
      //
      timestepKernel<<<gridSize, BLOCK_SIZE>>>
	(toKernel(c->cuStream->cuda_particles),
	 toKernel(c->cuStream->indx1), ctr[0], ctr[1], ctr[2],
	 mfirst[mdrft], c->dim, stride, lohi, switching, notFirstCall);
    }
  }

#ifdef VERBOSE_TIMING
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::micro> duration = finish - start;
  time1 += duration.count()*1.0e-6;
#endif
  
  //
  // Finish the update
  //
  for (auto c : comp->components) {
#ifdef VERBOSE_TIMING
    start = std::chrono::high_resolution_clock::now();
#endif
    cudaGetDeviceProperties(&deviceProp, c->cudaDevice);
    cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);
    
    if (not c->force->NoCoefs()) c->force->multistep_update_cuda();

#ifdef VERBOSE_TIMING
    finish = std::chrono::high_resolution_clock::now();
    duration = finish - start;
    timeADJ += duration.count()*1.0e-6;
    start = std::chrono::high_resolution_clock::now();
#endif

    // Compute grid
    //
    PII lohi = {0, c->cuStream->cuda_particles.size()};
    if (multistep) lohi = c->CudaGetLevelRange(mfirst[mdrft], multistep);
      
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
    if (N>0) {

      if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

      timestepFinalizeKernel<<<gridSize, BLOCK_SIZE>>>
	(toKernel(c->cuStream->cuda_particles),
	 toKernel(c->cuStream->indx1), stride, lohi);

#ifdef VERBOSE_TIMING
      finish = std::chrono::high_resolution_clock::now();
      duration = finish - start;
      time2 += duration.count()*1.0e-6;
      start = std::chrono::high_resolution_clock::now();
#endif

      // Resort for next substep; this makes indx1
      //
      c->CudaSortByLevel();
    }

#ifdef VERBOSE_TIMING
    finish = std::chrono::high_resolution_clock::now();
    duration = finish - start;
    timeSRT += duration.count()*1.0e-6;
    start = std::chrono::high_resolution_clock::now();
#endif

    c->fix_positions_cuda();

#ifdef VERBOSE_TIMING
    finish = std::chrono::high_resolution_clock::now();
    duration = finish - start;
    timeCOM += duration.count()*1.0e-6;
#endif

    c->force->multistep_update_finish();
  }

#ifdef VERBOSE_TIMING
  auto finish0 = std::chrono::high_resolution_clock::now();
  duration = finish0 - start0;
  auto timeTOT = 1.0e-6*duration.count();

  std::cout << std::string(60, '-') << std::endl
	    << "Time in timestep  =" << std::setw(16) << time1
	    << std::setw(16) << time1/timeTOT   << std::endl
	    << "Time in timelevl  =" << std::setw(16) << time2
	    << std::setw(16) << time2/timeTOT << std::endl
	    << "Time in sort      =" << std::setw(16) << timeSRT
	    << std::setw(16) << timeSRT/timeTOT << std::endl
	    << "Time in adjust    =" << std::setw(16) << timeADJ
	    << std::setw(16) << timeADJ/timeTOT << std::endl
	    << "Time in COM       =" << std::setw(16) << timeCOM
	    << std::setw(16) << timeCOM/timeTOT << std::endl
	    << "Total time in adj =" << std::setw(16) << timeTOT
	    << std::setw(16) << 1.0 << std::endl
	    << std::string(60, '-') << std::endl;
#endif
}
