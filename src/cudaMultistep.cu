// -*- C++ -*-

#include <limits>

#include "expand.H"
#include <Component.H>
#include <cudaReduce.cuH>

// Define this for deep time step check
//
// #define OVER_UNDER

// Define this to see per operation sub timings
//
// #define VERBOSE_TIMING

// Global symbols for time step selection
//
__device__ __constant__
cuFP_t cuDynfracS, cuDynfracD, cuDynfracV, cuDynfracA, cuDynfracP, cuDtime;

__device__ __constant__
cuFP_t cuMaxDbl;

__device__ __constant__
int cuMultistep, cuShiftlev;

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

  cuda_safe_call(cudaMemcpyToSymbol(cuMaxDbl, &(z=1.0e32), sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuMaxDbl");
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
  printf("   MaxFL  = %e\n", cuMaxDbl   );
  printf("------------------------\n"   );
}

__global__ void
timestepSetKernel(dArray<cudaParticle> P, int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int npart = tid*stride + n;
    if (npart < P._s) P._v[npart].dtreq = cuMaxDbl;
  }
  // END: stride loop

}

__global__ void
timestepKernel(dArray<cudaParticle> P,
	       dArray<int> I,
	       dArray<int> loLev,
	       dArray<int> hiLev,
	       dArray<float> mindt,
	       dArray<float> maxdt,
	       cuFP_t cx, cuFP_t cy, cuFP_t cz,
	       int minactlev, int dim, int stride, PII lohi,
	       bool noswitch, bool apply)
{
  extern __shared__ float cache[];

  const int tid    = blockDim.x * blockIdx.x + threadIdx.x;
  const cuFP_t eps = 1.0e-10;
  const int cIndex = threadIdx.x;

  int lo = 0;
  int hi = 0;
  float minDT = 1.0e32;
  float maxDT = 0.0;

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
      
      cuFP_t dtd=1.0/eps, dtv=1.0/eps, dta=1.0/eps, dtA=1.0/eps, dts=1.0/eps;
      cuFP_t atot = 0.0, vtot = 0.0;

      // dtd = eps* rscale/v_i    -- char. drift time scale
      // dtv = eps* min(v_i/a_i)  -- char. force time scale
      // dta = eps* phi/(v * a)   -- char. work time scale
      // dtA = eps* sqrt(phi/a^2) -- char. "escape" time scale
      
      cuFP_t dtr  = 0.0;
      
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
      
      // Smallest time step
      //
      cuFP_t dt = dts;
      if (dt > dtd) dt = dtd;
      if (dt > dtv) dt = dtv;
      if (dt > dta) dt = dta;
      if (dt > dtA) dt = dtA;
      if (dt < eps) dt = eps;
      
      if (minDT > dt) minDT = dt;
      if (maxDT < dt) maxDT = dt;

      if (noswitch) {
	if (dt < p.dtreq) p.dtreq = dt;
      } else {
	p.dtreq = dt;
      }

      // Assign new levels?
      //
      if (apply) {

	// Time step wants to be LARGER than the maximum
	//
	if (p.dtreq >= cuDtime) {
	  p.lev[1] = 0;
	  hi++;
	} else
	  p.lev[1] = (int)floor(log(cuDtime/p.dtreq)/log(2.0));
	
	// Enforce n-level shifts at a time
	//
	if (cuShiftlev) {
	  if (p.lev[1] > p.lev[0]) {
	    if (p.lev[1] - p.lev[0] > cuShiftlev)
	      p.lev[1] = p.lev[0] + cuShiftlev;
	  } else if (p.lev[0] > p.lev[1]) {
	    if (p.lev[0] - p.lev[1] > cuShiftlev)
	      p.lev[1] = p.lev[0] - cuShiftlev;
	  }
	}

	// Time step wants to be SMALLER than the maximum
	//
	if (p.lev[1]>cuMultistep) {
	  p.lev[1] = cuMultistep;
	  lo++;
	}
	
	// Limit new level to minimum active level
	//
	if (p.lev[1]<minactlev) p.lev[1] = minactlev;
      }
      // Apply level selction
	
    } // Particle index block
    
  } // END: stride loop

  if (apply) {

    // set the cache values
    //
    cache[cIndex + 0*blockDim.x] = lo;
    cache[cIndex + 1*blockDim.x] = hi;
    cache[cIndex + 2*blockDim.x] = minDT;
    cache[cIndex + 3*blockDim.x] = maxDT;

    __syncthreads();
  
    // perform parallel reduction, threadsPerBlock must be 2^m
    //
    int ib = blockDim.x / 2;
    while (ib != 0) {
      
      if (cIndex < ib) {
	
	// Sums
	cache[cIndex + 0*blockDim.x] += cache[cIndex + ib + 0*blockDim.x]; 
	cache[cIndex + 1*blockDim.x] += cache[cIndex + ib + 1*blockDim.x]; 
	
	// Min
	if (cache[cIndex + ib + 2*blockDim.x] < cache[cIndex + 2*blockDim.x])
	  cache[cIndex + 2*blockDim.x] = cache[cIndex + ib + 2*blockDim.x]; 
	
	// Max
	if (cache[cIndex + ib + 3*blockDim.x] > cache[cIndex + 3*blockDim.x])
	  cache[cIndex + 3*blockDim.x] = cache[cIndex + ib + 3*blockDim.x]; 
      }

      __syncthreads();
      
      ib /= 2;
    }
    
    if (cIndex == 0) {
      loLev._v[blockIdx.x] = (int)cache[0];
      hiLev._v[blockIdx.x] = (int)cache[1];
      mindt._v[blockIdx.x] = cache[2];
      maxdt._v[blockIdx.x] = cache[3];
    }
  }
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
      
      p.lev[0] = p.lev[1];

    } // Particle index block
    
  } // END: stride loop

}


#ifdef OVER_UNDER

// Convert linear index to row index for column reduction
//
template <typename T>
struct to_row_index : public thrust::unary_function<T,T> {

  T Ncols; // --- Number of columns
  
  __host__ __device__ to_row_index(T Ncols) : Ncols(Ncols) {}
  
  __host__ __device__ T operator()(T i) { return i / Ncols; }
};


__global__ void
DTKernel(dArray<cudaParticle> P, dArray<int> I, dArray<int> ret,
	 int stride, PII lohi)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Index into the sorted array

    if (npart < lohi.second) {
      
      cudaParticle & p = P._v[I._v[npart]];
      
      if (p.dtreq >10.0) ret._v[i*2+0] += 1;
      if (p.dtreq<=10.0) ret._v[i*2+1] += 1;

    } // Particle index block
    
  } // END: stride loop

}

#endif

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

  bool firstCall = this_step==0 and mstep==0;

  //
  // Begin the update
  //
  for (auto c : comp->components) {
    c->force->multistep_update_begin();
    if (not c->force->cudaAware()) c->ParticlesToCuda();
  }

  std::map<Component*, int> offlo, offhi;
  float mindt = 1.0e32, maxdt = 0.0;

  cudaDeviceProp deviceProp;

#ifdef VERBOSE_TIMING
  double time1 = 0.0, time2 = 0.0, timeSRT = 0.0, timeADJ = 0.0, timeCOM = 0.0;
  auto start0 = std::chrono::high_resolution_clock::now();
  auto start  = std::chrono::high_resolution_clock::now();
#endif

  static bool firsttime = true;

  for (auto c : comp->components) {
    
    cudaGetDeviceProperties(&deviceProp, c->cudaDevice);
    cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

    // Initiate new level computation
    //
    bool apply = not c->NoSwitch() or mdrft==Mstep or firstCall;
    //                        ^            ^            ^
    //                        |            |            |
    // at every substep-------+            |            |
    //                                     |            |
    // otherwise: at end of full step------+            |
    //                                                  |
    // or on the very first call to initialize levels---+


    // Reset minimum time step field at the beginning of the master
    // step.  Only need to do this if the no-switch algorithm is on...
    //
    if (c->NoSwitch()) {

      // Zero dtreq?
      //
      if ((c->DTreset() and mstep==0) or firstCall) {

	// Compute grid
	//
	unsigned int N         = c->cuStream->cuda_particles.size();
	unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
	unsigned int gridSize  = N/BLOCK_SIZE/stride;
	
	if (N>0) {
	  if (N > gridSize*BLOCK_SIZE*stride) gridSize++;
	  
	  timestepSetKernel<<<gridSize, BLOCK_SIZE>>>
	    (toKernel(c->cuStream->cuda_particles), stride);
	}
      }
    }

    // Sort by levels
    //
    PII lohi = {0, c->cuStream->cuda_particles.size()};
    if (multistep) lohi = c->CudaGetLevelRange(mfirst[mdrft], multistep);
      
    // DEEP DEBUGGING
    if (false and multistep>0) {
      for (int n=0; n<numprocs; n++) {
	if (n==myid) testConstantsMultistep<<<1, 1>>>();
	std::cout << std::string(60, '-') << std::endl
		  << "[" << myid << ", " << c->name
		  << "]: T=" << tnow << " mlevel=" << mfirst[mdrft]
		  << " mstep=" << mstep << " mdrft=" << mdrft
		  << " (lo, hi) = (" << lohi.first << ", " << lohi.second << ")"
		  << std::endl << std::string(60, '-') << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    // END DEEP DEBUGGING

    // Compute grid
    //
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
    unsigned int sm        = BLOCK_SIZE*sizeof(float)*4;

    // We have particles at this level: get the new time-step level
    //
    if (N>0) {

      if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

      // Multilevel device storage
      //
      c->loLev.resize(gridSize);
      c->hiLev.resize(gridSize);
      c->minDT.resize(gridSize);
      c->maxDT.resize(gridSize);

      // Get the current expansion center
      //
      auto ctr = c->getCenter(Component::Local | Component::Centered);
      
      // Do the work
      //
      timestepKernel<<<gridSize, BLOCK_SIZE, sm>>>
	(toKernel(c->cuStream->cuda_particles),
	 toKernel(c->cuStream->indx1),
	 toKernel(c->loLev), toKernel(c->hiLev),
	 toKernel(c->minDT), toKernel(c->maxDT),
	 ctr[0], ctr[1], ctr[2],
	 mfirst[mdrft], c->dim, stride, lohi, c->NoSwitch(), apply);

      // Reductions
      //
      offlo[c] = thrust::reduce(c->loLev.begin(), c->loLev.end(),
				0, thrust::plus<int>());
      offhi[c] = thrust::reduce(c->hiLev.begin(), c->hiLev.end(),
				0, thrust::plus<int>());

      float minT = *thrust::min_element(c->minDT.begin(), c->maxDT.end());
      float maxT = *thrust::max_element(c->maxDT.begin(), c->maxDT.end());

      if (minT<mindt) mindt = minT;
      if (maxT>maxdt) maxdt = maxT;
    }
  }
  // END: component loop

  firsttime = false;

#ifdef VERBOSE_TIMING
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::micro> duration = finish - start;
  time1 += duration.count()*1.0e-6;
#endif
  
  // Finish by updating the coefficients and level list
  //
  for (auto c : comp->components) {

    bool firstCall = this_step==0 and mstep==0;
    bool apply = not c->NoSwitch() or mdrft==Mstep or firstCall;
    //                        ^             ^              ^
    //                        |             |              |
    // at every substep-------+             |              |
    //                                      |              |
    // otherwise: at end of full step-------+              |
    //                                                     |
    // or on the very first call to initialize levels------+
  
    // Freeze levels after first step
    //
    if (not firstCall and c->FreezeLev()) apply = false;

    // Call the multistep update for the force instance
    //
    if (apply) {
#ifdef VERBOSE_TIMING
      start = std::chrono::high_resolution_clock::now();
#endif
      cudaGetDeviceProperties(&deviceProp, c->cudaDevice);
      cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);
    
      // Update the coefficient level tableau
      //
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
    
      // Assign the newly computed, updated level
      //
      int dtover = 0, dtundr = 0;
      
      if (N>0) {

	if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

	timestepFinalizeKernel<<<gridSize, BLOCK_SIZE>>>
	  (toKernel(c->cuStream->cuda_particles),
	   toKernel(c->cuStream->indx1), stride, lohi);

	// Re-sort level array for next substep; this makes indx1
	//
	c->CudaSortByLevel();

#ifdef OVER_UNDER

	// Test uncomputed timesteps
	//
	const int Ncols = 2;
	thrust::device_vector<int> ret(Ncols*N);
	thrust::fill(ret.begin(), ret.end(), 0);

	// Allocate space for row sums and indices
	//
	thrust::device_vector<cuFP_t> d_sums   (Ncols);
	thrust::device_vector<int>    d_indices(Ncols);
	
	// Find big dtreq values
	// 
	DTKernel<<<gridSize, BLOCK_SIZE, 0, c->cuStream->stream>>>
	  (toKernel(c->cuStream->cuda_particles), toKernel(c->cuStream->indx1),
	   toKernel(ret), stride, lohi);
	
	// Perform sum over columns by summing values with equal column indices
	//
	thrust::reduce_by_key
	  (thrust::make_transform_iterator(thrust::counting_iterator<int>(0), to_row_index<int>(N)),
	   thrust::make_transform_iterator(thrust::counting_iterator<int>(0), to_row_index<int>(N)) + (N*Ncols),
	   thrust::make_permutation_iterator
	   (ret.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0),(thrust::placeholders::_1 % N) * Ncols + thrust::placeholders::_1 / N)),
	   d_indices.begin(),
	   d_sums.begin(),
	   thrust::equal_to<int>(),
	   thrust::plus<cuFP_t>());

	dtover += d_sums[0];
	dtundr += d_sums[1];
#endif

#ifdef VERBOSE_TIMING
	finish = std::chrono::high_resolution_clock::now();
	duration = finish - start;
	time2 += duration.count()*1.0e-6;
	start = std::chrono::high_resolution_clock::now();
#endif
      }

#ifdef OVER_UNDER
      // Sum errors
      //
      int totover = 0, totundr = 0;
      MPI_Reduce(&dtover, &totover, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&dtundr, &totundr, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      if (myid==0) {
	std::cout << "cudaMultistep: over=" << totover
		  << " under=" << totundr
		  << " T=" << tnow << std::endl;
      }
#endif

      // Commit the updates
      //
      c->force->multistep_update_finish();

    }
    // END: coefficient update stanza

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
  }
  // END: component loop

  // Check for multistep overrun
  //
  if (VERBOSE>0 && mdrft==Mstep) {

    unsigned sumlo=0, sumhi=0;
    for (auto c : comp->components) {
      sumlo += offlo[c];
      sumhi += offhi[c];
    }
      
    if (sumlo || sumhi) {
      std::cout << std::endl
		<< std::setw(70) << std::setfill('-') << '-'
		<< std::endl << std::setfill(' ')
		<< std::left << "--- Multistepping overrun" << std::endl;
      if (sumlo)
	std::cout << std::left << "--- Min DT=" << std::setw(16) << mindt  
		  << " < " << std::setw(16) << dtime/(1<<multistep) 
		  << " [" << sumlo << "]" <<  std::endl;
      if (sumhi)
	std::cout << std::left << "--- Max DT=" << std::setw(16) << maxdt  
		  << " > " << std::setw(16) << dtime 
		  << " [" << sumhi << "]" << std::endl;
      std::cout << std::setw(70) << std::setfill('-') << '-' << std::endl 
		<< std::setfill(' ') << std::right;
      
      if (sumlo) {
	for (auto c : comp->components) {
	  std::ostringstream sout;
	  sout << "Component <" << c->name << ">";
	  std::cout << std::setw(30) << sout.str() << " |   low: "
		    << offlo[c] << "/" << c->CurTotal() << std::endl;
	}
      }
      
      if (sumhi) {
	for (auto c : comp->components) {
	  std::ostringstream sout;
	  sout << "Component <" << c->name << ">";
	  std::cout << std::setw(30) << sout.str() << " |  high: "
		    << offhi[c] << "/" << c->CurTotal() << std::endl;
	}
      }
      
      std::cout << std::setw(70) << std::setfill('-') << '-'
		<< std::endl << std::setfill(' ');
    }
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
