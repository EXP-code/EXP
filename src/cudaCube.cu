// -*- C++ -*-

#include <thrust/tuple.h>

#include <cudaReduce.cuH>
#include <Component.H>
#include <Cube.H>

#include "expand.H"

// Define for debugging
//
// #define BOUNDS_CHECK
// #define VERBOSE_RPT
// #define VERBOSE_DBG
// #define NORECURSION

// Global symbols for cube construction
//
__device__ __constant__
int cubeNumX, cubeNumY, cubeNumZ, cubeNX, cubeNY, cubeNZ, cubeNdim;

__device__ __constant__
cuFP_t cubeDfac;

// Alias for Thrust complex type to make this code more readable
//
using CmplxT = thrust::complex<cuFP_t>;

// Index functions for coefficients based on Eigen Tensor packing order
//
__device__
int cubeIndex(int i, int j, int k)
{
  i += cubeNumX;
  j += cubeNumY;
  k += cubeNumZ;
  return k*cubeNX*cubeNY + j*cubeNX + i;
}

/*
// Index function for modulus coefficients
//
__device__
thrust::tuple<int, int, int> TensorIndices(int indx)
{
  int k = indx/(cubeNX*cubeNY);
  int j = (indx - k*cubeNX*cubeNY)/cubeNX;
  int i = indx - (j + k*cubeNY)*cubeNX;

  return {i, j, k};
}
*/

__device__
thrust::tuple<int, int, int> cubeWaveNumbers(int indx)
{
  int k = indx/(cubeNX*cubeNY);
  int j = (indx - k*cubeNX*cubeNY)/cubeNX;
  int i = indx - (j + k*cubeNY)*cubeNX;

  return {i-cubeNumX, j-cubeNumY, k-cubeNumZ};
}

__global__
void testConstantsCube()
{
  printf("-------------------------\n");
  printf("---CubeBasis constants---\n");
  printf("-------------------------\n");
  printf("   Ndim   = %d\n", cubeNdim );
  printf("   Numx   = %d\n", cubeNumX );
  printf("   Numy   = %d\n", cubeNumY );
  printf("   Numy   = %d\n", cubeNumZ );
  printf("   Nx     = %d\n", cubeNX   );
  printf("   Ny     = %d\n", cubeNY   );
  printf("   Nz     = %d\n", cubeNZ   );
  printf("   Dfac   = %e\n", cubeDfac );
  printf("-------------------------\n");
}

// Initialize anything cuda specific
//
void Cube::cuda_initialize()
{
  // Default
  //
  byPlanes = true;
  
  // Make method string lower case
  //
  std::transform(cuMethod.begin(), cuMethod.end(), cuMethod.begin(),
                 [](unsigned char c){ return std::tolower(c); });

  // Parse cuMethods variable
  //
  // All dimensions at once
  //
  if (cuMethod.find("all")    != std::string::npos) byPlanes = false;
  if (cuMethod.find("full")   != std::string::npos) byPlanes = false;
  if (cuMethod.find("3d")     != std::string::npos) byPlanes = false;
  //
  // Only one dimension at a time
  //
  if (cuMethod.find("planes") != std::string::npos) byPlanes = true;
  if (cuMethod.find("axes")   != std::string::npos) byPlanes = true;
  if (cuMethod.find("1d")     != std::string::npos) byPlanes = true;

  std::cout << "---- Cube::cuda_initialize: byPlanes="
	    << std::boolalpha << byPlanes << std::endl;
}

// Copy constants to device
//
void Cube::initialize_constants()
{
  cuda_safe_call(cudaMemcpyToSymbol(cubeNumX, &nmaxx, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cubeNumX");

  cuda_safe_call(cudaMemcpyToSymbol(cubeNumY, &nmaxy, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cubeNumY");

  cuda_safe_call(cudaMemcpyToSymbol(cubeNumZ, &nmaxz, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cubeNumZ");

  cuda_safe_call(cudaMemcpyToSymbol(cubeNX, &imx, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cubeNX");

  cuda_safe_call(cudaMemcpyToSymbol(cubeNY, &imy, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cubeNY");

  cuda_safe_call(cudaMemcpyToSymbol(cubeNZ, &imz, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cubeNZ");

  cuda_safe_call(cudaMemcpyToSymbol(cubeNdim, &osize, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cubeNdim");

  cuFP_t dfac = 2.0*M_PI;

  cuda_safe_call(cudaMemcpyToSymbol(cubeDfac, &dfac, sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cubeDfac");
}

__global__ void coefKernelCube
(dArray<cudaParticle> P, dArray<int> I, dArray<CmplxT> coef,
 int stride, PII lohi)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int N   = lohi.second - lohi.first;

  for (int n=0; n<stride; n++) {

    // Particle counter
    //
    int i     = tid*stride + n;
    int npart = i + lohi.first;

    if (npart < lohi.second) {	// Check that particle index is in
				// range for consistency with other
				// kernels

#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      
      cuFP_t pos[3] = {p.pos[0], p.pos[1], p.pos[2]};
      cuFP_t mm     = p.mass;

      bool good = true;
      for (int k=0; k<3; k++) {	// Are particles inside the unit cube?
	if (pos[k]<0.0 or pos[k]>1.0) good = false;
      }

      // Yes: do the expansion
      //
      if (good) {

#ifdef NORECURSION
	// Index loop
	//
	for (int s=0; s<cubeNdim; s++) {

	  // Get the wave numbers
	  int ii, jj, kk;
	  thrust::tie(ii, jj, kk) = cubeWaveNumbers(s);

	  // Skip the constant term and the divide by zero
	  //
	  if (ii!=0 or jj!=0 or kk!=0) {
	    cuFP_t expon = cubeDfac*(pos[0]*ii + pos[1]*jj + pos[2]*kk);
	    cuFP_t norm  = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));
	    
	    coef._v[s*N + i] = -mm*thrust::exp(CmplxT(0.0, -expon))*norm;
	  }
	}
	// END: index loop
#else
	// Wave number loop
	//
	const auto xx = CmplxT(0.0, cubeDfac*pos[0]); // Phase values
	const auto yy = CmplxT(0.0, cubeDfac*pos[1]);
	const auto zz = CmplxT(0.0, cubeDfac*pos[2]);

	// Recursion increments and initial values
	const auto sx = thrust::exp(-xx), cx = thrust::exp(xx*cubeNumX);
	const auto sy = thrust::exp(-yy), cy = thrust::exp(yy*cubeNumY);
	const auto sz = thrust::exp(-zz), cz = thrust::exp(zz*cubeNumZ);
	
	CmplxT X, Y, Z;		// Will contain the incremented basis
	
	X = cx;			// Assign the min X wavenumber conjugate
	for (int ii=-cubeNumX; ii<=cubeNumX; ii++, X*=sx) {
	  Y = cy;			// Assign the min Y wavenumber conjugate
	  for (int jj=-cubeNumY; jj<=cubeNumY; jj++, Y*=sy) {
	    Z = cz;		// Assign the min Z wavenumber conjugate
	    for (int kk=-cubeNumZ; kk<=cubeNumZ; kk++, Z*=sz) {
	      int l2 = ii*ii + jj*jj + kk*kk;
	      if (l2) {		// Only compute for non-zero l-index
		cuFP_t norm = -mm/sqrt(M_PI*l2);
		coef._v[cubeIndex(ii, jj, kk)*N + i] = X*Y*Z*norm;
	      }
	    }
	  }
	}
	// END: wave number loop
#endif
      }
      // END: particle in unit cube
    }
    // END: particle index limit
  }
  // END: stride loop

}

__global__ void coefKernelCubeX
(dArray<cudaParticle> P, dArray<int> I, dArray<CmplxT> coef,
 int jj, int kk, int stride, PII lohi)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int N   = lohi.second - lohi.first;

  for (int n=0; n<stride; n++) {

    // Particle counter
    //
    int i     = tid*stride + n;
    int npart = i + lohi.first;

    if (npart < lohi.second) {	// Check that particle index is in
				// range for consistency with other
				// kernels

#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      
      cuFP_t pos[3] = {p.pos[0], p.pos[1], p.pos[2]};
      cuFP_t mm     = p.mass;

      bool good = true;
      for (int k=0; k<3; k++) {	// Are particles inside the unit cube?
	if (pos[k]<0.0 or pos[k]>1.0) good = false;
      }

      // Yes: do the expansion
      //
      if (good) {

#ifdef NORECURSION
	// Index loop
	//
	for (int s=0; s<cubeNX; s++) {
	  
	  // Get the wave numbers
	  int kk = s - cubeNumX;

	  // Skip the constant term and the divide by zero
	  //
	  if (ii!=0 or jj!=0 or kk!=0) {
	    cuFP_t expon = cubeDfac*(pos[0]*ii + pos[1]*jj + pos[2]*kk);
	    cuFP_t norm  = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));
	    
	    coef._v[s*N + i] = -mm*thrust::exp(CmplxT(0.0, -expon))*norm;
	  }
	}
	// END: index loop
#else
	// Wave number loop
	//
	const auto xx = CmplxT(0.0, cubeDfac*pos[0]); // Phase values
	const auto yy = CmplxT(0.0, cubeDfac*pos[1]);
	const auto zz = CmplxT(0.0, cubeDfac*pos[2]);
	
	// Recursion increments and initial values
	const auto sx = thrust::exp(-xx), cx = thrust::exp(xx*cubeNumX);
	const auto Y = thrust::exp(-yy*jj), Z = thrust::exp(-zz*kk);
      
	CmplxT X;		// Will contain the incremented basis

	X = cx;			// Assign the min X wavenumber conjugate
	for (int ii=-cubeNumX; ii<=cubeNumX; ii++, X*=sx) {
	  int l2 = ii*ii + jj*jj + kk*kk;
	  if (l2) {		// Only compute for non-zero l-index
	    cuFP_t norm = -mm/sqrt(M_PI*l2);
	    coef._v[(ii+cubeNumX)*N + i] = X*Y*Z*norm;
	  }
	}
	// END: wave number loop
#endif
      }
      // END: particle in unit cube
    }
    // END: particle index limit
  }
  // END: stride loop

}


__global__ void
forceKernelCube(dArray<cudaParticle> P, dArray<int> I,
		dArray<CmplxT> coef, int stride, PII lohi)
{
  // Thread ID
  //
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {	// Check that particle index is in
				// range
      
#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      
      CmplxT acc[3] = {0.0, 0.0, 0.0}, pot = 0.0;
      cuFP_t pos[3] = {p.pos[0], p.pos[1], p.pos[2]};
      cuFP_t mm = p.mass;
      int ind[3];

#ifdef NORECURSION
      // Index loop
      //
      for (int s=0; s<cubeNdim; s++) {

	thrust::tie(ind[0], ind[1], ind[2]) = cubeWaveNumbers(s);

	// Exponent and normalization
	//
	cuFP_t expon = 0.0, norm = 0.0;
	for (int k=0; k<3; k++) {
	  expon += pos[k]*ind[k];
	  norm  += ind[k]*ind[k];
	}
	norm = 1.0/sqrt(M_PI*norm);

	// No force from the constant term
	//
	if (ind[0]!=0 or ind[1]!=0 or ind[2]!=0) {

	  auto pfac = coef._v[s] *
	    thrust::exp(CmplxT(0.0, cubeDfac*expon)) * norm;

	  pot += pfac;
	  for (int k=0; k<3; k++) {
	    acc[k] += CmplxT(0.0, -cubeDfac*ind[k]) * pfac;
	  }
	}
      }
      // END: index loop
#else
      // Wave number loop
      //
      const auto xx = CmplxT(0.0, cubeDfac*pos[0]); // Phase values
      const auto yy = CmplxT(0.0, cubeDfac*pos[1]);
      const auto zz = CmplxT(0.0, cubeDfac*pos[2]);

      // Recursion increments and initial values
      const auto sx = thrust::exp(xx), cx = thrust::exp(-xx*cubeNumX);
      const auto sy = thrust::exp(yy), cy = thrust::exp(-yy*cubeNumY);
      const auto sz = thrust::exp(zz), cz = thrust::exp(-zz*cubeNumZ);

      CmplxT X, Y, Z;		// Will contain the incremented basis

      X = cx;			// Assign the min X wavenumber
      for (int ii=-cubeNumX; ii<=cubeNumX; ii++, X*=sx) {
	Y = cy;			// Assign the min Y wavenumber
	for (int jj=-cubeNumY; jj<=cubeNumY; jj++, Y*=sy) {
	  Z = cz;		// Assign the min Z wavenumber
	  for (int kk=-cubeNumZ; kk<=cubeNumZ; kk++, Z*=sz) {
	    int l2 = ii*ii + jj*jj + kk*kk;
	    if (l2) {		// Only compute the non-constant terms
	      cuFP_t norm  = 1.0/sqrt(M_PI*l2);
	      auto pfac = coef._v[cubeIndex(ii, jj, kk)] * X*Y*Z*norm;

	      pot    += pfac;	// Potential and force vector
	      acc[0] += CmplxT(0.0, -cubeDfac*ii) * pfac;
	      acc[1] += CmplxT(0.0, -cubeDfac*jj) * pfac;
	      acc[2] += CmplxT(0.0, -cubeDfac*kk) * pfac;
	    }
	  }
	  // END: z wavenumber loop
	}
	// END: y wave number loop
      }
      // END: x wavenumber loop
#endif
      // Particle assignment
      //
      p.pot = pot.real();
      for (int k=0; k<3; k++) p.acc[k] = acc[k].real();
    }
    // END: particle index limit
  }
  // END: stride loop
}

template<typename T>
class LessAbs : public std::binary_function<bool, T, T>
{
public:
  bool operator()( const T &a, const T &b ) const
  {
    return (thrust::abs(a) < thrust::abs(b));
  }
};

void Cube::cudaStorage::resize_coefs(int N, int osize, int gridSize, int stride)
{
  // Reserve space for coefficient reduction
  //
  if (dN_coef.capacity() < osize*N)
    dN_coef.reserve(osize*N);
  
  if (dc_coef.capacity() < osize*gridSize)
    dc_coef.reserve(osize*gridSize);
  
  // Set space for current step
  //
  dN_coef.resize(osize*N);
  dc_coef.resize(osize*gridSize);
  dw_coef.resize(osize);	// This will stay fixed
}


void Cube::cuda_zero_coefs()
{
  auto cr = component->cuStream;
  
  cuS.df_coef.resize(osize);
    
  // Zero output array
  //
  thrust::fill(thrust::cuda::par.on(cr->stream),
	       cuS.df_coef.begin(), cuS.df_coef.end(), 0.0);
}

void Cube::determine_coefficients_cuda()
{
  // Only do this once but copying mapping coefficients and textures
  // must be done every time
  //
  if (initialize_cuda_cube) {
    initialize_cuda();
    initialize_cuda_cube = false;
  }

  // Copy coordinate mapping
  //
  initialize_constants();


  std::cout << std::scientific;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, component->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  // This will stay fixed for the entire run
  //
  host_coefs.resize(osize);

  // Get the stream for this component
  //
  auto cs = component->cuStream;

  // VERBOSE diagnostic output on first call
  //
  static bool firstime = true;

  if (firstime and myid==0 and VERBOSE>4) {
    testConstantsCube<<<1, 1, 0, cs->stream>>>();
    cudaDeviceSynchronize();
    cuda_check_last_error_mpi("cudaDeviceSynchronize", __FILE__, __LINE__, myid);
    firstime = false;
  }
  
  // Zero counter and coefficients
  //
  thrust::fill(host_coefs.begin(), host_coefs.end(), 0.0);

  // Zero out coefficient storage
  //
  cuda_zero_coefs();

  // Get sorted particle range for mlevel
  //
  PII lohi = component->CudaGetLevelRange(mlevel, mlevel), cur;

  if (false) {
    for (int n=0; n<numprocs; n++) {
      if (myid==n) std::cout << "[" << myid << "] mlevel=" << mlevel
			     << " coef check (lo, hi) = (" << lohi.first << ", "
			     << lohi.second << ")" << std::endl
			     << std::string(60, '-') << std::endl;
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  
  unsigned int Ntotal = lohi.second - lohi.first;
  unsigned int Npacks = Ntotal/component->bunchSize + 1;

  // Loop over bunches
  //
  for (int n=0; n<Npacks; n++) {

    // Current bunch
    //
    cur. first = lohi.first + component->bunchSize*n;
    cur.second = lohi.first + component->bunchSize*(n+1);
    cur.second = std::min<unsigned int>(cur.second, lohi.second);

    if (cur.second <= cur.first) break;
    
    // Compute grid
    //
    unsigned int N         = cur.second - cur.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
  
    if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

#ifdef VERBOSE_RPT
    static unsigned debug_max_count = 100;
    static unsigned debug_cur_count = 0;
    if (debug_cur_count++ < debug_max_count) {
      std::cout << std::endl
		<< "** -------------------------" << std::endl
		<< "** cudaCube coefficients" << std::endl
		<< "** -------------------------" << std::endl
		<< "** N      = " << N            << std::endl
		<< "** Npacks = " << Npacks       << std::endl
		<< "** I low  = " << cur.first    << std::endl
		<< "** I high = " << cur.second   << std::endl
		<< "** Stride = " << stride       << std::endl
		<< "** Block  = " << BLOCK_SIZE   << std::endl
		<< "** Grid   = " << gridSize     << std::endl
		<< "** Level  = " << mlevel       << std::endl
		<< "** lo     = " << lohi.first   << std::endl
		<< "** hi     = " << lohi.second  << std::endl
		<< "**" << std::endl;
  }
#endif
  
    // Shared memory size for the reduction
    //
    int sMemSize = BLOCK_SIZE * sizeof(CmplxT);
    
    if (byPlanes) {

      // Adjust cached storage, if necessary
      //
      cuS.resize_coefs(N, imx, gridSize, stride);
    
      // Compute the coefficient contribution for each order
      //
      auto beg  = cuS.df_coef.begin();

      for (int kk=-nmaxz; kk<=nmaxz; kk++) {

	for (int jj=-nmaxy; jj<=nmaxy; jj++) {
	
	  coefKernelCubeX<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
	    (toKernel(cs->cuda_particles), toKernel(cs->indx1),
	     toKernel(cuS.dN_coef), jj, kk, stride, cur);
      
	  // Begin the reduction by blocks [perhaps this should use a
	  // stride?]
	  //
	  unsigned int gridSize1 = N/BLOCK_SIZE;
	  if (N > gridSize1*BLOCK_SIZE) gridSize1++;

	  reduceSum<CmplxT, BLOCK_SIZE>
	    <<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
	    (toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), imx, N);
      
	  // Finish the reduction for this order in parallel
	  //
	  thrust::counting_iterator<int> index_begin(0);
	  thrust::counting_iterator<int> index_end(gridSize1*imx);
      
	  // The key_functor indexes the sum reduced series by array index
	  //
	  thrust::reduce_by_key
	    (
	     thrust::cuda::par.on(cs->stream),
	     thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
	     thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
	     cuS.dc_coef.begin(), thrust::make_discard_iterator(), cuS.dw_coef.begin()
	     );
      
	  thrust::transform(thrust::cuda::par.on(cs->stream),
			    cuS.dw_coef.begin(), cuS.dw_coef.end(),
			    beg, beg, thrust::plus<CmplxT>());
	  
	  thrust::advance(beg, imx);
	}
	// END: y wave number loop
      }
      // END: z wave number loop
      
    } else {

      // Adjust cached storage, if necessary
      //
      cuS.resize_coefs(N, osize, gridSize, stride);
    
      // Compute the coefficient contribution for each order
      //
      auto beg  = cuS.df_coef.begin();

      coefKernelCube<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
	(toKernel(cs->cuda_particles), toKernel(cs->indx1),
	 toKernel(cuS.dN_coef), stride, cur);
      
      // Begin the reduction by blocks [perhaps this should use a
      // stride?]
      //
      unsigned int gridSize1 = N/BLOCK_SIZE;
      if (N > gridSize1*BLOCK_SIZE) gridSize1++;

      reduceSum<CmplxT, BLOCK_SIZE>
	<<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
	(toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), osize, N);
      
      // Finish the reduction for this order in parallel
      //
      thrust::counting_iterator<int> index_begin(0);
      thrust::counting_iterator<int> index_end(gridSize1*osize);
      
      // The key_functor indexes the sum reduced series by array index
      //
      thrust::reduce_by_key
	(
	 thrust::cuda::par.on(cs->stream),
	 thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
	 thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
	 cuS.dc_coef.begin(), thrust::make_discard_iterator(), cuS.dw_coef.begin()
       );
      
      thrust::transform(thrust::cuda::par.on(cs->stream),
			cuS.dw_coef.begin(), cuS.dw_coef.end(),
			beg, beg, thrust::plus<CmplxT>());
      
      thrust::advance(beg, osize);
    }

    use1 += N;			// Increment particle count
  }

  // Accumulate the coefficients from the device to the host
  //
  host_coefs = cuS.df_coef;

  // Create a wavenumber tuple from a flattened index
  //
  auto indices = [&](int indx)
  {
    int NX = 2*this->nmaxx+1, NY = 2*this->nmaxy+1;
    int k  = indx/(NX*NY);
    int j  = (indx - k*NX*NY)/NX;
    int i  = indx - (j + k*NY)*NX;
    
    return std::tuple<int, int, int>{i, j, k};
  };

  // DEBUG
  //
  if (false) {
    std::cout << std::string(3*4+4*20, '-') << std::endl
	      << "---- Cube, T=" << tnow    << std::endl
	      << std::string(3*4+4*20, '-') << std::endl
	      << std::setprecision(10);

    std::cout << std::setw(4)  << "i"
	      << std::setw(4)  << "j"
	      << std::setw(4)  << "k"
	      << std::setw(20) << "GPU"
	      << std::setw(20) << "CPU"
	      << std::setw(20) << "diff"
	      << std::setw(20) << "rel diff"
	      << std::endl;
    
    auto cmax = std::max_element(host_coefs.begin(), host_coefs.begin()+osize,
				 LessAbs<CmplxT>());

    for (int n=0; n<osize; n++) {
      auto [i, j, k] = indices(n);
      auto a = static_cast<std::complex<double>>(host_coefs[n]);
      auto b = expcoef[0](i, j, k);
      auto c = std::abs(a - b);
      std::cout << std::setw(4)  << i-nmaxx
		<< std::setw(4)  << j-nmaxy
		<< std::setw(4)  << k-nmaxz
		<< std::setw(20) << a
		<< std::setw(20) << b
		<< std::setw(20) << c
		<< std::setw(20) << c/thrust::abs(*cmax)
		<< std::endl;
    }

    std::cout << std::string(3*4+4*20, '-') << std::endl;
  }


  // Dump the coefficients for sanity checking (set to false for
  // production)
  //
  if (false) {

    coefType test;
    if (myid==0) test.resize(2*nmaxx+1, 2*nmaxy+1, 2*nmaxz+1);

    MPI_Reduce (thrust::raw_pointer_cast(&host_coefs[0]), test.data(), osize,
		MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid==0) {

      std::string ofile = "cube_dump." + runtag + ".dat";
      std::ofstream out(ofile, ios::app | ios::out);

      if (out) {
	out << std::string(3*4+3*20, '-') << std::endl
	    << "---- Cube, T=" << tnow    << std::endl
	    << std::string(3*4+3*20, '-') << std::endl
	    << std::setprecision(10);
	
	out << std::setw(4)  << "i"
	    << std::setw(4)  << "j"
	    << std::setw(4)  << "k"
	    << std::setw(20) << "Real"
	    << std::setw(20) << "Imag"
	    << std::setw(20) << "Abs"
	    << std::endl;
	

	for (int n=0; n<osize; n++) {
	  auto [i, j, k] = indices(n);
	  auto a = test(i, j, k);
	  out << std::setw(4)  << i-nmaxx
	      << std::setw(4)  << j-nmaxy
	      << std::setw(4)  << k-nmaxz
	      << std::setw(20) << std::real(a)
	      << std::setw(20) << std::imag(a)
	      << std::setw(20) << std::abs(a)
	      << std::endl;
	}
	out << std::string(3*4+4*20, '-') << std::endl;
      } else {
	std::cout << "Error opening <" << ofile << ">" << std::endl;
      }
    }
  }

  // Deep debug for checking a single wave number from cubeics
  //
  if (false) {

    coefType test;
    if (myid==0) test.resize(2*nmaxx+1, 2*nmaxy+1, 2*nmaxz+1);

    MPI_Reduce (thrust::raw_pointer_cast(&host_coefs[0]), test.data(), osize,
		MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid==0) {

      std::string ofile = "cube_test." + runtag + ".dat";
      std::ofstream out(ofile, ios::app | ios::out);

      if (out) {
	std::multimap<double, int> biggest;

	for (int n=0; n<osize; n++)
	  biggest.insert({std::abs(test.data()[n]), n});

	out << std::string(3*4+3*20, '-') << std::endl
	    << "---- Cube, T=" << tnow    << std::endl
	    << std::string(3*4+3*20, '-') << std::endl
	    << std::setprecision(10);
	
	out << std::setw(4)  << "i"
	    << std::setw(4)  << "j"
	    << std::setw(4)  << "k"
	    << std::setw(20) << "Real"
	    << std::setw(20) << "Imag"
	    << std::setw(20) << "Abs"
	    << std::endl;
	
	int cnt = 0;
	for (auto it = biggest.rbegin(); it!=biggest.rend() and cnt<20; it++, cnt++) {
	  auto [i, j, k] = indices(it->second);
	  auto a = test(i, j, k);

	  out << std::setw(4)  << i-nmaxx
	      << std::setw(4)  << j-nmaxy
	      << std::setw(4)  << k-nmaxz
	      << std::setw(20) << std::real(a)
	      << std::setw(20) << std::imag(a)
	      << std::setw(20) << std::abs(a)
	      << std::endl;
	}
	out << std::string(3*4+4*20, '-') << std::endl;
      } else {
	std::cout << "Error opening <" << ofile << ">" << std::endl;
      }
    }
  }


  //
  // TEST comparison of coefficients for debugging
  //
  if (false and myid==0) {

    struct Element
    {
      std::complex<double> d;
      std::complex<double> f;
      
      int    i;
      int    j;
      int    k;
    }
    elem;

    std::multimap<double, Element> compare;

    std::string ofile = "test_cube." + runtag + ".dat";
    std::ofstream out(ofile);

    if (out) {

      // m loop
      for (int n=0; n<osize; n++) {
	
	std::tie(elem.i, elem.j, elem.k) = indices(n);

	elem.d = expcoef[0](elem.i, elem.j, elem.k);
	elem.f = static_cast<std::complex<double>>(host_coefs[n]);
	  
	double test = std::abs(elem.d - elem.f);
	if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);
	  
	compare.insert(std::make_pair(test, elem));
	
	out << std::setw( 5) << elem.i - nmaxx
	    << std::setw( 5) << elem.j - nmaxy
	    << std::setw( 5) << elem.k - nmaxz
	    << std::setw( 5) << n
	    << std::setw(20) << elem.d
	    << std::setw(20) << elem.f
	    << std::endl;
      }
    
      std::map<double, Element>::iterator best = compare.begin();
      std::map<double, Element>::iterator midl = best;
      std::advance(midl, compare.size()/2);
      std::map<double, Element>::reverse_iterator last = compare.rbegin();
      
      std::cout << std::string(3*3 + 3*20 + 20, '-') << std::endl
		<< "---- Cube coefficients" << std::endl
		<< std::string(3*3 + 3*20 + 20, '-') << std::endl;
      
      std::cout << "Best case: ["
		<< std::setw( 3) << best->second.i << ", "
		<< std::setw( 3) << best->second.j << ", "
		<< std::setw( 3) << best->second.k << "] = "
		<< std::setw(20) << best->second.d
		<< std::setw(20) << best->second.f
		<< std::setw(20) << fabs(best->second.d - best->second.f)
		<< std::endl;
  
      std::cout << "Mid case:  ["
		<< std::setw( 3) << midl->second.i << ", "
		<< std::setw( 3) << midl->second.j << ", "
		<< std::setw( 3) << midl->second.k << "] = "
		<< std::setw(20) << midl->second.d
		<< std::setw(20) << midl->second.f
		<< std::setw(20) << fabs(midl->second.d - midl->second.f)
		<< std::endl;
    
      std::cout << "Last case: ["
		<< std::setw( 3) << last->second.i << ", "
		<< std::setw( 3) << last->second.j << ", "
		<< std::setw( 3) << last->second.k << "] = "
		<< std::setw(20) << last->second.d
		<< std::setw(20) << last->second.f
		<< std::setw(20) << fabs(last->second.d - last->second.f)
		<< std::endl;
    }
  }

}


void Cube::determine_acceleration_cuda()
{
  // Only do this once but copying mapping coefficients and textures
  // must be done every time
  //
  if (initialize_cuda_cube) {
    initialize_cuda();
    initialize_cuda_cube = false;
  }

  // Copy coordinate mapping
  //
  initialize_constants();

  std::cout << std::scientific;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cC->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  auto cs = cC->cuStream;

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

#ifdef VERBOSE_RPT
    static unsigned debug_max_count = 100;
    static unsigned debug_cur_count = 0;
    if (debug_cur_count++ < debug_max_count) {
      std::cout << std::endl
		<< "** -------------------------" << std::endl
		<< "** cudaCube acceleration" << std::endl
		<< "** -------------------------" << std::endl
		<< "** N      = " << N            << std::endl
		<< "** Stride = " << stride       << std::endl
		<< "** Block  = " << BLOCK_SIZE   << std::endl
		<< "** Grid   = " << gridSize     << std::endl
		<< "** Level  = " << mlevel       << std::endl
		<< "** lo     = " << lohi.first   << std::endl
		<< "** hi     = " << lohi.second  << std::endl
		<< "**" << std::endl;
    }
#endif
    
    // Shared memory size for the reduction
    //
    int sMemSize = BLOCK_SIZE * sizeof(CmplxT);
      
    forceKernelCube<<<gridSize, BLOCK_SIZE, sMemSize, cs->stream>>>
      (toKernel(cs->cuda_particles), toKernel(cs->indx1),
       toKernel(dev_coefs), stride, lohi);
  }
}

void Cube::HtoD_coefs()
{
  // Check size
  host_coefs.resize(osize);

  // Copy from Cube
  for (int i=0; i<host_coefs.size(); i++)
    host_coefs[i] = expcoef[0].data()[i];

  // Copy to device
  dev_coefs = host_coefs;
}


void Cube::DtoH_coefs(unsigned M)
{
  // Copy from host device to Cube
  for (int i=0; i<expcoef[0].size(); i++)
    expcoef[0].data()[i] = host_coefs[i];
}

void Cube::multistep_update_cuda()
{
  // The plan: for the current active level search above and below for
  // particles for correction to coefficient matrix
  //

  //! Sort the device vector by level changes
  auto chg = component->CudaSortLevelChanges();

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, component->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);
  auto cs = component->cuStream;
  
  // Step through all levels
  //
  for (int olev=mfirst[mstep]; olev<=multistep; olev++) {

    for (int nlev=0; nlev<=multistep; nlev++) {

      if (olev == nlev) continue;

      // Get range of update block in particle index
      //
      unsigned int Ntotal = chg[olev][nlev].second - chg[olev][nlev].first;

      if (Ntotal==0) continue; // No particles [from, to]=[olev, nlev]

      unsigned int Npacks = Ntotal/component->bunchSize + 1;

      // Zero out coefficient storage
      //
      cuda_zero_coefs();

#ifdef VERBOSE_DBG
      std::cout << "[" << myid << ", " << tnow
		<< "] Adjust cube: Ntotal=" << Ntotal << " Npacks=" << Npacks
		<< " for (m, d)=(" << olev << ", " << nlev << ")" << std::endl;
#endif

      // Loop over bunches
      //
      for (int n=0; n<Npacks; n++) {

	PII cur;
	
	// Current bunch
	//
	cur. first = chg[olev][nlev].first + component->bunchSize*n;
	cur.second = chg[olev][nlev].first + component->bunchSize*(n+1);
	cur.second = std::min<unsigned int>(cur.second, chg[olev][nlev].second);

	if (cur.second <= cur.first) break;
    
	// Compute grid
	//
	unsigned int N         = cur.second - cur.first;
	unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
	unsigned int gridSize  = N/BLOCK_SIZE/stride;
	
	if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

	// Shared memory size for the reduction
	//
	int sMemSize = BLOCK_SIZE * sizeof(CmplxT);

	if (byPlanes) {
	  // Adjust cached storage, if necessary
	  //
	  cuS.resize_coefs(N, imx, gridSize, stride);
	
	  // Compute the coefficient contribution for each order
	  //
	  auto beg  = cuS.df_coef.begin();

	  for (int kk=-nmaxz; kk<=nmaxz; kk++) {
	    
	    for (int jj=-nmaxy; jj<=nmaxy; jj++) {

	      // Do the work!
	      //
	      coefKernelCubeX<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
		(toKernel(cs->cuda_particles), toKernel(cs->indx1),
		 toKernel(cuS.dN_coef), jj, kk, stride, cur);
      
	      unsigned int gridSize1 = N/BLOCK_SIZE;
	      if (N > gridSize1*BLOCK_SIZE) gridSize1++;
	  
	      reduceSum<CmplxT, BLOCK_SIZE>
		<<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
		(toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), imx, N);
	  
	      // Finish the reduction for this order in parallel
	      //
	      thrust::counting_iterator<int> index_begin(0);
	      thrust::counting_iterator<int> index_end(gridSize1*imx);
	      
	      // The key_functor indexes the sum reduced series by array index
	      //
	      thrust::reduce_by_key
		(
		 thrust::cuda::par.on(cs->stream),
		 thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
		 thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
		 cuS.dc_coef.begin(), thrust::make_discard_iterator(), cuS.dw_coef.begin()
		 );

	      thrust::transform(thrust::cuda::par.on(cs->stream),
				cuS.dw_coef.begin(), cuS.dw_coef.end(),
				beg, beg, thrust::plus<CmplxT>());
	  
	      thrust::advance(beg, imx);
	    }
	    // END: z wave numbers
	  }
	  // END: y wave numbers
	}
	// END: bunch by planes
	else {
	  // Adjust cached storage, if necessary
	  //
	  cuS.resize_coefs(N, osize, gridSize, stride);
	
	  // Shared memory size for the reduction
	  //
	  int sMemSize = BLOCK_SIZE * sizeof(CmplxT);
    
	  // Compute the coefficient contribution for each order
	  //
	  auto beg  = cuS.df_coef.begin();

	  // Do the work!
	  //
	  coefKernelCube<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
	    (toKernel(cs->cuda_particles), toKernel(cs->indx1),
	     toKernel(cuS.dN_coef), stride, cur);
      
	  unsigned int gridSize1 = N/BLOCK_SIZE;
	  if (N > gridSize1*BLOCK_SIZE) gridSize1++;
	  
	  reduceSum<CmplxT, BLOCK_SIZE>
	    <<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
	    (toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), osize, N);
	  
	  // Finish the reduction for this order in parallel
	  //
	  thrust::counting_iterator<int> index_begin(0);
	  thrust::counting_iterator<int> index_end(gridSize1*osize);
	  
	  // The key_functor indexes the sum reduced series by array index
	  //
	  thrust::reduce_by_key
	    (
	     thrust::cuda::par.on(cs->stream),
	     thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
	     thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
	     cuS.dc_coef.begin(), thrust::make_discard_iterator(), cuS.dw_coef.begin()
	     );

	  thrust::transform(thrust::cuda::par.on(cs->stream),
			    cuS.dw_coef.begin(), cuS.dw_coef.end(),
			    beg, beg, thrust::plus<CmplxT>());
	  
	  thrust::advance(beg, osize);
	}
	// END: all wave numbers at once
      }
      // END: bunches

      // Accumulate the coefficients from the device to the host
      //
      thrust::host_vector<CmplxT> ret = cuS.df_coef;

      // Decrement current level and increment new level using the
      // Cube update matricies
      //
      for (int i=0; i<osize; i++) {
	std::complex<double> val = ret[i];
	differ1[0][olev].data()[i] -= val;
	differ1[0][nlev].data()[i] += val;
      }
    }
    // DONE: Inner loop
  }
  // DONE: Outer loop
}


void Cube::destroy_cuda()
{
  // Nothing
}


