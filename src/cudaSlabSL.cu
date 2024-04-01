// -*- C++ -*-

#include <thrust/tuple.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#if CUDART_VERSION < 12000
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#endif

#include <cudaUtil.cuH>
#include <cudaReduce.cuH>
#include <Component.H>
#include <SlabSL.H>

#include "expand.H"

// Define for debugging
//
// #define BOUNDS_CHECK
// #define VERBOSE_RPT
// #define VERBOSE_DBG

// Global symbols for slab construction
//
__device__ __constant__
int slabNumX, slabNumY, slabNumZ, slabNX, slabNY, slabNZ, slabNdim, slabNum;


__device__ __constant__
int slabCmap;

__device__ __constant__
cuFP_t slabDfac, slabHscl, slabXmin, slabXmax, slabDxi;

// Alias for Thrust complex type to make this code more readable
//
using CmplxT = thrust::complex<cuFP_t>;

// Index functions for coefficients based on Eigen Tensor packing order
//
__device__
int slabIndex(int i, int j, int k)
{
  i += slabNumX;
  j += slabNumY;
  return k*slabNX*slabNY + j*slabNX + i;
}

// Index function for modulus coefficients
//
__device__
thrust::tuple<int, int, int> slabTensorIndices(int indx)
{
  int k = indx/(slabNX*slabNY);
  int j = (indx - k*slabNX*slabNY)/slabNX;
  int i = indx - (j + k*slabNY)*slabNX;

  return {i, j, k};
}

__device__
thrust::tuple<int, int, int> slabWaveNumbers(int indx)
{
  int k = indx/(slabNX*slabNY);
  int j = (indx - k*slabNX*slabNY)/slabNX;
  int i = indx - (j + k*slabNY)*slabNX;

  return {i-slabNumX, j-slabNumY, k};
}

__global__
void testConstantsSlab()
{
  printf("-------------------------\n");
  printf("---SlabBasis constants---\n");
  printf("-------------------------\n");
  printf("   Ndim   = %d\n", slabNdim );
  printf("   Numx   = %d\n", slabNumX );
  printf("   Numy   = %d\n", slabNumY );
  printf("   Numy   = %d\n", slabNumZ );
  printf("   Nx     = %d\n", slabNX   );
  printf("   Ny     = %d\n", slabNY   );
  printf("   Nz     = %d\n", slabNZ   );
  printf("   Dfac   = %e\n", slabDfac );
  printf("   Hscl   = %e\n", slabHscl );
  printf("   Cmap   = %d\n", slabCmap );
  printf("   Xmin   = %e\n", slabXmin );
  printf("   Xmax   = %e\n", slabXmax );
  printf("   Dxi    = %e\n", slabDxi  );
  printf("-------------------------\n");
}

__global__
void testFetchSlab(dArray<cudaTextureObject_t> T, dArray<cuFP_t> f,
		   int kx, int ky, int j, int nmax, int numz)
{
  const int n = blockDim.x * blockIdx.x + threadIdx.x;
  const int l = 1 + kx*(kx+1)/2*(slabNumY+1) + j;

#if cuREAL == 4
  if (n < numz) f._v[n] = tex1D<float>(T._v[l], n);
#else
  if (n < numz) f._v[n] = int2_as_double(tex1D<int2>(T._v[l], n));
#endif
}


thrust::host_vector<cuFP_t> returnTestSlab
(thrust::host_vector<cudaTextureObject_t>& tex,
 int kx, int ky, int j, int nmax, int numz)
{
  thrust::device_vector<cudaTextureObject_t> t_d = tex;
  
  unsigned int gridSize  = numz/BLOCK_SIZE;
  if (numz > gridSize*BLOCK_SIZE) gridSize++;
  
  thrust::device_vector<cuFP_t> f_d(numz);

  testFetchSlab<<<gridSize, BLOCK_SIZE>>>(toKernel(t_d), toKernel(f_d),
					  kx, ky, j, nmax, numz);

  cudaDeviceSynchronize();

  return f_d;
}

__device__
cuFP_t cu_z_to_xi(cuFP_t z)
{
  cuFP_t ret;

  if (slabCmap==0) {
    ret = tanh(z/slabHscl);
  } else if (slabCmap==1) {
    ret = z/sqrt(z*z + slabHscl*slabHscl);
  } else {
    ret = z;
  }    

  return ret;
}
    
__device__
cuFP_t cu_xi_to_z(cuFP_t xi)
{
  cuFP_t ret;

  if (slabCmap==0) {
    ret = slabHscl*atanh(xi);
  } else if (slabCmap==1) {
    ret = xi*slabHscl/sqrt(1.0 - xi*xi);
  } else {
    ret = xi;
  }

  return ret;
}

__device__
cuFP_t cu_d_xi_to_z(cuFP_t xi)
{
  cuFP_t ret;

  if (slabCmap==0) {
    ret = (1.0 - xi*xi)/slabHscl;
  } else if (slabCmap==1) {
    ret = pow(1.0 - xi*xi, 1.5)/slabHscl;
  } else {
    ret = 1.0;
  }

  return ret;
}

// Initialize anything cuda specific
//
void SlabSL::cuda_initialize()
{
  // Nothing
}

// Copy constants to device
//
void SlabSL::initialize_constants()
{
  auto f = grid->getCudaMappingConstants();
  cuFP_t z;

  cuda_safe_call(cudaMemcpyToSymbol(slabNumX, &nmaxx, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabNumX");

  cuda_safe_call(cudaMemcpyToSymbol(slabNumY, &nmaxy, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabNumY");

  cuda_safe_call(cudaMemcpyToSymbol(slabNumZ, &nmaxz, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabNumZ");

  cuda_safe_call(cudaMemcpyToSymbol(slabNX, &imx, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabNX");

  cuda_safe_call(cudaMemcpyToSymbol(slabNY, &imy, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabNY");

  cuda_safe_call(cudaMemcpyToSymbol(slabNZ, &imz, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabNZ");

  cuda_safe_call(cudaMemcpyToSymbol(slabNdim, &jmax, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabNdim");

  cuda_safe_call(cudaMemcpyToSymbol(slabNum, &f.numr, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabNum");

  int Cmap = 0;

  cuda_safe_call(cudaMemcpyToSymbol(slabCmap, &Cmap, sizeof(int),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabCmap");

  cuda_safe_call(cudaMemcpyToSymbol(slabDfac, &(z=2.0*M_PI), sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabDfac");

  cuda_safe_call(cudaMemcpyToSymbol(slabHscl, &(z=SLGridSlab::H), sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabHscl");

  cuda_safe_call(cudaMemcpyToSymbol(slabXmin, &(z=f.xmin), sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabXmin");

  cuda_safe_call(cudaMemcpyToSymbol(slabXmax, &(z=f.xmax), sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabXmax");

  cuda_safe_call(cudaMemcpyToSymbol(slabDxi, &(z=f.dxi), sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying slabDxi");
}


__global__ void coefKernelSlab
(dArray<cudaParticle> P, dArray<int> I, dArray<CmplxT> coef,
 dArray<cudaTextureObject_t> tex, int stride, PII lohi)
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


      // Restore particles to the unit slab
      //
      for (int k=0; k<2; k++) {
	if (pos[k]<0.0)
	  pos[k] += floor(-pos[k]) + 1.0;
	else
	  pos[k] += -floor(pos[k]);
      }
    
      // Wave number loop
      //
      const auto xx = CmplxT(0.0, slabDfac*pos[0]); // Phase values
      const auto yy = CmplxT(0.0, slabDfac*pos[1]);

      // Recursion increments and initial values
      //
      const auto sx = thrust::exp(-xx), cx = thrust::exp(xx*slabNumX);
      const auto sy = thrust::exp(-yy), cy = thrust::exp(yy*slabNumY);
      
      // Vertical interpolation
      //
      cuFP_t  x = cu_z_to_xi(pos[2]);
      cuFP_t xi = (x - slabXmin)/slabDxi;

      int ind   = floor(xi);
      int in0   = ind;

      if (in0 < 0) in0 = 0;
      if (in0 > slabNum-2) in0 = slabNum - 2;

      if (ind < 1) ind = 1;
      if (ind > slabNum-2) ind = slabNum - 2;

      cuFP_t  a = (cuFP_t)(in0+1) - xi;
      cuFP_t  b = 1.0 - a;

      // Flip sign for antisymmetric basis functions
      //
      int sign = 1;
      if (x<0 && 2*(n/2)!=n) sign = -1;


      cuFP_t p0 =
#if cuREAL == 4
        a*tex1D<float>(tex._v[0], ind  ) +
        b*tex1D<float>(tex._v[0], ind+1) ;
#else
        a*int2_as_double(tex1D<int2>(tex._v[0], ind  )) +
        b*int2_as_double(tex1D<int2>(tex._v[0], ind+1)) ;
#endif

      // Will contain the incremented basis
      //
      CmplxT X, Y;		
      
      X = cx;			// Assign the min X wavenumber conjugate

      for (int ii=-slabNumX; ii<=slabNumX; ii++, X*=sx) {

	Y = cy;			// Assign the min Y wavenumber conjugate

	int kx = abs(ii);

	for (int jj=-slabNumY; jj<=slabNumY; jj++, Y*=sy) {

	  int ky = abs(jj);

				// The vertical basis iteration
	  for (int n=0; n<slabNumZ; n++) {

	    int k = 1 + kx*(kx+1)/2*(slabNumY+1) + n;

#ifdef BOUNDS_CHECK
	    if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	    cuFP_t v = (
#if cuREAL == 4
			a*tex1D<float>(tex._v[k], ind  ) +
			b*tex1D<float>(tex._v[k], ind+1)
#else
			a*int2_as_double(tex1D<int2>(tex._v[k], ind  )) +
			b*int2_as_double(tex1D<int2>(tex._v[k], ind+1))
#endif
			) * p0 * sign;
	  
	    coef._v[slabIndex(ii, jj, n)] = -2.0*slabDfac * X * Y * v * mm;

	  }
	}
      }
      // END: wave number loop
    }
    // END: particle index limit
  }
  // END: stride loop
}


__global__ void
forceKernelSlab(dArray<cudaParticle> P, dArray<int> I,
		dArray<CmplxT> coef, dArray<cudaTextureObject_t> tex, int stride, PII lohi)
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
      
      CmplxT acc[3] = {0.0, 0.0, 0.0}, pot = 0.0, fac, facf;
      cuFP_t pos[3] = {p.pos[0], p.pos[1], p.pos[2]};
      cuFP_t mm = p.mass;

      // Wave number loop
      //
      const auto xx = CmplxT(0.0, slabDfac*pos[0]); // Phase values
      const auto yy = CmplxT(0.0, slabDfac*pos[1]);

      // Recursion increments and initial values
      const auto sx = thrust::exp(xx), cx = thrust::exp(-xx*slabNumX);
      const auto sy = thrust::exp(yy), cy = thrust::exp(-yy*slabNumY);

      // Vertical interpolation
      //
      cuFP_t  x = cu_z_to_xi(pos[2]);
      cuFP_t xi = (x - slabXmin)/slabDxi;
      cuFP_t dx = cu_d_xi_to_z(xi)/slabDxi;

      int ind   = floor(xi);
      int in0   = ind;

      if (in0 < 0) in0 = 0;
      if (in0 > slabNum-2) in0 = slabNum - 2;

      if (ind < 1) ind = 1;
      if (ind > slabNum-2) ind = slabNum - 2;

      cuFP_t  a = (cuFP_t)(in0+1) - xi;
      cuFP_t  b = 1.0 - a;
      

      // For 3-pt formula

      int jn0 = floor(xi);
      if (jn0 < 1) jn0 = 1;
      if (jn0 > slabNum-2) jn0 = slabNum - 2;

      cuFP_t s = (x - slabXmin - slabDxi*jn0)/slabDxi;

      cuFP_t p0 =
#if cuREAL == 4
        a*tex1D<float>(tex._v[0], ind  ) +
        b*tex1D<float>(tex._v[0], ind+1) ;
#else
        a*int2_as_double(tex1D<int2>(tex._v[0], ind  )) +
        b*int2_as_double(tex1D<int2>(tex._v[0], ind+1)) ;
#endif

      // Flip sign for antisymmetric basis functions
      int sign = 1;
      if (pos[2]<0 && 2*(n/2)!=n) sign = -1;

      CmplxT X, Y;		// Will contain the incremented basis

      X = cx;			// Assign the min X wavenumber
      for (int ii=-slabNumX; ii<=slabNumX; ii++, X*=sx) {
	Y = cy;			// Assign the min Y wavenumber
	for (int jj=-slabNumY; jj<=slabNumY; jj++, Y*=sy) {

	  for (int n=0; n<slabNumZ; n++) {

	      int kx = ii + slabNumX;
	      int ky = jj + slabNumY;
	      int k  = 1 + kx*(kx+1)/2*(slabNumY+1) + n;

#ifdef BOUNDS_CHECK
	      if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	      cuFP_t v = (
#if cuREAL == 4
			  a*tex1D<float>(tex._v[k], ind  ) +
			  b*tex1D<float>(tex._v[k], ind+1)
#else
			  a*int2_as_double(tex1D<int2>(tex._v[k], ind  )) +
			  b*int2_as_double(tex1D<int2>(tex._v[k], ind+1))
#endif
			  ) * p0 * sign;
	  
	      cuFP_t f = (
#if cuREAL == 4
			  (s - 0.5)*tex1D<float>(tex._v[k], jn0-1)*tex1D<float>(tex._v[0], jn0-1)
			  -2.0*tex1D<float>(tex._v[k], jnd)*tex1D<float>(tex._v[0], jn0) +
			  (s + 0.5)*tex1D<float>(tex._v[k], jn0+1)*tex1D<float>(tex._v[0], jn0+1)
#else
			  (s - 0.5)*int2_as_double(tex1D<int2>(tex._v[k], jn0-1))*
			  int2_as_double(tex1D<int2>(tex._v[0], jn0-1))
			  -2.0*int2_as_double(tex1D<int2>(tex._v[k], jn0))*
			  int2_as_double(tex1D<int2>(tex._v[0], jn0)) + 
			  (s + 0.5)*int2_as_double(tex1D<int2>(tex._v[k], jn0+1))*
			  int2_as_double(tex1D<int2>(tex._v[0], jn0+1))
#endif
			  ) * sign * dx;

	      fac  = X * Y * v * coef._v[slabIndex(ii, jj, n)];
	      facf = X * Y * f * coef._v[slabIndex(ii, jj, n)];

	      acc[0] += CmplxT(0.0, -slabDfac*ii) * fac;
	      acc[1] += CmplxT(0.0, -slabDfac*jj) * fac;
	      acc[2] += -facf;
	  }
	  // END: z wavenumber loop
	}
	// END: y wave number loop
      }
      // END: x wavenumber loop

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

void SlabSL::cudaStorage::resize_coefs(int N, int osize, int gridSize, int stride)
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


void SlabSL::cuda_zero_coefs()
{
  auto cr = component->cuStream;
  
  cuS.df_coef.resize(jmax);
    
  // Zero output array
  //
  thrust::fill(thrust::cuda::par.on(cr->stream),
	       cuS.df_coef.begin(), cuS.df_coef.end(), 0.0);
}

void SlabSL::determine_coefficients_cuda()
{
  // Only do this once but copying mapping coefficients and textures
  // must be done every time
  //
  if (initialize_cuda_slab) {
    initialize_cuda();
    initialize_cuda_slab = false;
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
  host_coefs.resize(jmax);

  // Get the stream for this component
  //
  auto cs = component->cuStream;

  // VERBOSE diagnostic output on first call
  //
  static bool firstime = true;

  if (firstime and myid==0 and VERBOSE>4) {
    testConstantsSlab<<<1, 1, 0, cs->stream>>>();
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
		<< "** cudaSlab coefficients" << std::endl
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
    
    // Adjust cached storage, if necessary
    //
    cuS.resize_coefs(N, jmax, gridSize, stride);
    
    // Compute the coefficient contribution for each order
    //
    auto beg  = cuS.df_coef.begin();
    
    coefKernelSlab<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
      (toKernel(cs->cuda_particles), toKernel(cs->indx1),
       toKernel(cuS.dN_coef), toKernel(t_d) ,stride, cur);
      
    // Begin the reduction by blocks [perhaps this should use a
    // stride?]
    //
    unsigned int gridSize1 = N/BLOCK_SIZE;
    if (N > gridSize1*BLOCK_SIZE) gridSize1++;

    reduceSum<CmplxT, BLOCK_SIZE>
      <<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
      (toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), jmax, N);
    
    // Finish the reduction for this order in parallel
    //
    thrust::counting_iterator<int> index_begin(0);
    thrust::counting_iterator<int> index_end(gridSize1*jmax);
    
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
    
    thrust::advance(beg, jmax);
  }

  // use1 += N;			// Increment particle count


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
	      << "---- Slab, T=" << tnow    << std::endl
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
    
    auto cmax = std::max_element(host_coefs.begin(), host_coefs.begin()+jmax,
				 LessAbs<CmplxT>());

    for (int n=0; n<jmax; n++) {
      auto [i, j, k] = indices(n);
      auto a = static_cast<std::complex<double>>(host_coefs[n]);
      auto b = expccof[0](i, j, k);
      auto c = std::abs(a - b);
      std::cout << std::setw(4)  << i-nmaxx
		<< std::setw(4)  << j-nmaxy
		<< std::setw(4)  << k
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
    if (myid==0) test.resize(2*nmaxx+1, 2*nmaxy+1, nmaxz);

    MPI_Reduce (thrust::raw_pointer_cast(&host_coefs[0]), test.data(), jmax,
		MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid==0) {

      std::string ofile = "slab_dump." + runtag + ".dat";
      std::ofstream out(ofile, ios::app | ios::out);

      if (out) {
	out << std::string(3*4+3*20, '-') << std::endl
	    << "---- Slab, T=" << tnow    << std::endl
	    << std::string(3*4+3*20, '-') << std::endl
	    << std::setprecision(10);
	
	out << std::setw(4)  << "i"
	    << std::setw(4)  << "j"
	    << std::setw(4)  << "k"
	    << std::setw(20) << "Real"
	    << std::setw(20) << "Imag"
	    << std::setw(20) << "Abs"
	    << std::endl;
	

	for (int n=0; n<jmax; n++) {
	  auto [i, j, k] = indices(n);
	  auto a = test(i, j, k);
	  out << std::setw(4)  << i-nmaxx
	      << std::setw(4)  << j-nmaxy
	      << std::setw(4)  << k
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

  // Deep debug for checking a single wave number from slabics
  //
  if (false) {

    coefType test;
    if (myid==0) test.resize(2*nmaxx+1, 2*nmaxy+1, nmaxz);

    MPI_Reduce (thrust::raw_pointer_cast(&host_coefs[0]), test.data(), jmax,
		MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid==0) {

      std::string ofile = "slab_test." + runtag + ".dat";
      std::ofstream out(ofile, ios::app | ios::out);

      if (out) {
	std::multimap<double, int> biggest;

	for (int n=0; n<jmax; n++)
	  biggest.insert({std::abs(test.data()[n]), n});

	out << std::string(3*4+3*20, '-') << std::endl
	    << "---- Slab, T=" << tnow    << std::endl
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
	      << std::setw(4)  << k
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

    std::string ofile = "test_slab." + runtag + ".dat";
    std::ofstream out(ofile);

    if (out) {

      // m loop
      for (int n=0; n<jmax; n++) {
	
	std::tie(elem.i, elem.j, elem.k) = indices(n);

	elem.d = expccof[0](elem.i, elem.j, elem.k);
	elem.f = static_cast<std::complex<double>>(host_coefs[n]);
	  
	double test = std::abs(elem.d - elem.f);
	if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);
	  
	compare.insert(std::make_pair(test, elem));
	
	out << std::setw( 5) << elem.i - nmaxx
	    << std::setw( 5) << elem.j - nmaxy
	    << std::setw( 5) << elem.k
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
		<< "---- Slab coefficients" << std::endl
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


void SlabSL::determine_acceleration_cuda()
{
  // Only do this once but copying mapping coefficients and textures
  // must be done every time
  //
  if (initialize_cuda_slab) {
    initialize_cuda();
    initialize_cuda_slab = false;
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
		<< "** cudaSlab acceleration" << std::endl
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
      
    forceKernelSlab<<<gridSize, BLOCK_SIZE, sMemSize, cs->stream>>>
      (toKernel(cs->cuda_particles), toKernel(cs->indx1),
       toKernel(dev_coefs), toKernel(t_d), stride, lohi);
  }
}

void SlabSL::HtoD_coefs()
{
  // Check size
  host_coefs.resize(jmax);

  // Copy from Slab
  for (int i=0; i<host_coefs.size(); i++)
    host_coefs[i] = expccof[0].data()[i];

  // Copy to device
  dev_coefs = host_coefs;
}


void SlabSL::DtoH_coefs(unsigned M)
{
  // Copy from host device to Slab
  for (int i=0; i<expccof[0].size(); i++)
    expccof[0].data()[i] = host_coefs[i];
}

void SlabSL::multistep_update_cuda()
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
		<< "] Adjust slab: Ntotal=" << Ntotal << " Npacks=" << Npacks
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

	// Adjust cached storage, if necessary
	//
	cuS.resize_coefs(N, jmax, gridSize, stride);
	
	// Compute the coefficient contribution for each order
	//
	auto beg  = cuS.df_coef.begin();
	
	// Do the work!
	//
	coefKernelSlab<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
	  (toKernel(cs->cuda_particles), toKernel(cs->indx1),
	   toKernel(cuS.dN_coef), toKernel(t_d), stride, cur);
	
	unsigned int gridSize1 = N/BLOCK_SIZE;
	if (N > gridSize1*BLOCK_SIZE) gridSize1++;
	
	reduceSum<CmplxT, BLOCK_SIZE>
	  <<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
	  (toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), jmax, N);
	
	// Finish the reduction for this order in parallel
	//
	thrust::counting_iterator<int> index_begin(0);
	thrust::counting_iterator<int> index_end(gridSize1*jmax);
	
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
	
	thrust::advance(beg, jmax);
      }
      // END: bunches

      // Accumulate the coefficients from the device to the host
      //
      thrust::host_vector<CmplxT> ret = cuS.df_coef;

      // Decrement current level and increment new level using the
      // Slab update matricies
      //
      for (int i=0; i<jmax; i++) {
	std::complex<double> val = ret[i];
	differ1[0][olev].data()[i] -= val;
	differ1[0][nlev].data()[i] += val;
      }
    }
    // DONE: Inner loop
  }
  // DONE: Outer loop
}


void SlabSL::destroy_cuda()
{
  // Nothing
}


