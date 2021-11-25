// -*- C++ -*-

#ifndef _CUDAHELP_H
#define _CUDAHELP_H

#define BLOCK_SIZE 256

#include <sstream>
#include <string>

#include <mpi.h>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <thrust/system/cuda/error.h>
#include <thrust/system_error.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>

//! Wrapper function for a cuda call that processes error
inline void
cuda_safe_call(cudaError_t code, const char *file, int line,
	       const std::string& message = "")
{
  if (code != cudaSuccess) {
    const char * errorMessage = cudaGetErrorString(code);
    std::stringstream ss;
    ss << "CUDA error returned from " << file << ":" << line
       << " Error code: " << code << "(" << errorMessage << ")";
    if (message.size()>0) ss << " [" << message << "]";
    
    throw thrust::system_error(code, thrust::cuda_category(), ss.str());
  }
}

//! Wrapper function for a cuda call that processes error [MPI version]
inline void
cuda_safe_call_mpi(cudaError_t code, const char *file, int line, int myid,
		   const std::string& message = "")
{
  std::stringstream ss;
  int bad = 0, tot = 0;

  if (code != cudaSuccess) {
    const char * errorMessage = cudaGetErrorString(code);
    ss << "CUDA error returned from " << file << ":" << line
       << ", mpi rank=" << myid
       << " Error code: " << code << "(" << errorMessage << ")";
    if (message.size()>0) ss << " [" << message << "]";
    bad = 1;
  }
  
  // All processes get the total number of errors
  MPI_Allreduce(&bad, &tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Print local error
  if (bad) {
    std::cerr << "[" << myid << "] " << ss.str() << std::endl << std::flush;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // All processes throw if any process has an error
  if (tot) {
    throw thrust::system_error(code, thrust::cuda_category(), ss.str());
  }
}

// For error reporting with error code
//
inline void
cuda_check_error(cudaError_t code,
		 char const * func, const char *FILE, const int LINE)
{
  if (code != cudaSuccess) {
    const char * errorMessage = cudaGetErrorString(code);
    std::cerr << "CUDA error returned from \"" << func << "\" at "
	      << FILE << ":" << LINE
	      << " Error code: " << code << "(" << errorMessage << ")"
	      << std::endl;
  }
}
 
// For error reporting without error code
//
inline void
cuda_check_last_error(char const * func, const char *FILE, const int LINE)

{
  cudaError_t code = cudaGetLastError();
  if (code != cudaSuccess) {
    const char * errorMessage = cudaGetErrorString(code);
    std::cerr << "CUDA error returned from \"" << func << "\" at "
	      << FILE << ":" << LINE
	      << " Error code: " << code << "(" << errorMessage << ")"
	      << std::endl;
  }
}

// For error reporting without error code
//
inline void
cuda_check_last_error_mpi(char const * func, const char *FILE, const int LINE,
			  const int myid)
{
  cudaError_t code = cudaGetLastError();
  if (code != cudaSuccess) {
    const char * errorMessage = cudaGetErrorString(code);
    std::cerr << "CUDA error returned from \"" << func
	      << ", mpi rank=" << myid
	      << "\" at " << FILE << ":" << LINE
	      << " Error code: " << code << "(" << errorMessage << ")"
	      << std::endl;
  }
}

// use CUDA's high-resolution timers when possible
//
struct cudaTimer
{
  cudaEvent_t start;
  cudaEvent_t end;

  cudaTimer(void)
  {
    cuda_safe_call(cudaEventCreate(&start), __FILE__, __LINE__);
    cuda_safe_call(cudaEventCreate(&end),   __FILE__, __LINE__);
    restart();
  }

  ~cudaTimer(void)
  {
    cuda_safe_call(cudaEventDestroy(start), __FILE__, __LINE__);
    cuda_safe_call(cudaEventDestroy(end),   __FILE__, __LINE__);
  }

  void restart(void)
  {
    cuda_safe_call(cudaEventRecord(start, 0), __FILE__, __LINE__);
  }

  double elapsed(void)
  {
    cuda_safe_call(cudaEventRecord(end, 0),   __FILE__, __LINE__);
    cuda_safe_call(cudaEventSynchronize(end), __FILE__, __LINE__);

    float ms_elapsed;
    cuda_safe_call(cudaEventElapsedTime(&ms_elapsed, start, end),
		   __FILE__, __LINE__);
    return ms_elapsed / 1e3;
  }

  double epsilon(void)
  {
    return 0.5e-6;
  }
};

inline
int ConvertSMVer2Cores(int major, int minor)
{
  // Defines for GPU Architecture types (using the SM version to
  // determine the # of cores per SM
  typedef struct
  {
    // 0xMm (hexidecimal notation), M = SM Major version, and m = SM
    // minor version
    int SM;
    int Cores;
  }
  sSMtoCores;
  
  sSMtoCores nGpuArchCoresPerSM[] =
    {
      { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
      { 0x32, 192}, // Kepler Generation (SM 3.2) GK10x class
      { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
      { 0x37, 192}, // Kepler Generation (SM 3.7) GK21x class
      { 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
      { 0x52, 128}, // Maxwell Generation (SM 5.2) GM20x class
      { 0x53, 128}, // Maxwell Generation (SM 5.3) GM20x class
      { 0x60, 64 }, // Pascal Generation (SM 6.0) GP100 class
      { 0x61, 128}, // Pascal Generation (SM 6.1) GP10x class
      { 0x62, 128}, // Pascal Generation (SM 6.2) GP10x class
      { 0x62, 128}, // Pascal Generation (SM 6.2) GP10x class
      { 0x70, 64 }, // Volta Generation (SM 7.0) GV100 class
      {   -1, -1 }
    };
  
  int index = 0;

  while (nGpuArchCoresPerSM[index].SM != -1) {
    
    if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
      return nGpuArchCoresPerSM[index].Cores;
    }

    index++;
  }

  // If we don't find the values, we default use the previous one to
  // run properly
  std::cout << "MapSMtoCores for SM " << major << "." << minor
	    << " is undefined.  The default is to use "
	    << nGpuArchCoresPerSM[index-1].Cores << "Cores/SM" << std::endl;
  
  return nGpuArchCoresPerSM[index-1].Cores;
}

// Convert CUDA int2 to double inline
//
__forceinline__ __device__ double int2_as_double (int2 a)
{
  double res;
  asm ("mov.b64 %0, {%1,%2};" : "=d"(res) : "r"(a.x), "r"(a.y));
  return res;
}

// Template structure to pass to kernel to avoid all of those pesky
// raw pointer casts
//
template <typename T>
struct dArray
{
  T*     _v;			// Pointer to underlying array data
  size_t _s;			// Number of elements
};
 
// Template function to convert device_vector to structure for passing
// from host to device in a cuda kernel
//
template <typename T>
dArray<T> toKernel(thrust::device_vector<T>& dev)
{
  return {thrust::raw_pointer_cast(&dev[0]), dev.size()};
}

// Template function to convert device_vector to structure for passing
// from host to device in a cuda kernel
//
template <typename T>
dArray<T> toKernelS(thrust::device_vector<T>& dev, int start, int size)
{
  return {thrust::raw_pointer_cast(&dev[start]), size};
}


// Select floating point type size for single/double precision switch
//
#ifdef O_SINGLE
#define cuREAL 4
#else
#define cuREAL 8
#endif

#if cuREAL == 4
typedef float  cuFP_t;
#else
typedef double cuFP_t;
#endif

// Work around for long-standing Thrust bug that is finally being
// addressed.  See https://github.com/NVIDIA/thrust/pull/1104
//
#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#elif (CUDART_VERSION >= 11010)
const bool thrust_binary_search_workaround = false;
#else
const bool thrust_binary_search_workaround = true;
#endif

#endif
