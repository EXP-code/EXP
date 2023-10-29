// -*- C++ -*-

#ifndef CUDA_REDUCE_CUH
#define CUDA_REDUCE_CUH

#include <cudaUtil.cuH>
#include <thrust/iterator/discard_iterator.h>

// Machine constants
//
constexpr cuFP_t FSMALL = 20.0 * std::numeric_limits<cuFP_t>::min();
constexpr cuFP_t MINEPS = 20.0 * std::numeric_limits<cuFP_t>::min();

typedef std::pair<unsigned int, unsigned int> PII;
// Loop unroll
//
// Assumes blocksize is 64 or larger power of 2 up to max of 512 (or
// 1024 see below)
//
// https://developer.nvidia.com/blog/using-cuda-warp-level-primitives/
//
template <typename T>
__device__ void warpReduce(T *sdata, unsigned int tid) {

  T v = sdata[tid]; 
  v += sdata[tid+32]; 
  sdata[tid] = v;     __syncwarp();
  v += sdata[tid+16]; __syncwarp();
  sdata[tid] = v;     __syncwarp();
  v += sdata[tid+8];  __syncwarp();
  sdata[tid] = v;     __syncwarp();
  v += sdata[tid+4];  __syncwarp();
  sdata[tid] = v;     __syncwarp();
  v += sdata[tid+2];  __syncwarp();
  sdata[tid] = v;     __syncwarp();
  v += sdata[tid+1];  __syncwarp();
  sdata[tid] = v;
}

template <typename T, unsigned int blockSize>
__global__ void reduceSum(dArray<T> out, dArray<T> in,
			  unsigned int dim, unsigned int n)
{
  extern __shared__ T sdata[];

  unsigned int tid      = threadIdx.x;
  unsigned int gridSize = blockSize*gridDim.x*2;
    
  for (unsigned j=0; j<dim; j++) {

    sdata[tid] = 0;
    
    unsigned int i = blockIdx.x*blockSize*2 + tid;

    while (i < n) {
      // Sanity checks
      //
      if (i+n*j>=in._s) {
	printf("out of bounds: %d >= %d @ %s:%d\n",
	       i+n*j, in._s, __FILE__, __LINE__);
      }

      if (i+blockSize<n) {
	if (i+n*j+blockSize>=in._s)
	  printf("out of bounds: block %d >= %d @ %s:%d\n",
		 i+n*j+blockSize, in._s, __FILE__, __LINE__);
      }

      sdata[tid] +=
	in._v[i + n*j] + (i+blockSize<n ? in._v[i + blockSize + n*j] : T(0));
      i += gridSize;
    }
  
    __syncthreads();
  
    if (blockSize >= 1024) {
      if (tid < 512) {
	sdata[tid] += sdata[tid + 512];
	__syncthreads();
      }
    }
    
    if (blockSize >= 512) {
      if (tid < 256) {
	sdata[tid] += sdata[tid + 256];
      }
      __syncthreads();
    }

    if (blockSize >= 256) {
      if (tid < 128) {
	sdata[tid] += sdata[tid + 128];
      }
      __syncthreads();
    }

    if (blockSize >= 128) {
      if (tid < 64) {
	sdata[tid] += sdata[tid + 64];
      }
      __syncthreads();
    }

    if (tid < 32) {
      warpReduce<T>(&sdata[0], tid);
    }

    if (tid == 0) {
      if (blockIdx.x + j*gridDim.x>=out._s) {
	printf("reduceSum: out of bounds, b=%d/%d j=%d\n", blockIdx.x, gridDim.x, j);
      }
      
     out._v[blockIdx.x + j*gridDim.x] = sdata[tid];
    }
    __syncthreads();
  }
}

template <typename T, unsigned int blockSize>
__global__ void reduceSumS(dArray<T> out, dArray<T> in,
			   unsigned int dim, unsigned int n,
			   unsigned int beg, unsigned int e)
{
  extern __shared__ T sdata[];

  unsigned int tid      = threadIdx.x;
  unsigned int gridSize = blockSize*gridDim.x*2;
    
  for (unsigned j=0; j<dim; j++) {

    sdata[tid] = 0;
    
    unsigned int i = blockIdx.x*blockSize*2 + beg + tid;

    while (i < e) {
      // Sanity checks
      //
      if (i+n*j>=in._s) {
	printf("reduceSumS out of bounds: %d >= %d @ %s:%d\n",
	       i+n*j, in._s, __FILE__, __LINE__);
      }

      sdata[tid] +=
	in._v[i + n*j] + (i+blockSize<e ? in._v[i + blockSize + n*j] : T(0));
      i += gridSize;
    }
  
    __syncthreads();
  
    if (blockSize >= 1024) {
      if (tid < 512) {
	sdata[tid] += sdata[tid + 512];
	__syncthreads();
      }
    }
    
    if (blockSize >= 512) {
      if (tid < 256) {
	sdata[tid] += sdata[tid + 256];
      }
      __syncthreads();
    }

    if (blockSize >= 256) {
      if (tid < 128) {
	sdata[tid] += sdata[tid + 128];
      }
      __syncthreads();
    }

    if (blockSize >= 128) {
      if (tid < 64) {
	sdata[tid] += sdata[tid + 64];
      }
      __syncthreads();
    }

    if (tid < 32) {
      warpReduce<T>(&sdata[0], tid);
    }

    if (tid == 0) {
      if (blockIdx.x + j*gridDim.x>=out._s) {
	printf("reduceSumS: out of bounds, b=%d/%d j=%d\n", blockIdx.x, gridDim.x, j);
      }
      
     out._v[blockIdx.x + j*gridDim.x] = sdata[tid];
    }
    __syncthreads();
  }
}

struct key_functor : public thrust::unary_function<int, int>
{
  const int Ngrid;
  key_functor(int _Ngrid) : Ngrid(_Ngrid) {}
 
  __host__ __device__
  int operator()(int x) { return x / Ngrid; }
};


#endif
