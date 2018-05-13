#ifndef CUDA_REDUCE_CUH
#define CUDA_REDUCE_CUH

#include <cudaUtil.cuH>
#include <thrust/iterator/discard_iterator.h>

#define FSMALL 1.0e-16
#define MINEPS 1.0e-10

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

typedef std::pair<unsigned int, unsigned int> PII;

// Loop unroll
//
template <typename T, unsigned int blockSize>
__device__
void warpReduce(volatile T *sdata, unsigned int tid)
{
  if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
  if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
  if (blockSize >= 16) sdata[tid] += sdata[tid +  8];
  if (blockSize >=  8) sdata[tid] += sdata[tid +  4];
  if (blockSize >=  4) sdata[tid] += sdata[tid +  2];
  if (blockSize >=  2) sdata[tid] += sdata[tid +  1];
}

template <typename T, unsigned int blockSize>
__global__ void reduceSum(dArray<float> out, dArray<float> in,
			  unsigned int dim, unsigned int n)
{
  extern __shared__ T sdata[];

  unsigned int tid      = threadIdx.x;
  unsigned int gridSize = blockSize*gridDim.x*2;
    
  for (unsigned j=0; j<dim; j++) {

    sdata[tid] = 0;
    
    unsigned int i = blockIdx.x*blockSize*2 + tid;

    while (i < n) {
      // Sanity check
      if (i+n*j>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
      if (i+blockSize<n)
	if (i+n*j+blockSize>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

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
      warpReduce<T, blockSize>(&sdata[0], tid);
    }

    if (tid == 0) {
      if (blockIdx.x + j*gridDim.x>=out._s)
	printf("reduceSum: out of bounds, b=%d/%d j=%d\n", blockIdx.x, gridDim.x, j);
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


__device__ __constant__ float cuRscale, cuHscale, cuXmin, cuXmax, cuYmin, cuYmax, cuDxi, cuDyi;

__device__ __constant__ int   cuNumr, cuNumx, cuNumy, cuCmap;


#endif
