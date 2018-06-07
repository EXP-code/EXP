#ifndef _BOUNDS_CUH_
#define _BOUNDS_CUH_

#include <limits>

#include <thrust/tuple.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>

template<typename T>
__device__ __host__
int binary_search(T A[], T val, int imin, int imax, int flag)
{
  bool lo, hi;

  if (flag<0) {			// Lower bound
    lo = val - A[imin] >  0;
    hi = A[imax] - val >= 0;
  } else {			// Upper bound
    lo = val - A[imin] >= 0;
    hi = A[imax] - val >  0;
  }

  // We've got it
  if (lo and hi) {

    bool test;
    while (imax > imin+1) {
      int imid = (imax + imin)/2;
      if (flag<0) test = val - A[imid] >  0;
      else        test = val - A[imid] >= 0;
      if (test) {
	imin = imid;
      } else {
	imax = imid;
      }
    }

    return imax;
  }

  return -1;
}

enum class BoundType {Lower, Upper};

template<typename T, BoundType B>
__global__
void search(T * source, int * result, int Nin, int Nout, T value)
{
  extern __shared__ T buff[];
  int tid = threadIdx.x + blockIdx.x*(blockDim.x-1);
  
  T val = std::numeric_limits<T>::max();
  if (tid < Nin) val = source[tid];
  buff[threadIdx.x]  = val;
  __syncthreads();

  if (threadIdx.x==0) {

    int res;

    switch (B) {
    case BoundType::Lower:
      res = binary_search<T>(buff, value, 0, blockDim.x-1, -1);
      break;
    
    case BoundType::Upper:
      res = binary_search<T>(buff, value, 0, blockDim.x-1, 1);
      break;
    }
    
    if (blockIdx.x < Nout) {
      if (res>=0)
	result[blockIdx.x] = res + blockIdx.x*(blockDim.x-1);
      else
	result[blockIdx.x] = -1;
    }
  }
}

// return the biggest of two tuples
//
template <class T>
struct bigger_tuple
{
  __device__ __host__
  thrust::tuple<T,int> operator()(const thrust::tuple<T,int> &a, const thrust::tuple<T,int> &b) 
  {
    if (a > b) return a;
    else return b;
  } 
};

// Return the index of the largest element
//
template <class T>
thrust::tuple<T,int> maxVec(thrust::device_vector<T>& vec, cudaStream_t stream)
{
  // create implicit index sequence [0, 1, 2, ... )
  thrust::counting_iterator<int> begin(0);
  thrust::counting_iterator<int> end(vec.size());
  thrust::tuple<T,int> init(vec[0], 0); 
  thrust::tuple<T,int> smallest;
  
  smallest = thrust::reduce
    (thrust::cuda::par.on(stream),
     thrust::make_zip_iterator(thrust::make_tuple(vec.begin(), begin)),
     thrust::make_zip_iterator(thrust::make_tuple(vec.end(), end)),
     init, bigger_tuple<T>());

  return smallest;
}


template<typename T>
int getBound(T value, thrust::device_vector<T>& sorted, cudaStream_t stream, BoundType B)
{
  const int inputLength = sorted.size();
  
  T beg = sorted[0];
  T end = sorted[inputLength-1];

  if (value < beg) {
    return 0;
  }

  if (value > end) {
    return inputLength;
  }

  const int blockSize   = 256;
  const int outputLength = inputLength/(blockSize-1) +
    (inputLength % (blockSize-1) > 0 ? 1 : 0);

  thrust::device_vector<int> devOutput(outputLength);

  dim3 DimBlock(blockSize, 1, 1);
  dim3 DimGrid(1, 1, 1);

  DimGrid.x = (inputLength / (DimBlock.x-1)) + ((inputLength % (DimBlock.x-1) > 0) ? 1 : 0); 

  size_t shmsz = DimBlock.x * sizeof(T);
  
  switch (B) {
  case BoundType::Upper:
    search<T, BoundType::Upper><<<DimGrid, DimBlock, shmsz, stream>>>
      (thrust::raw_pointer_cast(sorted.data()),
       thrust::raw_pointer_cast(devOutput.data()),
       inputLength, outputLength, value);
    break;
  case BoundType::Lower:
    search<T, BoundType::Lower><<<DimGrid, DimBlock, shmsz, stream>>>
      (thrust::raw_pointer_cast(sorted.data()),
       thrust::raw_pointer_cast(devOutput.data()),
       inputLength, outputLength, value);
  }    
  
  thrust::tuple<int, int> answer = maxVec(devOutput, stream);
  
  int pos = thrust::get<0>(answer);
  if (pos < 0) pos = 0;
  if (pos > inputLength) pos = inputLength;

  return pos;
}


#endif
