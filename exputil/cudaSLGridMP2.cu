#include <SLGridMP2.h>

#include <iostream>
#include <iomanip>
#include <map>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>

__constant__ float cuRscale, cuXmin, cuXmax, cuDxi;
__constant__ int   cuNumr, cuCmap;

void SLGridSph::initialize_cuda(cudaChannelFormatDesc& channelDesc,
				std::vector<cudaArray*>& cuArray,
				std::vector<cudaResourceDesc>& resDesc,
				struct cudaTextureDesc& texDesc,
				thrust::host_vector<cudaTextureObject_t>& tex
				)
{
  // Number of texture arrays
  //
  int ndim = (lmax+1)*nmax;

  // Allocate CUDA array in device memory
  channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

  cuArray.resize(ndim);
  size_t tsize = BLOCK_SIZE*sizeof(float);

  // Create texture objects
  tex.resize(ndim);
  thrust::fill(tex.begin(), tex.end(), 0);

  resDesc.resize(ndim);

  // Specify texture object parameters
  //
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.filterMode = cudaFilterModePoint;
  texDesc.readMode = cudaReadModeElementType;
  texDesc.normalizedCoords = 0;

  thrust::host_vector<float> tt(numr);

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      int i = l*nmax + n;
      cuda_safe_call(cudaMallocArray(&cuArray[i], &channelDesc, BLOCK_SIZE, 1), "malloc cuArray");

      // Copy to device memory some data located at address h_data
      // in host memory
      for (int j=0; j<numr; j++) tt[j] = table[l].ef[n+1][j];

      cuda_safe_call(cudaMemcpyToArray(cuArray[i], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), "copy texture to array");

      // Specify texture
      memset(&resDesc[i], 0, sizeof(resDesc));
      resDesc[i].resType = cudaResourceTypeArray;
      resDesc[i].res.array.array = cuArray[i];

      cuda_safe_call(cudaCreateTextureObject(&tex[i], &resDesc[i], &texDesc, NULL), "create texture object");
    }
  }
}
