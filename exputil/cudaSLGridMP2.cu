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

void SLGridSph::initialize_cuda(std::vector<cudaArray_t>& cuArray,
				thrust::host_vector<cudaTextureObject_t>& tex)
{
  // Number of texture arrays
  //
  int ndim = (lmax+1)*nmax + 1;

  // Allocate CUDA array in device memory (a one-dimension 'channel')
  //
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  // cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

  // Interpolation data array
  //
  cuArray.resize(ndim);

  // Size of interpolation array
  //
  size_t tsize = numr*sizeof(float);

  // Create texture objects
  //
  tex.resize(ndim);
  thrust::fill(tex.begin(), tex.end(), 0);

  cudaResourceDesc resDesc;

  // Specify texture object parameters
  //
  cudaTextureDesc texDesc;

  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.filterMode = cudaFilterModePoint;
  texDesc.readMode = cudaReadModeElementType;
  texDesc.normalizedCoords = 0;

  thrust::host_vector<float> tt(numr);

  cuda_safe_call(cudaMallocArray(&cuArray[0], &channelDesc, numr, 1),
		 __FILE__, __LINE__, "malloc cuArray");

  // Copy to device memory some data located at address h_data
  // in host memory
  for (int j=0; j<numr; j++) tt[j] = p0[j];

  cuda_safe_call(cudaMemcpyToArray(cuArray[0], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

  // Specify texture
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeArray;
  resDesc.res.array.array = cuArray[0];

  cuda_safe_call(cudaCreateTextureObject(&tex[0], &resDesc, &texDesc, NULL),
		 __FILE__, __LINE__, "create texture object");

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      int i = 1 + l*nmax + n;
      cuda_safe_call(cudaMallocArray(&cuArray[i], &channelDesc, numr, 1),
		     __FILE__, __LINE__, "malloc cuArray");

      // Copy to device memory some data located at address h_data
      // in host memory
      float fac = sqrt(table[l].ev[n+1]);
      for (int j=0; j<numr; j++) tt[j] = table[l].ef[n+1][j] / fac;

      cuda_safe_call(cudaMemcpyToArray(cuArray[i], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

      // Specify texture
      resDesc.res.array.array = cuArray[i];

      cuda_safe_call(cudaCreateTextureObject(&tex[i], &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");
    }
  }

  if (false) {
    printf("**HOST** Texture compare\n");
    {
      for (int l : {0, 1, 2}) {
	for (int j=0; j<10; j++) {
	  for (int i : {3980, 3990, 3995, 3999}) 
	    printf("%5d %5d %5d %13.7e\n", l, j, i, table[l].ef[j+1][i]);
	}
      }
    }
  }

  if (false) {
    std::cout << "cuInterpArray size = " << cuArray.size() << std::endl;
    unsigned cnt = 0;
    for (size_t i=0; i<cuArray.size(); i++) {
      std::ostringstream sout;
      sout << "trying to free cuArray [" << cnt++ << "]";
      cuda_safe_call(cudaFreeArray(cuArray[i]), __FILE__, __LINE__, sout.str());
    }
    
    std::cout << "texture object array size = " << tex.size() << std::endl;
    for (size_t i=0; i<tex.size(); i++) {
      std::ostringstream sout;
      sout << "trying to free TextureObject [" << cnt++ << "]";
      cuda_safe_call(cudaDestroyTextureObject(tex[i]), __FILE__, __LINE__, sout.str());
    }
    
    std::cout << "cuda memory freed" << std::endl;
    exit(-1);
  }
}

