#include <SLGridMP2.h>

#include <iostream>
#include <iomanip>
#include <map>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>

#include <cudaUtil.cuH>

__constant__ double cuRscale, cuXmin, cuXmax, cuDxi;
__constant__ int   cuNumr, cuCmap;


__global__
void testFetchSph(dArray<cudaTextureObject_t> T, dArray<double> f,
		  int l, int j, int nmax, int numr)
{
  const int n = blockDim.x * blockIdx.x + threadIdx.x;
  const int k = l*nmax + 1;
  if (n < numr) f._v[n] = int2_as_double(tex1D<int2>(T._v[k+j], n));
}


thrust::host_vector<double> returnTestSph
(thrust::host_vector<cudaTextureObject_t>& tex,
 int l, int j, int nmax, int numr)
{
  thrust::device_vector<cudaTextureObject_t> t_d = tex;
  
  unsigned int gridSize  = numr/BLOCK_SIZE;
  if (numr > gridSize*BLOCK_SIZE) gridSize++;

  thrust::device_vector<double> f_d(numr);

  testFetchSph<<<gridSize, BLOCK_SIZE>>>(toKernel(t_d), toKernel(f_d),
					 l, j, nmax, numr);

  cudaDeviceSynchronize();

  return f_d;
}

static std::vector<cudaResourceDesc> resDesc;
static std::vector<cudaTextureDesc>  texDesc;


void SLGridSph::initialize_cuda(std::vector<cudaArray_t>& cuArray,
				thrust::host_vector<cudaTextureObject_t>& tex)
{
  // Number of texture arrays
  //
  int ndim = (lmax+1)*nmax + 1;

  // Allocate CUDA array in device memory (a one-dimension 'channel')
  //
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();

  // Interpolation data array
  //
  cuArray.resize(ndim);

  // Size of interpolation array
  //
  size_t tsize = numr*sizeof(double);

  // Create texture objects
  //
  tex.resize(ndim);
  thrust::fill(tex.begin(), tex.end(), 0);

  // cudaResourceDesc resDesc;
  resDesc.resize(ndim);

  // Specify texture object parameters
  //
  texDesc.resize(ndim);

  memset(&texDesc[0], 0, sizeof(cudaTextureDesc));
  texDesc[0].addressMode[0] = cudaAddressModeClamp;
  texDesc[0].filterMode = cudaFilterModePoint;
  texDesc[0].readMode = cudaReadModeElementType;
  texDesc[0].normalizedCoords = 0;

  thrust::host_vector<double> tt(numr);

  cuda_safe_call(cudaMallocArray(&cuArray[0], &channelDesc, numr), __FILE__, __LINE__, "malloc cuArray");

  // Copy to device memory some data located at address h_data
  // in host memory
  for (int j=0; j<numr; j++) tt[j] = p0[j];

  cuda_safe_call(cudaMemcpyToArray(cuArray[0], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

  // Specify texture
  memset(&resDesc[0], 0, sizeof(cudaResourceDesc));
  resDesc[0].resType = cudaResourceTypeArray;
  resDesc[0].res.array.array = cuArray[0];

  cuda_safe_call(cudaCreateTextureObject(&tex[0], &resDesc[0], &texDesc[0], NULL), __FILE__, __LINE__, "create texture object");

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      int i = 1 + l*nmax + n;
      cuda_safe_call(cudaMallocArray(&cuArray[i], &channelDesc, numr), __FILE__, __LINE__, "malloc cuArray");

      // Copy to device memory some data located at address h_data
      // in host memory
      double fac = sqrt(table[l].ev[n+1]);
      for (int j=0; j<numr; j++) tt[j] = table[l].ef[n+1][j] / fac;

      cuda_safe_call(cudaMemcpyToArray(cuArray[i], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");
      
      // Specify texture
      memset(&resDesc[i], 0, sizeof(cudaResourceDesc));
      resDesc[i].resType = cudaResourceTypeArray;
      resDesc[i].res.array.array = cuArray[i];

      memset(&texDesc[i], 0, sizeof(cudaTextureDesc));
      texDesc[i].addressMode[0] = cudaAddressModeClamp;
      texDesc[i].filterMode = cudaFilterModePoint;
      texDesc[i].readMode = cudaReadModeElementType;
      texDesc[i].normalizedCoords = 0;

      cuda_safe_call(cudaCreateTextureObject(&tex[i], &resDesc[i], &texDesc[i], NULL), __FILE__, __LINE__, "create texture object");
    }
  }

  if (true) {
    thrust::host_vector<double> ret(numr);
    std::cout << "**HOST** Texture compare" << std::endl;
    unsigned tot = 0, bad = 0;
    for (int l=0; l<=lmax; l++) {
      for (int j=0; j<nmax; j++) {
	ret = returnTestSph(tex, l, j, nmax, numr);
	for (int i=0; i<numr; i++) {
	  double a = table[l].ef[j+1][i]/sqrt(table[l].ev[j+1]);
	  double b = ret[i];
	  if (a>1.0e-18) {
	    if ( fabs((a - b)/a ) > 1.0e-7) {
	      std::cout << std::setw( 5) << l << std::setw( 5) << j
			<< std::setw( 5) << i << std::setw(15) << a
			<< std::setw(15) << (a-b)/a << std::endl;
	      bad++;
	    }
	  }
	  tot++;
	}
      }
    }
    std::cout << "**Found " << bad << "/" << tot << " bad values" << std::endl
	      << "**" << std::endl;
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

