// -*- C++ -*-

#include <SLGridMP2.H>

#include <iostream>
#include <iomanip>
#include <limits>
#include <map>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#if CUDART_VERSION < 12000
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#endif

#include <cudaUtil.cuH>

__global__
void testFetchSph(dArray<cudaTextureObject_t> T, dArray<cuFP_t> f,
		  int l, int j, int nmax, int numr)
{
  const int n = blockDim.x * blockIdx.x + threadIdx.x;
  const int k = l*nmax + 1;
#if cuREAL == 4
  if (n < numr) f._v[n] = tex1D<float>(T._v[k+j], n);
#else
  if (n < numr) f._v[n] = int2_as_double(tex1D<int2>(T._v[k+j], n));
#endif
}


thrust::host_vector<cuFP_t> returnTestSph
(thrust::host_vector<cudaTextureObject_t>& tex,
 int l, int j, int nmax, int numr)
{
  thrust::device_vector<cudaTextureObject_t> t_d = tex;
  
  unsigned int gridSize  = numr/BLOCK_SIZE;
  if (numr > gridSize*BLOCK_SIZE) gridSize++;
  
  thrust::device_vector<cuFP_t> f_d(numr);

  testFetchSph<<<gridSize, BLOCK_SIZE>>>(toKernel(t_d), toKernel(f_d),
					 l, j, nmax, numr);

  cudaDeviceSynchronize();

  return f_d;
}

struct Element {
  int l;
  double a;
  double b;
};

void SLGridSph::initialize_cuda(std::vector<cudaArray_t>& cuArray,
				thrust::host_vector<cudaTextureObject_t>& tex)
{
  // Number of texture arrays
  //
  int ndim = (lmax+1)*nmax + 1;

  // Allocate CUDA array in device memory (a one-dimension 'channel')
  //
#if cuREAL == 4
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
#else
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
#endif

  // Interpolation data array
  //
  cuArray.resize(ndim);

  // Size of interpolation array
  //
  size_t tsize = numr*sizeof(cuFP_t);

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

  thrust::host_vector<cuFP_t> tt(numr);

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
      cuFP_t fac = sqrt(table[l].ev[n]);
      for (int j=0; j<numr; j++) tt[j] = table[l].ef(n, j) / fac;

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

  // This is for debugging: compare texture table fetches to original
  // tables
  //
  if (false) {
    const cuFP_t tol = 10.0*std::numeric_limits<cuFP_t>::epsilon();

    std::multimap<double, Element> compare;

    thrust::host_vector<cuFP_t> ret(numr);
    std::cout << "**HOST** Texture compare" << std::endl << std::scientific;
    unsigned tot = 0, bad = 0;
    for (int l=0; l<=lmax; l++) {
      for (int j=0; j<nmax; j++) {
	ret = returnTestSph(tex, l, j, nmax, numr);
	for (int i=0; i<numr; i++) {
	  cuFP_t a = table[l].ef(j, i)/sqrt(table[l].ev(j));
	  cuFP_t b = ret[i];
	  if (a>1.0e-18) {

	    Element comp = {l, a, b};
	    compare.insert(std::make_pair(fabs((a - b)/a), comp));

	    if ( fabs((a - b)/a ) > tol) {
	      std::cout << std::setw( 5) << l << std::setw( 5) << j
			<< std::setw( 5) << i << std::setw(20) << a
			<< std::setw(20) << (a-b)/a << std::endl;
	      bad++;
	    }
	  }
	  tot++;
	}
      }
    }

    std::multimap<double, Element>::iterator beg = compare.begin();
    std::multimap<double, Element>::iterator end = compare.end();
    std::multimap<double, Element>::iterator
      lo1=beg, lo9=beg, mid=beg, hi9=end, hi1=end;

    std::advance(lo9, 9);
    std::advance(mid, compare.size()/2);
    std::advance(hi1, -1);
    std::advance(hi9, -10);

    std::cout << "**Found " << bad << "/" << tot << " values > eps" << std::endl
	      << "**Low[1] : " << lo1->first << " (" << lo1->second.l << ", " << lo1->second.a << ", " << lo1->second.b << ", " << lo1->second.a - lo1->second.b << ")" << std::endl
	      << "**Low[9] : " << lo9->first << " (" << lo9->second.l << ", " << lo9->second.a << ", " << lo9->second.b << ", " << lo9->second.a - lo9->second.b << ")" << std::endl
	      << "**Middle : " << mid->first << " (" << mid->second.l << ", " << mid->second.a << ", " << mid->second.b << ", " << mid->second.a - mid->second.b << ")" << std::endl
	      << "**Hi [9] : " << hi9->first << " (" << hi9->second.l << ", " << hi9->second.a << ", " << hi9->second.b << ", " << hi9->second.a - hi9->second.b << ")" << std::endl
	      << "**Hi [1] : " << hi1->first << " (" << hi1->second.l << ", " << hi1->second.a << ", " << hi1->second.b << ", " << hi1->second.a - hi1->second.b << ")" << std::endl
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

