// -*- C++ -*-

#include <EmpCyl2D.H>

#include <iostream>
#include <iomanip>
#include <limits>
#include <map>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>

#include <cudaUtil.cuH>

__global__
void testFetchCyl2d(dArray<cudaTextureObject_t> T, dArray<cuFP_t> f,
		    int m, int j, int nmax, int numr)
{
  const int n = blockDim.x * blockIdx.x + threadIdx.x;
  const int k = m*nmax;
#if cuREAL == 4
  if (n < numr) f._v[n] = tex1D<float>(T._v[k+j], n);
#else
  if (n < numr) f._v[n] = int2_as_double(tex1D<int2>(T._v[k+j], n));
#endif
}


thrust::host_vector<cuFP_t> returnTestCyl2d
(thrust::host_vector<cudaTextureObject_t>& tex,
 int m, int j, int nmax, int numr)
{
  thrust::device_vector<cudaTextureObject_t> t_d = tex;
  
  unsigned int gridSize  = numr/BLOCK_SIZE;
  if (numr > gridSize*BLOCK_SIZE) gridSize++;
  
  thrust::device_vector<cuFP_t> f_d(numr);

  testFetchCyl2d<<<gridSize, BLOCK_SIZE>>>(toKernel(t_d), toKernel(f_d),
					   m, j, nmax, numr);

  cudaDeviceSynchronize();

  return f_d;
}

static std::vector<cudaResourceDesc> resDesc;
static std::vector<cudaTextureDesc>  texDesc;

struct Element {
  int m;
  double a;
  double b;
};

void EmpCyl2D::initialize_cuda(std::vector<cudaArray_t>& cuArray,
			       thrust::host_vector<cudaTextureObject_t>& tex)
{
  // Number of texture arrays
  //
  int ndim = (mmax+1)*nmax + 1;

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

  for (int m=0; m<=mmax; l++) {
    for (int n=0; n<nmax; n++) {
      int i = m*nmax + n;
      cuda_safe_call(cudaMallocArray(&cuArray[i], &channelDesc, numr), __FILE__, __LINE__, "malloc cuArray");

      // Copy to device memory some data located at address h_data
      // in host memory
      for (int j=0; j<numr; j++) tt[j] = potl_array[m](j, n);

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
    for (int m=0; m<=mmax; m++) {
      for (int j=0; j<nmax; j++) {
	ret = returnTestCyl2d(tex, m, j, nmax, numr);
	for (int i=0; i<numr; i++) {
	  cuFP_t a = potl_array[m](j, n);
	  cuFP_t b = ret[i];
	  if (a>1.0e-18) {

	    Element comp = {m, a, b};
	    compare.insert(std::make_pair(fabs((a - b)/a), comp));

	    if ( fabs((a - b)/a ) > tol) {
	      std::cout << std::setw( 5) << m << std::setw( 5) << j
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
	      << "**Low[1] : " << lo1->first << " (" << lo1->second.m << ", " << lo1->second.a << ", " << lo1->second.b << ", " << lo1->second.a - lo1->second.b << ")" << std::endl
	      << "**Low[9] : " << lo9->first << " (" << lo9->second.m << ", " << lo9->second.a << ", " << lo9->second.b << ", " << lo9->second.a - lo9->second.b << ")" << std::endl
	      << "**Middle : " << mid->first << " (" << mid->second.m << ", " << mid->second.a << ", " << mid->second.b << ", " << mid->second.a - mid->second.b << ")" << std::endl
	      << "**Hi [9] : " << hi9->first << " (" << hi9->second.m << ", " << hi9->second.a << ", " << hi9->second.b << ", " << hi9->second.a - hi9->second.b << ")" << std::endl
	      << "**Hi [1] : " << hi1->first << " (" << hi1->second.m << ", " << hi1->second.a << ", " << hi1->second.b << ", " << hi1->second.a - hi1->second.b << ")" << std::endl
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

