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

  // std::vector<cudaResourceDesc> resDesc;
  // std::vector<cudaTextureDesc>  texDesc;


  cudaTextureDesc texDesc;

  memset(&texDesc, 0, sizeof(cudaTextureDesc));
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.filterMode = cudaFilterModePoint;
  texDesc.readMode = cudaReadModeElementType;
  texDesc.normalizedCoords = 0;

  thrust::host_vector<cuFP_t> tt(numr);

  cuda_safe_call(cudaMallocArray(&cuArray[0], &channelDesc, numr), __FILE__, __LINE__, "malloc cuArray");

  // Copy to device memory some data located at address h_data
  // in host memory
  for (int j=0; j<numr; j++) tt[j] = p0[j];

  cuda_safe_call(cudaMemcpyToArray(cuArray[0], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

  // Specify texture
  cudaResourceDesc resDesc;

  memset(&resDesc, 0, sizeof(cudaResourceDesc));
  resDesc.resType = cudaResourceTypeArray;
  resDesc.res.array.array = cuArray[0];

  cuda_safe_call(cudaCreateTextureObject(&tex[0], &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");

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
      cudaResourceDesc resDesc;

      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array = cuArray[i];

      cudaTextureDesc texDesc;

      memset(&texDesc, 0, sizeof(cudaTextureDesc));
      texDesc.addressMode[0] = cudaAddressModeClamp;
      texDesc.filterMode = cudaFilterModePoint;
      texDesc.readMode = cudaReadModeElementType;
      texDesc.normalizedCoords = 0;

      cuda_safe_call(cudaCreateTextureObject(&tex[i], &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");
    }
  }

  // This is for debugging: compare texture table fetches to original
  // tables
  //
  if (false) {
    const cuFP_t tol = 10.0*std::numeric_limits<cuFP_t>::epsilon();

    struct Element {
      int l;
      double a;
      double b;
    };

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


__global__
void testFetchSlab(dArray<cudaTextureObject_t> T, dArray<cuFP_t> f,
		  int l, int j, int nmax, int numz)
{
  const int n = blockDim.x * blockIdx.x + threadIdx.x;
  const int k = l*nmax + 1;
#if cuREAL == 4
  if (n < numz) f._v[n] = tex1D<float>(T._v[k+j], n);
#else
  if (n < numz) f._v[n] = int2_as_double(tex1D<int2>(T._v[k+j], n));
#endif
}


thrust::host_vector<cuFP_t> returnTestSlab
(thrust::host_vector<cudaTextureObject_t>& tex,
 int kx, int ky, int j, int numk, int nmax, int numz)
{
  thrust::device_vector<cudaTextureObject_t> t_d = tex;
  
  unsigned int gridSize  = numz/BLOCK_SIZE;
  if (numz > gridSize*BLOCK_SIZE) gridSize++;
  
  thrust::device_vector<cuFP_t> f_d(numz);

  int l = kx*(kx+1)/2 + ky;

  testFetchSlab<<<gridSize, BLOCK_SIZE>>>(toKernel(t_d), toKernel(f_d),
					  l, j, nmax, numz);
  
  cudaDeviceSynchronize();

  return f_d;
}

void SLGridSlab::initialize_cuda(std::vector<cudaArray_t>& cuArray,
				 thrust::host_vector<cudaTextureObject_t>& tex)
{
  // Number of texture arrays
  //
  int ndim = (numk+1)*(numk+2)/2*nmax + 1;

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
  size_t tsize = numz*sizeof(cuFP_t);

  // Create texture objects
  //
  tex.resize(ndim);
  thrust::fill(tex.begin(), tex.end(), 0);

  // std::vector<cudaResourceDesc> resDesc;
  // std::vector<cudaTextureDesc>  texDesc;


  cudaTextureDesc texDesc;

  memset(&texDesc, 0, sizeof(cudaTextureDesc));
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.filterMode = cudaFilterModePoint;
  texDesc.readMode = cudaReadModeElementType;
  texDesc.normalizedCoords = 0;

  thrust::host_vector<cuFP_t> tt(numz);

  cuda_safe_call(cudaMallocArray(&cuArray[0], &channelDesc, numz), __FILE__, __LINE__, "malloc cuArray");

  // Copy to device memory some data located at address h_data
  // in host memory
  for (int j=0; j<numz; j++) tt[j] = p0[j];

  cuda_safe_call(cudaMemcpyToArray(cuArray[0], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

  // Specify texture
  cudaResourceDesc resDesc;

  memset(&resDesc, 0, sizeof(cudaResourceDesc));
  resDesc.resType = cudaResourceTypeArray;
  resDesc.res.array.array = cuArray[0];

  cuda_safe_call(cudaCreateTextureObject(&tex[0], &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");

  int i = 1;
  for (int kx=0; kx<=numk; kx++) {
    for (int ky=0; ky<=kx; ky++) {
      for (int n=0; n<nmax; n++, i++) {
	cuda_safe_call(cudaMallocArray(&cuArray[i], &channelDesc, numz), __FILE__, __LINE__, "malloc cuArray");

	// Copy to device memory some data located at address h_data
	// in host memory
	cuFP_t fac = sqrt(table[kx][ky].ev[n]);
	for (int j=0; j<numz; j++) tt[j] = table[kx][ky].ef(n, j) / fac;

	cuda_safe_call(cudaMemcpyToArray(cuArray[i], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");
      
	// Specify texture
	cudaResourceDesc resDesc;

	memset(&resDesc, 0, sizeof(cudaResourceDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = cuArray[i];

	cudaTextureDesc texDesc;

	memset(&texDesc, 0, sizeof(cudaTextureDesc));
	texDesc.addressMode[0] = cudaAddressModeClamp;
	texDesc.filterMode = cudaFilterModePoint;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = 0;
	
	cuda_safe_call(cudaCreateTextureObject(&tex[i], &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");
      }
      // vertical loop
    }
    // y axis loop
  }
  // x axis loop

  // This is for debugging: compare texture table fetches to original
  // tables
  //
  if (false) {
    const cuFP_t tol = 10.0*std::numeric_limits<cuFP_t>::epsilon();

    struct Element {
      int kx;
      int ky;
      int j;
      double a;
      double b;
    };

    std::multimap<double, Element> compare;

    thrust::host_vector<cuFP_t> ret(numz);
    std::cout << "**HOST** Texture compare" << std::endl << std::scientific;
    unsigned tot = 0, bad = 0;

    for (int kx=0; kx<=numk; kx++) {

      for (int ky=0; ky<=kx; ky++) {

	for (int j=0; j<nmax; j++) {

	  ret = returnTestSlab(tex, kx, ky, j, numk, nmax, numz);
	  for (int i=0; i<numz; i++) {
	    cuFP_t a = table[kx][ky].ef(j, i)/sqrt(table[kx][ky].ev(j));
	    cuFP_t b = ret[i];
	    if (a>1.0e-18) {

	      Element comp = {kx, ky, j, a, b};
	      compare.insert(std::make_pair(fabs((a - b)/a), comp));

	      double test = (a - b)/a;

	      if ( fabs(test) > tol) {
		std::cout << std::setw( 5) << kx << std::setw( 5) << ky
			  << std::setw( 5) << j  << std::setw( 5) << i
			  << std::setw(20) << a
			  << std::setw(20) << test << std::endl;
		bad++;
	      }
	    }
	    tot++;
	  }
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
	      << "**Low[1] : " << lo1->first << " (" << lo1->second.kx << ", " << lo1->second.ky << ", " << lo1->second.j << ", " << lo1->second.a << ", " << lo1->second.b << ", " << lo1->second.a - lo1->second.b << ")" << std::endl
	      << "**Low[9] : " << lo9->first << " (" << lo9->second.kx << ", " << lo9->second.ky << ", " << lo9->second.j << ", " << lo9->second.a << ", " << lo9->second.b << ", " << lo9->second.a - lo9->second.b << ")" << std::endl
	      << "**Middle : " << mid->first << " (" << mid->second.kx << ", " << mid->second.ky << ", " << mid->second.j << ", " << mid->second.a << ", " << mid->second.b << ", " << mid->second.a - mid->second.b << ")" << std::endl
	      << "**Hi [9] : " << hi9->first << " (" << hi9->second.kx << ", " << hi9->second.ky << ", " << hi9->second.j << ", " << hi9->second.a << ", " << hi9->second.b << ", " << hi9->second.a - hi9->second.b << ")" << std::endl
	      << "**Hi [1] : " << hi1->first << " (" << hi1->second.kx << ", " << hi1->second.ky << ", " << hi1->second.j << ", " << hi1->second.a << ", " << hi1->second.b << ", " << hi1->second.a - hi1->second.b << ")" << std::endl
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

