// -*- C++ -*-

#include <EmpCylSL.H>

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
void testFetchCyl(dArray<cudaTextureObject_t> T, dArray<cuFP_t> f,
		  int m, int n, int k, int nmax, int NUMX, int NUMY)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int l = m*nmax + n;
  const int j = tid/NUMX;
  const int i = tid - NUMX*j;
  if (i<NUMX and j<NUMY) {
#if cuREAL == 4
    f._v[tid] = tex3D<float>(T._v[l], i, j, k);
#else
    f._v[tid] = int2_as_double(tex3D<int2>(T._v[l], i, j, k));
#endif
  }
}

thrust::host_vector<cuFP_t> returnTestCyl
  (thrust::host_vector<cudaTextureObject_t>& tex,
   int m, int n, int k, int nmax, int NUMX, int NUMY)
{
  thrust::device_vector<cudaTextureObject_t> t_d = tex;
  
  unsigned int        N  = NUMX*NUMY;
  unsigned int gridSize  = N/BLOCK_SIZE;
  if (N > gridSize*BLOCK_SIZE) gridSize++;

  thrust::device_vector<cuFP_t> f_d(N);

  testFetchCyl<<<gridSize, BLOCK_SIZE>>>(toKernel(t_d), toKernel(f_d),
					 m, n, k, nmax, NUMX, NUMY);

  cudaDeviceSynchronize();

  return f_d;
}

void EmpCylSL::initialize_cuda
(std::vector<cudaArray_t>& cuArray,
 thrust::host_vector<cudaTextureObject_t>& tex
 )
{
  // Number of texture arrays
  //
  size_t ndim = (MMAX+1)*rank3;

  // Interpolation data array
  //
  cuArray.resize(ndim);

  // Create texture objects
  //
  tex.resize(ndim);
  thrust::fill(tex.begin(), tex.end(), 0);

  cudaTextureDesc texDesc;

  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.readMode = cudaReadModeElementType;
  texDesc.filterMode = cudaFilterModePoint;
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.addressMode[1] = cudaAddressModeClamp;
  texDesc.addressMode[2] = cudaAddressModeClamp;
  
  // Temporary storage
  //
  cuFP_t *d_Interp;
  cuda_safe_call(cudaMalloc((void **)&d_Interp, NUMX*NUMY*6*sizeof(cuFP_t)),
		 __FILE__, __LINE__,
		 "Error allocating d_Interp for texture construction");
  
  std::vector<cuFP_t> h_buffer(NUMX*NUMY*6, 0.0);
  size_t k = 0;

  for (size_t mm=0; mm<=MMAX; mm++) {

    for (size_t n=0; n<rank3; n++) {

      // Copy table to flat array
      //
      for (int j=0; j<NUMY; j++) {
	for (int i=0; i<NUMX; i++) {
	  h_buffer[i+j*NUMX              ]   = potC   [mm][n](i, j);
	  h_buffer[i+j*NUMX + NUMX*NUMY  ]   = rforceC[mm][n](i, j);
	  h_buffer[i+j*NUMX + NUMX*NUMY*2]   = zforceC[mm][n](i, j);
	  if (mm) {
	    h_buffer[i+j*NUMX + NUMX*NUMY*3] = potS   [mm][n](i, j);
	    h_buffer[i+j*NUMX + NUMX*NUMY*4] = rforceS[mm][n](i, j);
	    h_buffer[i+j*NUMX + NUMX*NUMY*5] = zforceS[mm][n](i, j);
	  }
	}
      }
      
      // Copy data to device
      cuda_safe_call(cudaMemcpy(d_Interp, &h_buffer[0], NUMX*NUMY*6*sizeof(cuFP_t), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying texture table to device");

      // cudaArray Descriptor
      //
#if cuREAL == 4
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
#else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
#endif
      // cuda Array
      //
      cuda_safe_call(cudaMalloc3DArray(&cuArray[k], &channelDesc, make_cudaExtent(NUMX, NUMY, 6), 0), __FILE__, __LINE__, "Error allocating cuArray for 3d texture");

      // Array creation
      //
      cudaMemcpy3DParms copyParams = {0};
      
      copyParams.srcPtr   = make_cudaPitchedPtr(d_Interp, NUMX*sizeof(cuFP_t), NUMX, NUMY);
      copyParams.dstArray = cuArray[k];
      copyParams.extent   = make_cudaExtent(NUMX, NUMY, 6);
      copyParams.kind     = cudaMemcpyDeviceToDevice;

      cuda_safe_call(cudaMemcpy3D(&copyParams), __FILE__, __LINE__, "Error in copying 3d pitched array");

      cudaResourceDesc resDesc;

      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array  = cuArray[k];

      cuda_safe_call
	(cudaCreateTextureObject(&tex[k], &resDesc, &texDesc, NULL),
	 __FILE__, __LINE__, "Failure in 2d texture creation");
      
      k++;

    } // radial order loop

  } // harmonic subspace loop


  assert(k == ndim);

  cuda_safe_call(cudaFree(d_Interp), __FILE__, __LINE__, "Failure freeing device memory");

  // This is for debugging: compare texture table fetches to original
  // tables
  //
  if (false) {
    constexpr cuFP_t tol = 10.0*std::numeric_limits<cuFP_t>::epsilon();

    struct Element {
      int m;
      int k;
      double a;
      double b;
    };

    std::multimap<double, Element> compare;

    thrust::host_vector<cuFP_t> xyg;
    std::cout << "**HOST** Texture 2D compare" << std::endl;
    unsigned tot = 0, bad = 0;
    for (int mm=0; mm<=MMAX; mm++) {
      for (size_t n=0; n<rank3; n++) {
	
	std::vector<Eigen::MatrixXd*> orig =
	  {&potC[mm][n], &rforceC[mm][n], &zforceC[mm][n],
	   &potS[mm][n], &rforceS[mm][n], &zforceS[mm][n]};

	int kmax = 3;
	if (mm) kmax += 3;

	for (int k=0; k<kmax; k++) {
	  xyg = returnTestCyl(tex, mm, n, k, rank3, NUMX, NUMY);
	
	  for (int j=0; j<NUMY; j++) {
	    for (int i=0; i<NUMX; i++) {
	      double a = (*orig[k])(i, j);
	      double b = xyg[j*NUMX + i];
	      if (a>1.0e-18) {

		Element comp = {mm, k, a, b};
		compare.insert(std::make_pair(fabs((a - b)/a), comp));

		if ( fabs((a - b)/a ) > tol) {
		  std::cout << std::setw( 5) << mm << std::setw( 5) << n
			    << std::setw( 5) << i  << std::setw( 5) << j
			    << std::setw( 5) << k  << std::setw(20) << a
			    << std::setw(20) << b  << std::setw(20) << (a-b)/a
			    << std::endl;
		  bad++;
		}
	      }
	      tot++;
	    }
	  }
	}
      }
    }

    std::multimap<double, Element>::iterator beg = compare.begin();
    std::multimap<double, Element>::iterator end = compare.end();
    std::multimap<double, Element>::iterator lo1=beg, lo9=beg, mid=beg, hi9=end, hi1=end;

    std::advance(lo9, 9);
    std::advance(mid, compare.size()/2);
    std::advance(hi1, -1);
    std::advance(hi9, -10);

    std::cout << "**Found " << bad << "/" << tot << " values > eps" << std::endl
	      << "**Low[1] : " << lo1->first << " (" << lo1->second.m << ", " << lo1->second.k << ", " << lo1->second.a << ", " << lo1->second.b << ", " << lo1->second.a - lo1->second.b << ")" << std::endl
	      << "**Low[9] : " << lo9->first << " (" << lo9->second.m << ", " << lo9->second.k << ", " << lo9->second.a << ", " << lo9->second.b << ", " << lo9->second.a - lo9->second.b << ")" << std::endl
	      << "**Middle : " << mid->first << " (" << mid->second.m << ", " << mid->second.k << ", " << mid->second.a << ", " << mid->second.b << ", " << mid->second.a - mid->second.b << ")" << std::endl
	      << "**Hi [9] : " << hi9->first << " (" << hi9->second.m << ", " << hi9->second.k << ", " << hi9->second.a << ", " << lo1->second.b << ", " << hi9->second.a - hi9->second.b << ")" << std::endl
	      << "**Hi [1] : " << hi1->first << " (" << hi1->second.m << ", " << hi1->second.k << ", " << hi1->second.a << ", " << hi1->second.b << ", " << hi1->second.a - hi1->second.b << ")" << std::endl
	      << "**" << std::endl;
  }
}


cudaMappingConstants EmpCylSL::getCudaMappingConstants()
{
  cudaMappingConstants ret;

  ret.rscale = ASCALE;
  ret.hscale = HSCALE;
  ret.xmin   = XMIN;
  ret.xmax   = XMAX;
  ret.ymin   = YMIN;
  ret.ymax   = YMAX;
  ret.numr   = 0;
  ret.numx   = NUMX;
  ret.numy   = NUMY;
  ret.dxi    = dX;
  ret.dyi    = dY;
  ret.cmapR  = CMAPR;
  ret.cmapZ  = CMAPZ;

  return ret;
}
