#include <EmpOrth9thd.h>

#include <iostream>
#include <iomanip>
#include <map>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>

__global__
void testFetchCyl(cudaTextureObject_t* T, double* f, int m, int n, int i, int j, int k, int nmax)
{
  int l = m*nmax + n;
  *f = int2_as_double(tex3D<int2>(T[l], i, j, k));
}

double returnTestCyl(thrust::host_vector<cudaTextureObject_t>& tex, int m, int n, int i, int j, int k, int nmax)
{
  thrust::device_vector<cudaTextureObject_t> t_d = tex;
  cudaTextureObject_t *T = thrust::raw_pointer_cast(t_d.data());
  
  double* f;
  cuda_safe_call(cudaMalloc(&f, sizeof(double)),  __FILE__, __LINE__, "cudaMalloc");

  testFetchCyl<<<1, 1>>>(T, f, m, n, i, j, k, nmax);

  cudaDeviceSynchronize();

  double ret;
  cuda_safe_call(cudaMemcpy(&ret, f, sizeof(double), cudaMemcpyDeviceToHost), __FILE__, __LINE__, "cudaMemcpy");

  cuda_safe_call(cudaFree(f), __FILE__, __LINE__, "cudaFree");

  return ret;
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
  int2 *d_Interp;
  cuda_safe_call(cudaMalloc((void **)&d_Interp, NUMX*NUMY*6*sizeof(int2)),
		 __FILE__, __LINE__,
		 "Error allocating d_Interp for texture construction");
  
  std::vector<double> h_buffer(NUMX*NUMY*6, 0.0);
  size_t k = 0;

  for (size_t mm=0; mm<=MMAX; mm++) {

    for (size_t n=0; n<rank3; n++) {

      // Copy table to flat array
      //
      for (int j=0; j<NUMY; j++) {
	for (int i=0; i<NUMX; i++) {
	  h_buffer[i+j*NUMX              ]   = potC   [mm][n][i][j];
	  h_buffer[i+j*NUMX + NUMX*NUMY  ]   = rforceC[mm][n][i][j];
	  h_buffer[i+j*NUMX + NUMX*NUMY*2]   = zforceC[mm][n][i][j];
	  if (mm) {
	    h_buffer[i+j*NUMX + NUMX*NUMY*3] = potS   [mm][n][i][j];
	    h_buffer[i+j*NUMX + NUMX*NUMY*4] = rforceS[mm][n][i][j];
	    h_buffer[i+j*NUMX + NUMX*NUMY*5] = zforceS[mm][n][i][j];
	  }
	}
      }
      
      // Copy data to device
      cuda_safe_call(cudaMemcpy(d_Interp, &h_buffer[0], NUMX*NUMY*6*sizeof(int2), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying texture table to device");

      // cudaArray Descriptor
      //
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();

      // cuda Array
      //
      cuda_safe_call(cudaMalloc3DArray(&cuArray[k], &channelDesc, make_cudaExtent(NUMX, NUMY, 6), 0), __FILE__, __LINE__, "Error allocating cuArray for 3d texture");

      // Array creation
      //
      cudaMemcpy3DParms copyParams = {0};
      
      copyParams.srcPtr   = make_cudaPitchedPtr(d_Interp, NUMX*sizeof(int2), NUMX, NUMY);
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

  if (true) {
    std::cout << "**HOST** Texture 2D compare" << std::endl;
    unsigned tot = 0, bad = 0;
    double a, b;
    for (size_t mm=0; mm<=MMAX; mm++) {
      for (size_t n=0; n<rank3; n++) {
	for (int j=0; j<NUMY; j++) {
	  for (int i=0; i<NUMX; i++) {
	    {
	      double a = potC[mm][n][i][j];
	      double b = returnTestCyl(tex, mm, n, i, j, 0, rank3);
	      if (a>1.0e-18) {
		if ( fabs((a - b)/a ) > 1.0e-7) {
		  std::cout << std::setw( 5) << mm << std::setw( 5) << n
			    << std::setw( 5) << i  << std::setw( 5) << j
			    << std::setw(15) << a  << std::setw(15) << (a-b)/a << std::endl;
		  bad++;
		}
	      }
	      tot++;
	    }

	    {
	      double a = rforceC[mm][n][i][j];
	      double b = returnTestCyl(tex, mm, n, i, j, 1, rank3);
	      if (a>1.0e-18) {
		if ( fabs((a - b)/a ) > 1.0e-7) {
		  std::cout << std::setw( 5) << mm << std::setw( 5) << n
			    << std::setw( 5) << i  << std::setw( 5) << j
			    << std::setw(15) << a  << std::setw(15) << (a-b)/a << std::endl;
		  bad++;
		}
	      }
	      tot++;
	    }

	    {
	      double a = zforceC[mm][n][i][j];
	      double b = returnTestCyl(tex, mm, n, i, j, 2, rank3);
	      if (a>1.0e-18) {
		if ( fabs((a - b)/a ) > 1.0e-7) {
		  std::cout << std::setw( 5) << mm << std::setw( 5) << n
			    << std::setw( 5) << i  << std::setw( 5) << j
			    << std::setw(15) << a  << std::setw(15) << (a-b)/a << std::endl;
		  bad++;
		}
	      }
	      tot++;
	    }

	    if (mm) {
	      {
		double a = potS[mm][n][i][j];
		double b = returnTestCyl(tex, mm, n, i, j, 3, rank3);
		if (a>1.0e-18) {
		  if ( fabs((a - b)/a ) > 1.0e-7) {
		    std::cout << std::setw( 5) << mm << std::setw( 5) << n
			      << std::setw( 5) << i  << std::setw( 5) << j
			      << std::setw(15) << a  << std::setw(15) << (a-b)/a << std::endl;
		    bad++;
		  }
		}
		tot++;
	      }

	      {
		double a = rforceS[mm][n][i][j];
		double b = returnTestCyl(tex, mm, n, i, j, 4, rank3);
		if (a>1.0e-18) {
		  if ( fabs((a - b)/a ) > 1.0e-7) {
		    std::cout << std::setw( 5) << mm << std::setw( 5) << n
			      << std::setw( 5) << i  << std::setw( 5) << j
			      << std::setw(15) << a  << std::setw(15) << (a-b)/a << std::endl;
		    bad++;
		  }
		}
		tot++;
	      }

	      {
		double a = zforceS[mm][n][i][j];
		double b = returnTestCyl(tex, mm, n, i, j, 5, rank3);
		if (a>1.0e-18) {
		  if ( fabs((a - b)/a ) > 1.0e-7) {
		    std::cout << std::setw( 5) << mm << std::setw( 5) << n
			      << std::setw( 5) << i  << std::setw( 5) << j
			      << std::setw(15) << a  << std::setw(15) << (a-b)/a << std::endl;
		    bad++;
		  }
		}
		tot++;
	      }
	    }
	  }
	}
      }
    }

    std::cout << "**Found " << bad << "/" << tot << " bad values" << std::endl
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
  ret.cmap   = (CMAP ? 1 : 0);

  return ret;
}
