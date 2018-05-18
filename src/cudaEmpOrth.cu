#include <EmpOrth9thd.h>

#include <iostream>
#include <iomanip>
#include <map>

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>

void EmpCylSL::initialize_cuda
(std::vector<float*>& cuArray,
 std::vector<cudaResourceDesc>& resDesc,
 struct cudaTextureDesc& texDesc,
 thrust::host_vector<cudaTextureObject_t>& tex
 )
{
  // Number of texture arrays
  //
  size_t ndim = 3*(2*MMAX+1)*rank3;

  // Interpolation data array
  //
  cuArray.resize(ndim);

  // Create texture objects
  //
  tex.resize(ndim);
  thrust::fill(tex.begin(), tex.end(), 0);

  resDesc.resize(ndim);

  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.readMode = cudaReadModeElementType;
  texDesc.filterMode = cudaFilterModePoint;
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.addressMode[1] = cudaAddressModeClamp;
  
  // Temporary storage
  //
  std::vector<float> h_buffer(NUMX*NUMY);
  size_t k = 0, pitch;

  for (size_t mm=0; mm<=MMAX; mm++) {

    for (size_t n=0; n<rank3; n++) {
    
      size_t  ntab = 3;		// m=0
      if (mm) ntab = 6;		// m>0

      for (size_t t=0; t<ntab; t++) {
	cuda_safe_call
	  (cudaMallocPitch(&cuArray[k], &pitch, sizeof(float)*NUMX, NUMY),
	   __FILE__, __LINE__, "malloc pitch");
  
	cuda_safe_call
	  (cudaMemset2D(cuArray[k], pitch, 0, pitch, NUMY),
	   __FILE__, __LINE__, "memset 2d");
	
	// Copy table to flat array
	//
	for (int j=0; j<NUMY; j++) {
	  for (int i=0; i<NUMX; i++) {
	    if (t==0) h_buffer[i+j*NUMX] = potC   [mm][n][i][j];
	    if (t==1) h_buffer[i+j*NUMX] = rforceC[mm][n][i][j];
	    if (t==2) h_buffer[i+j*NUMX] = zforceC[mm][n][i][j];
	    if (t==3) h_buffer[i+j*NUMX] = potS   [mm][n][i][j];
	    if (t==4) h_buffer[i+j*NUMX] = rforceS[mm][n][i][j];
	    if (t==5) h_buffer[i+j*NUMX] = zforceS[mm][n][i][j];
	  }
	}
      
	cuda_safe_call
	  (cudaMemcpy2D(cuArray[k], pitch, &h_buffer[0], sizeof(float)*NUMX, sizeof(float)*NUMX, NUMY, cudaMemcpyHostToDevice),
	   __FILE__, __LINE__, "memcpy 2d");
	
	memset(&resDesc[k], 0, sizeof(resDesc[k]));
	resDesc[k].resType                   = cudaResourceTypePitch2D;
	resDesc[k].res.pitch2D.devPtr        = cuArray[k];
	resDesc[k].res.pitch2D.pitchInBytes  = pitch;
	resDesc[k].res.pitch2D.width         = NUMX;
	resDesc[k].res.pitch2D.height        = NUMY;
	resDesc[k].res.pitch2D.desc          = cudaCreateChannelDesc<float>();
	
	cuda_safe_call
	  (cudaCreateTextureObject(&tex[k], &resDesc[k], &texDesc, NULL),
	   __FILE__, __LINE__, "2d texture creation");
	
	k++;

      } // table loop

    } // radial order loop

  } // harmonic subspace loop


  assert(k == ndim);

  /*
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
  */
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
