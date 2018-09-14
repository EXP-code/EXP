// -*- C++ -*-

#include <boost/make_shared.hpp>

/** CUDA strategy

(1) All cross sections with lookup tables to be copied to cuda texture objects

(2) Copy necessary constants to cuda symbol data

(3) Cross sections are then evaluated on the device so we need cuda
__device__ functions for each desired cross section

  Four grids
  ----------
  free free    2-dimensional
  rad recomb   1-dimensional
  col excite   1-dimensional (2 values)
  col ionize   1-dimensional

*/

// Define for debugging
//
// #define OFF_GRID_ALERT
// #define BOUNDS_CHECK
// #define VERBOSE

// Global symbols for coordinate transformation
//
__device__ __constant__
cuFP_t ionEminGrid, ionEmaxGrid, ionDeltaEGrid;

__device__ __constant__
cuFP_t ionCollideEmin, ionEmaxGrid, ionDeltaEGrid;

__global__
void testConstantsIon()
{
  printf("** Rscale = %f\n", cylRscale);
  printf("** Hscale = %f\n", cylHscale);
  printf("** Xmin   = %f\n", cylXmin);
  printf("** Xmax   = %f\n", cylXmax);
  printf("** Ymin   = %f\n", cylYmin);
  printf("** Ymax   = %f\n", cylYmax);
  printf("** Dxi    = %f\n", cylDxi);
  printf("** Dyi    = %f\n", cylDyi);
  printf("** Numx   = %d\n", cylNumx);
  printf("** Numy   = %d\n", cylNumy);
  printf("** Cmap   = %d\n", cylCmap);
}

void chdata::cuda_initialize()
{
  for (auto v : IonList) {
    // Indecies of ion state
    //
    cuZ.push_back(v.first.first );
    cuC.push_back(v.first.second);
    
    // Initialize Ion cross-section grids
    //
    if (not freeFreeGridComputed)
      v.second->freeFreeMakeEvGrid(0);

    if (v.first.second>1 and not radRecombGridComputed)
      radRecombMakeEvGrid(0);

    if (v.first.second[k]<=v.first.first[k]) {
      if (not exciteGridComputed) collExciteMakeGrid(id);
      if (not ionizeGridComputed) directIonMakeGrid(id);
    }
  }

}

void chdata:cuda_initialize_textures()
{
  size_t ionSize = cuZ.size();

  // Interpolation data array
  //
  cuFFarray.resize(ionSize);
  cuRCarray.resize(ionSize);
  cuCEarray.resize(ionSize);
  cuCIarray.resize(ionSize);

  std::vector<cudaArray_t> cuFFarray, cuRCarray, cuCEarray, ciCIarray;
  std::vector<cudaTextureObject_t> ff_d, rc_d, ce_d, ci_d;
  std::vector<int> cuZ, cuC;

  thrust::host_vector  <double> hEminArray, hEmaxArray;
  thrust::device_vector<double> dEminArray, dEmaxArray;


  // Texture array
  //
  ff_d.resize(ionSize);
  rc_d.resize(ionSize);
  ce_d.resize(ionSize);
  ci_d.resize(ionSize);

  thrust::fill(ff_d.begin(), ff_d.end(), 0);
  thrust::fill(rc_d.begin(), rc_d.end(), 0);
  thrust::fill(ce_d.begin(), ce_d.end(), 0);
  thrust::fill(ci_d.begin(), ci_d.end(), 0);

  size_t k = 0;

  for (auto v : IonList) {

    IonPtr I = v.second;

    // The free-free array
    if (true) {
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
      cuda_safe_call(cudaMalloc((void **)&d_Interp, I->NfreeFreeGrid*I->kffsteps*sizeof(cuFP_t)),
		   __FILE__, __LINE__,
		     "Error allocating d_Interp for texture construction");
  
      std::vector<cuFP_t> h_buffer(NfreeFreeGrid*kffsteps, 0.0);

      // Copy vectors to buffer
    //
      for (int i = 0; i < I->NfreeFreeGrid; i++) {
	for (int j = 0; j < kffsteps; j++) {
	  h_buffer[i + j*I->NfreeFreeGrid] = I->freeFreeGrid[i][j];
	}
      }
    
    // Copy data to device
    cuda_safe_call(cudaMemcpy(d_Interp, &h_buffer[0], I->NfreeFreeGrid*sizeof(cuFP_t), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying texture table to device");
    
    // cudaArray Descriptor
    //
#if cuREAL == 4
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
#else
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
#endif
    // cuda Array
    //
    cuda_safe_call(cudaMalloc3DArray(&cuFFarray[k], &channelDesc, make_cudaExtent(I->NfreeFreeGrid, I->kffsteps, 1), 0), __FILE__, __LINE__, "Error allocating cuArray for 3d texture");
    
    // Array creation
    //
    cudaMemcpy3DParms copyParams = {0};
  
    copyParams.srcPtr   = make_cudaPitchedPtr(d_Interp, NfreeFreeGrid*sizeof(cuFP_t), I->NfreeFreeGrid, I->kffsteps);
    copyParams.dstArray = cuFFarray[k];
    copyParams.extent   = make_cudaExtent(I->NfreeFreeGrid, I->kffsteps, 1);
    copyParams.kind     = cudaMemcpyDeviceToDevice;
    
    cuda_safe_call(cudaMemcpy3D(&copyParams), __FILE__, __LINE__, "Error in copying 3d pitched array");

    cudaResourceDesc resDesc;

    memset(&resDesc, 0, sizeof(cudaResourceDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array  = cuFFarray[k];
    
    cuda_safe_call
      (cudaCreateTextureObject(&ff_d[k], &resDesc, &texDesc, NULL),
       __FILE__, __LINE__, "Failure in 2d texture creation");
    
    cuda_safe_call(cudaFree(d_Interp), __FILE__, __LINE__, "Failure freeing device memory");
  }

  // Radiative recombination texture (1-d)
  //
  if (cuC[k]>1) {
    // Allocate CUDA array in device memory (a one-dimension 'channel')
    //
#if cuREAL == 4
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
#else
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
#endif
    
    // Size of interpolation array
    //
    size_t tsize = I->NradRecombGrid*sizeof(cuFP_t);

    cudaTextureDesc texDesc;

    memset(&texDesc, 0, sizeof(cudaTextureDesc));
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.readMode = cudaReadModeElementType;
    texDesc.normalizedCoords = 0;
    
    thrust::host_vector<cuFP_t> tt(I->NradRecombGrid);

    cuda_safe_call(cudaMallocArray(&cuRCarray[k], &channelDesc, I->NradRecombGrid), __FILE__, __LINE__, "malloc cuArray");

    // Copy to device memory some data located at address h_data
    // in host memory
    for (size_t n = 0; n < I->NradRecombGrid; n++) tt[n] = I->radRecombGrid[n];
    
    cuda_safe_call(cudaMemcpyToArray(cuRCarray[k], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

    // Specify texture

    cudaResourceDesc resDesc;

    memset(&resDesc, 0, sizeof(cudaResourceDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = cuRCarray[k];

    cuda_safe_call(cudaCreateTextureObject(&rc_d[k], &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");
  }

  // The collisional excitation array

  if (cuC[k] <= cuZ[k]) {

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
    cuda_safe_call(cudaMalloc((void **)&d_Interp, I->NcollideGrid*2*sizeof(cuFP_t)),
		   __FILE__, __LINE__,
		   "Error allocating d_Interp for texture construction");
  
    std::vector<cuFP_t> h_buffer(I->NcollideGrid*2, 0.0);

    // Copy vectors to buffer
    //
    for (int i = 0; i < I->NcollideGrid; i++) {
      h_buffer[i                ] = I->collideDataGrid[i].first;
      h_buffer[i + NfreeFreeGrid] = I->collideDataGrid[i].second;
    }
    
    // Copy data to device
    cuda_safe_call(cudaMemcpy(d_Interp, &h_buffer[0], I->NcollideGrid*sizeof(cuFP_t), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying texture table to device");
    
    // cudaArray Descriptor
    //
#if cuREAL == 4
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
#else
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
#endif
    // cuda Array
    //
    cuda_safe_call(cudaMalloc3DArray(&cuCEarray[k], &channelDesc, make_cudaExtent(I->NcollideGrid, 2, 1), 0), __FILE__, __LINE__, "Error allocating cuArray for 3d texture");
    
    // Array creation
    //
    cudaMemcpy3DParms copyParams = {0};
  
    copyParams.srcPtr   = make_cudaPitchedPtr(d_Interp, NcollideGrid*sizeof(cuFP_t), NcollideGrid, 2);
    copyParams.dstArray = cuInterpArray[2];
    copyParams.extent   = make_cudaExtent(I->NcollideGrid, 2, 1);
    copyParams.kind     = cudaMemcpyDeviceToDevice;
    
    cuda_safe_call(cudaMemcpy3D(&copyParams), __FILE__, __LINE__, "Error in copying 3d pitched array");

    cudaResourceDesc resDesc;

    memset(&resDesc, 0, sizeof(cudaResourceDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array  = cuCEarray[k];
    
    cuda_safe_call
      (cudaCreateTextureObject(&ce_d[k], &resDesc, &texDesc, NULL),
       __FILE__, __LINE__, "Failure in 2d texture creation");
    
    cuda_safe_call(cudaFree(d_Interp), __FILE__, __LINE__, "Failure freeing device memory");
  }

  if (cuZ[k] <= cuC[k]) {

    // Allocate CUDA array in device memory (a one-dimension 'channel')
    //
#if cuREAL == 4
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
#else
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
#endif
    
    // Size of interpolation array
    //
    size_t tsize = I->NionizeGrid*sizeof(cuFP_t);

    cudaTextureDesc texDesc;

    memset(&texDesc, 0, sizeof(cudaTextureDesc));
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.readMode = cudaReadModeElementType;
    texDesc.normalizedCoords = 0;
    
    thrust::host_vector<cuFP_t> tt(NionizeGrid);

    cuda_safe_call(cudaMallocArray(&cuCIarray[k], &channelDesc, I->NionizeGrid), __FILE__, __LINE__, "malloc cuArray");

    // Copy to device memory some data located at address h_data
    // in host memory
    for (size_t n = 0; n < I->NionizeGrid; n++) tt[n] = I->ionizeDataGrid[n];
    
    cuda_safe_call(cudaMemcpyToArray(cuCIarray[k], 0, 0, &tt[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

    // Specify texture

    cudaResourceDesc resDesc;

    memset(&resDesc, 0, sizeof(cudaResourceDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = cuCIarray[k];

    cuda_safe_call(cudaCreateTextureObject(&ci_d[k], &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");
  }

  // Grid values for the various cross sections
  /*
  std::vector< double > kgrid;
  std::vector< double > kgr10;
  std::vector< double > nugrid;

  int kffsteps;
  int nusteps;
  */
}

void Ion::initialize_grid_constants()
{
  // Copy constants to device
  //
  
  cuda_safe_call(cudaMemcpyToSymbol(cylRscale, &f.rscale, sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylRscale");

  cuda_safe_call(cudaMemcpyToSymbol(cylHscale, &f.hscale, sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylHscale");

  cuda_safe_call(cudaMemcpyToSymbol(cylXmin,   &f.xmin,   sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylXmin");

  cuda_safe_call(cudaMemcpyToSymbol(cylXmax,   &f.xmax,   sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylXmax");

  cuda_safe_call(cudaMemcpyToSymbol(cylDxi,    &f.dxi,    sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylDxi");

  cuda_safe_call(cudaMemcpyToSymbol(cylNumx,   &f.numx,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylNumx");

  cuda_safe_call(cudaMemcpyToSymbol(cylYmin,   &f.ymin,   sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylYmin");

  cuda_safe_call(cudaMemcpyToSymbol(cylYmax,   &f.ymax,   sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylYmax");

  cuda_safe_call(cudaMemcpyToSymbol(cylDyi,    &f.dyi,    sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylDxi");

  cuda_safe_call(cudaMemcpyToSymbol(cylNumy,   &f.numy,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylNumy");

  cuda_safe_call(cudaMemcpyToSymbol(cylCmap,   &f.cmap,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylCmap");

}


__global__ void coordKernelCyl
(dArray<cudaParticle> in, dArray<cuFP_t> mass, dArray<cuFP_t> phi,
 dArray<cuFP_t> Xfac, dArray<cuFP_t> Yfac,
 dArray<int> IndX, dArray<int> IndY,
 unsigned int stride, PII lohi, cuFP_t rmax)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i = tid*stride + n;	// Particle counter
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {

#ifdef BOUNDS_CHECK
      if (npart>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle p = in._v[npart];
    
      cuFP_t xx=0.0, yy=0.0, zz=0.0;

      if (cylOrient) {
	for (int k=0; k<3; k++) xx += cylBody[0+k]*(p.pos[k] - cylCen[k]);
	for (int k=0; k<3; k++) yy += cylBody[3+k]*(p.pos[k] - cylCen[k]);
	for (int k=0; k<3; k++) zz += cylBody[6+k]*(p.pos[k] - cylCen[k]);
      } else {
	xx = p.pos[0] - cylCen[0];
	yy = p.pos[1] - cylCen[1];
	zz = p.pos[2] - cylCen[2];
      }
      
      cuFP_t R2 = xx*xx + yy*yy;
      cuFP_t r2 = R2 + zz*zz;
      cuFP_t R  = sqrt(R2);
      cuFP_t r  = sqrt(r2);
#ifdef BOUNDS_CHECK
      if (i>=mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      mass._v[i] = -1.0;
      
      if (r<=rmax) {
	
	mass._v[i] = p.mass;
	
	phi._v[i] = atan2(yy, xx);

#ifdef BOUNDS_CHECK
	if (i>=phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	// Interpolation indices
	//
	cuFP_t X  = (cu_r_to_xi_cyl(R) - cylXmin)/cylDxi;
	cuFP_t Y  = (cu_z_to_y_cyl(zz) - cylYmin)/cylDyi;

	int indX = floor(X);
	int indY = floor(Y);
	
	if (indX<0) indX = 0;
	if (indX>cylNumx-2) indX = cylNumx - 2;
	
	if (indY<0) indY = 0;
	if (indY>cylNumy-2) indY = cylNumy - 2;
	
	Xfac._v[i] = cuFP_t(indX+1) - X;
	IndX._v[i] = indX;

	Yfac._v[i] = cuFP_t(indY+1) - Y;
	IndY._v[i] = indY;

#ifdef OFF_GRID_ALERT
	if (Xfac._v[i]<-0.5 or Xfac._v[i]>1.5) printf("X off grid: x=%f\n", X);
	if (Yfac._v[i]<-0.5 or Yfac._v[i]>1.5) printf("Y off grid: y=%f\n", Y);
#endif
#ifdef BOUNDS_CHECK
	if (i>=Xfac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=IndX._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=Yfac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=IndY._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      }
    }
  }
}


__global__ void coefKernelCyl
(dArray<cuFP_t> coef, dArray<cudaTextureObject_t> tex,
 dArray<cuFP_t> Mass, dArray<cuFP_t> Phi,
 dArray<cuFP_t> Xfac, dArray<cuFP_t> Yfac,
 dArray<int> indX, dArray<int> indY,
 int stride, int m, unsigned int nmax, PII lohi)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of particles to be evaluated
  //
  const unsigned int N = lohi.second - lohi.first;

  const cuFP_t norm = -4.0*M_PI;	// Biorthogonality factor

  for (int istr=0; istr<stride; istr++) {

    int i = tid*stride + istr;	// Particle counter

    if (i<N) {			// Allow for grid padding

      cuFP_t mass = Mass._v[i];
      
#ifdef BOUNDS_CHECK
      if (i>=Mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif      
      if (mass>0.0) {
#ifdef BOUNDS_CHECK
	if (i>=Phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	cuFP_t phi  = Phi._v[i];
	cuFP_t cosp = cos(phi*m);
	cuFP_t sinp = sin(phi*m);
	
	// Do the interpolation
	//
	cuFP_t delx0 = Xfac._v[i];
	cuFP_t dely0 = Yfac._v[i];
	cuFP_t delx1 = 1.0 - delx0;
	cuFP_t dely1 = 1.0 - dely0;

#ifdef BOUNDS_CHECK
	if (i>=Xfac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=Yfac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	cuFP_t c00 = delx0*dely0;
	cuFP_t c10 = delx1*dely0;
	cuFP_t c01 = delx0*dely1;
	cuFP_t c11 = delx1*dely1;

	int   indx = indX._v[i];
	int   indy = indY._v[i];

#ifdef BOUNDS_CHECK
	if (i>=indX._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=indY._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	for (int n=0; n<nmax; n++) {

	  // Texture maps are packed in slices
	  // ---------------------------------
	  // potC, rforceC, zforceC, potS, rforceS, zforceS
	  // 0     1        2        3     4        5

	  int k = m*nmax + n;

#if cuREAL == 4
	  cuFP_t d00  = tex3D<float>(tex._v[k], indx,   indy  , 0);
	  cuFP_t d10  = tex3D<float>(tex._v[k], indx+1, indy  , 0);
	  cuFP_t d01  = tex3D<float>(tex._v[k], indx,   indy+1, 0);
	  cuFP_t d11  = tex3D<float>(tex._v[k], indx+1, indy+1, 0);

#else
	  cuFP_t d00  = int2_as_double(tex3D<int2>(tex._v[k], indx,   indy  , 0));
	  cuFP_t d10  = int2_as_double(tex3D<int2>(tex._v[k], indx+1, indy  , 0));
	  cuFP_t d01  = int2_as_double(tex3D<int2>(tex._v[k], indx,   indy+1, 0));
	  cuFP_t d11  = int2_as_double(tex3D<int2>(tex._v[k], indx+1, indy+1, 0));
#endif

#ifdef BOUNDS_CHECK
	  if (k>=tex._s)            printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	  if ((2*n+0)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	  coef._v[(2*n+0)*N + i] = (c00*d00 + c10*d10 + c01*d01 + c11*d11) * cosp * norm * mass;

	  if (m>0) {
	    // potS tables are offset from potC tables by +3
	    //
#if cuREAL == 4
	    d00  = tex3D<float>(tex._v[k], indx,   indy  , 3);
	    d10  = tex3D<float>(tex._v[k], indx+1, indy  , 3);
	    d01  = tex3D<float>(tex._v[k], indx,   indy+1, 3);
	    d11  = tex3D<float>(tex._v[k], indx+1, indy+1, 3);
#else
	    d00  = int2_as_double(tex3D<int2>(tex._v[k], indx,   indy  , 3));
	    d10  = int2_as_double(tex3D<int2>(tex._v[k], indx+1, indy  , 3));
	    d01  = int2_as_double(tex3D<int2>(tex._v[k], indx,   indy+1, 3));
	    d11  = int2_as_double(tex3D<int2>(tex._v[k], indx+1, indy+1, 3));
#endif

	    coef._v[(2*n+1)*N + i] = (c00*d00 + c10*d10 + c01*d01 + c11*d11) * sinp * norm * mass;

#ifdef BOUNDS_CHECK
	    if ((2*n+1)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	  } // m>0
	  else {
	    coef._v[(2*n+1)*N + i] = 0.0;
	  }

	} // norder loop

      } else {
	// No contribution from off-grid particles
	for (int n=0; n<nmax; n++) {
	  coef._v[(2*n+0)*N + i] = 0.0;
	  if (m) coef._v[(2*n+1)*N + i] = 0.0;
	}

      } // mass value check

    } // particle index check

  } // stride loop
}

__global__ void
forceKernelCyl(dArray<cudaParticle> in, dArray<cuFP_t> coef,
	       dArray<cudaTextureObject_t> tex,
	       int stride, unsigned int mmax, unsigned int nmax, PII lohi,
	       cuFP_t rmax, cuFP_t cylmass, bool external)
{
  // Thread ID
  //
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

  // Maximum radius squared
  //
  const cuFP_t rmax2 = rmax*rmax;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {
      
#ifdef BOUNDS_CHECK
      if (npart>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle p = in._v[npart];
      
      cuFP_t acc[3] = {0.0, 0.0, 0.0};
      cuFP_t xx=0.0, yy=0.0, zz=0.0;

      if (cylOrient) {
	for (int k=0; k<3; k++) xx += cylBody[0+k]*(p.pos[k] - cylCen[k]);
	for (int k=0; k<3; k++) yy += cylBody[3+k]*(p.pos[k] - cylCen[k]);
	for (int k=0; k<3; k++) zz += cylBody[6+k]*(p.pos[k] - cylCen[k]);
      } else {
	xx = p.pos[0] - cylCen[0];
	yy = p.pos[1] - cylCen[1];
	zz = p.pos[2] - cylCen[2];
      }

      cuFP_t phi = atan2(yy, xx);
      cuFP_t R2  = xx*xx + yy*yy;
      cuFP_t  R  = sqrt(R2) + FSMALL;
      
      const cuFP_t ratmin = 0.75;
      const cuFP_t maxerf = 3.0;
      const cuFP_t midpt  = ratmin + 0.5*(1.0 - ratmin);
      const cuFP_t rsmth  = 0.5*(1.0 - ratmin)/maxerf;

      cuFP_t ratio = sqrt( (R2 + zz*zz)/rmax2 );
      cuFP_t mfactor = 1.0, frac = 1.0, cfrac = 0.0;

      if (ratio >= 1.0) {
	// cfrac      = 1.0 - mfactor;
	cfrac = 1.0;
      } else if (ratio > ratmin) {
	frac  = 0.5*(1.0 - erf( (ratio - midpt)/rsmth )) * mfactor;
	cfrac = 1.0 - frac;
      } else {
	frac  = mfactor;
      }

      cuFP_t fr = 0.0;
      cuFP_t fz = 0.0;
      cuFP_t fp = 0.0;
      cuFP_t pp = 0.0;
      
      if (ratio < 1.0) {

	cuFP_t X  = (cu_r_to_xi_cyl(R) - cylXmin)/cylDxi;
	cuFP_t Y  = (cu_z_to_y_cyl(zz) - cylYmin)/cylDyi;

	int indX = floor(X);
	int indY = floor(Y);
	
	if (indX < 0) indX = 0;
	if (indY < 0) indY = 0;
	if (indX >= cylNumx) indX = cylNumx - 1;
	if (indY >= cylNumy) indY = cylNumy - 1;

	cuFP_t delx0 = cuFP_t(indX+1) - X;
	cuFP_t dely0 = cuFP_t(indY+1) - Y;

#ifdef OFF_GRID_ALERT
	if (delx0<0.0 or delx0>1.0) printf("X off grid: x=%f\n", delx0);
	if (dely0<0.0 or dely0>1.0) printf("Y off grid: y=%f\n", dely0);
#endif

	cuFP_t delx1 = 1.0 - delx0;
	cuFP_t dely1 = 1.0 - dely0;
      
	cuFP_t c00 = delx0*dely0;
	cuFP_t c10 = delx1*dely0;
	cuFP_t c01 = delx0*dely1;
	cuFP_t c11 = delx1*dely1;

	cuFP_t cos1 = cos(phi);
	cuFP_t sin1 = sin(phi);

	cuFP_t ccos = 1.0;
	cuFP_t ssin = 0.0;

	for (int mm=0; mm<=mmax; mm++) {

	  for (int n=0; n<nmax; n++) {
      
	    cuFP_t fac0 = coef._v[Imn(mm, 'c', n, nmax)];
	    cuFP_t fac1 = fac0 * ccos;
	    cuFP_t fac2 = fac0 * ssin;
      
	    // Texture table index
	    //
	    int k = mm*nmax + n;

	    pp += fac1 *
	      (
#if cuREAL == 4
	       tex3D<float>(tex._v[k], indX,   indY  , 0) * c00 +
	       tex3D<float>(tex._v[k], indX+1, indY  , 0) * c10 +
	       tex3D<float>(tex._v[k], indX,   indY+1, 0) * c01 +
	       tex3D<float>(tex._v[k], indX+1, indY+1, 0) * c11 
#else
	       int2_as_double(tex3D<int2>(tex._v[k], indX,   indY  , 0)) * c00 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY  , 0)) * c10 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX,   indY+1, 0)) * c01 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY+1, 0)) * c11 
#endif
	       );
	    
	    fr += fac1 *
	      (
#if cuREAL == 4
	       tex3D<float>(tex._v[k], indX,   indY  , 1) * c00 +
	       tex3D<float>(tex._v[k], indX+1, indY  , 1) * c10 +
	       tex3D<float>(tex._v[k], indX,   indY+1, 1) * c01 +
	       tex3D<float>(tex._v[k], indX+1, indY+1, 1) * c11 
#else
	       int2_as_double(tex3D<int2>(tex._v[k], indX,   indY  , 1)) * c00 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY  , 1)) * c10 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX,   indY+1, 1)) * c01 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY+1, 1)) * c11 
#endif
	       );
      
	    fz += fac1 *
	      (
#if cuREAL == 4
	       tex3D<float>(tex._v[k], indX,   indY  , 2) * c00 +
	       tex3D<float>(tex._v[k], indX+1, indY  , 2) * c10 +
	       tex3D<float>(tex._v[k], indX,   indY+1, 2) * c01 +
	       tex3D<float>(tex._v[k], indX+1, indY+1, 2) * c11 
#else
	       int2_as_double(tex3D<int2>(tex._v[k], indX,   indY  , 2)) * c00 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY  , 2)) * c10 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX,   indY+1, 2)) * c01 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY+1, 2)) * c11 
#endif
	       );
	    
	    fp += fac2 * mm *
	      (
#if cuREAL == 4
	       tex3D<float>(tex._v[k], indX,   indY  , 0) * c00 +
	       tex3D<float>(tex._v[k], indX+1, indY  , 0) * c10 +
	       tex3D<float>(tex._v[k], indX,   indY+1, 0) * c01 +
	       tex3D<float>(tex._v[k], indX+1, indY+1, 0) * c11 
#else
	       int2_as_double(tex3D<int2>(tex._v[k], indX,   indY  , 0)) * c00 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY  , 0)) * c10 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX,   indY+1, 0)) * c01 +
	       int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY+1, 0)) * c11 
#endif
	       );
      
      
	    if (mm) {
	
	      cuFP_t fac0 =  coef._v[Imn(mm, 's', n, nmax)];
	      cuFP_t fac1 =  fac0 * ssin;
	      cuFP_t fac2 = -fac0 * ccos;

	      pp += fac1 *
		(
#if cuREAL == 4
		 tex3D<float>(tex._v[k], indX,   indY  , 3) * c00 +
		 tex3D<float>(tex._v[k], indX+1, indY  , 3) * c10 +
		 tex3D<float>(tex._v[k], indX,   indY+1, 3) * c01 +
		 tex3D<float>(tex._v[k], indX+1, indY+1, 3) * c11
#else		 
 		 int2_as_double(tex3D<int2>(tex._v[k], indX,   indY  , 3)) * c00 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY  , 3)) * c10 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX,   indY+1, 3)) * c01 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY+1, 3)) * c11 
#endif
		 );
	      
	      fr += fac1 *
		(
#if cuREAL == 4
		 tex3D<float>(tex._v[k], indX,   indY  , 4) * c00 +
		 tex3D<float>(tex._v[k], indX+1, indY  , 4) * c10 +
		 tex3D<float>(tex._v[k], indX,   indY+1, 4) * c01 +
		 tex3D<float>(tex._v[k], indX+1, indY+1, 4) * c11 
#else
		 int2_as_double(tex3D<int2>(tex._v[k], indX,   indY  , 4)) * c00 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY  , 4)) * c10 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX,   indY+1, 4)) * c01 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY+1, 4)) * c11 
#endif
		 );
	      
	      fz += fac1 *
		(
#if cuREAL == 4
		 tex3D<float>(tex._v[k], indX,   indY  , 5) * c00 +
		 tex3D<float>(tex._v[k], indX+1, indY  , 5) * c10 +
		 tex3D<float>(tex._v[k], indX,   indY+1, 5) * c01 +
		 tex3D<float>(tex._v[k], indX+1, indY+1, 5) * c11 
#else
		 int2_as_double(tex3D<int2>(tex._v[k], indX,   indY  , 5)) * c00 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY  , 5)) * c10 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX,   indY+1, 5)) * c01 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY+1, 5)) * c11 
#endif
		 );
	      
	      fp += fac2 * mm *
		(
#if cuREAL == 4
		 tex3D<float>(tex._v[k], indX,   indY  , 3) * c00 +
		 tex3D<float>(tex._v[k], indX+1, indY  , 3) * c10 +
		 tex3D<float>(tex._v[k], indX,   indY+1, 3) * c01 +
		 tex3D<float>(tex._v[k], indX+1, indY+1, 3) * c11 
#else
		 int2_as_double(tex3D<int2>(tex._v[k], indX,   indY  , 3)) * c00 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY  , 3)) * c10 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX,   indY+1, 3)) * c01 +
		 int2_as_double(tex3D<int2>(tex._v[k], indX+1, indY+1, 3)) * c11 
#endif
		 );
	      
	    }
	  }
	  
	  // Trig recursion to squeeze avoid internal FP fct call
	  //
	  cuFP_t cosM = ccos;
	  cuFP_t sinM = ssin;

	  ccos = cosM * cos1 - sinM * sin1;
	  ssin = sinM * cos1 + cosM * sin1;
	}

	acc[0] += ( fr*xx/R - fp*yy/R2 ) * frac;
	acc[1] += ( fr*yy/R + fp*xx/R2 ) * frac;
	acc[2] += fz * frac;
      }

      if (ratio > ratmin) {

	cuFP_t r3 = R2 + zz*zz;
	pp = -cylmass/sqrt(r3);	// -M/r
	fr = pp/r3;		// -M/r^3

	acc[0] += xx*fr * cfrac;
	acc[1] += yy*fr * cfrac;
	acc[2] += zz*fr * cfrac;
      }

      if (cylOrient) {
	for (int j=0; j<3; j++) {
	  for (int k=0; k<3; k++) in._v[npart].acc[j] += cylOrig[3*j+k]*acc[k];
	}
      } else {
	for (int j=0; j<3; j++) in._v[npart].acc[j] += acc[j];
      }

      if (external)
	in._v[npart].potext += pp;
      else
	in._v[npart].pot    += pp;

    } // Particle index block

  } // END: stride loop

}



template<typename T>
class LessAbs : public std::binary_function<bool, T, T>
{
public:
  T operator()( const T &a, const T &b ) const
  {
    return (fabs(a) < fabs(b));
  }
};


void Cylinder::cudaStorage::resize_coefs
(int ncylorder, int mmax, int N, int gridSize)
{
  // Reserve space for coefficient reduction
  //
  if (dN_coef.capacity() < 2*ncylorder*N)
    dN_coef.reserve(2*ncylorder*N);
  
  if (dc_coef.capacity() < 2*ncylorder*gridSize)
    dc_coef.reserve(2*ncylorder*gridSize);
  
  if (m_d .capacity() < N) m_d .reserve(N);
  if (X_d .capacity() < N) X_d .reserve(N);
  if (Y_d .capacity() < N) Y_d .reserve(N);
  if (p_d .capacity() < N) p_d .reserve(N);
  if (iX_d.capacity() < N) iX_d.reserve(N);
  if (iY_d.capacity() < N) iY_d.reserve(N);
  
  
  // Set space for current step
  //
  dN_coef.resize(2*ncylorder*N);
  dc_coef.resize(2*ncylorder*gridSize);
  dw_coef.resize(2*ncylorder);	// This will stay fixed

  // Space for coordinate arrays on the current step
  //
  m_d .resize(N);
  X_d .resize(N);
  Y_d .resize(N);
  p_d .resize(N);
  iX_d.resize(N);
  iY_d.resize(N);
}

void Cylinder::zero_coefs()
{
  Component::cuRingType cr = *cC->cuRing.get();
  Cylinder ::cuRingType ar = *cuRing.get();
  
  for (int n=0; n<cuStreams; n++) {
    // Resize output array
    //
    ar->df_coef.resize((mmax+1)*2*ncylorder);
    
    // Zero output array
    //
    thrust::fill(thrust::cuda::par.on(cr->stream),
		 ar->df_coef.begin(), ar->df_coef.end(), 0.0);

    // Advance iterators
    //
    ++cr;			// Component stream
    ++ar;			// Force method storage
  }
}

static bool initialize_cuda_cyl = true;
static unsigned dbg_id = 0;

void Cylinder::determine_coefficients_cuda()
{
  if (initialize_cuda_cyl) {
    initialize_cuda();
    initialize_mapping_constants();
    initialize_cuda_cyl = false;
    // Only copy texture memory once
    t_d = tex;
  }

  std::cout << std::scientific;

  int deviceCount = 0;
  cuda_safe_call(cudaGetDeviceCount(&deviceCount),
		 __FILE__, __LINE__, "could not get device count");

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, deviceCount-1);

  // This will stay fixed for the entire run
  //
  host_coefs.resize((2*mmax+1)*ncylorder);

  // Set component center and orientation
  //
  std::vector<cuFP_t> ctr;
  for (auto v : cC->getCenter(Component::Local | Component::Centered)) ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(cylCen, &ctr[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylCen");

  bool orient = (cC->EJ & Orient::AXIS) && !cC->EJdryrun;

  cuda_safe_call(cudaMemcpyToSymbol(cylOrient, &orient,   sizeof(bool),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylOrient");

  if (orient) {
    std::vector<cuFP_t> trans(9);
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++) trans[i*3+j] = cC->orient->transformBody()[i][j];
  
    cuda_safe_call(cudaMemcpyToSymbol(cylBody, &trans[0], sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying cylBody");
  }

  Component::cuRingType cr = *cC->cuRing.get();
  Cylinder::cuRingType  ar = *cuRing.get();

  // For debugging (set to false to disable)
  //
  static bool firstime = false;

  if (firstime) {
    testConstantsCyl<<<1, 1, 0, cr->stream>>>();
    cudaDeviceSynchronize();
    firstime = false;
  }
  
  // Zero counter and coefficients
  //
  unsigned Ntot = 0;
  use[0]      = 0.0;
  cylmass0[0] = 0.0;
  thrust::fill(host_coefs.begin(), host_coefs.end(), 0.0);

  // Zero out coefficient storage
  //
  zero_coefs();

  // Maximum radius on grid
  //
  cuFP_t rmax = rcylmax * acyl * M_SQRT1_2;

  // Copy particles to host vector
  //
  cC->ParticlesToCuda();

  // Loop over bunches
  //
  size_t psize  = cC->host_particles.size();

  Component::hostPartItr begin = cC->host_particles.begin();
  Component::hostPartItr first = begin;
  Component::hostPartItr last  = begin;
  Component::hostPartItr end   = cC->host_particles.end();

  if (psize <= cC->bunchSize) last = end;
  else std::advance(last, cC->bunchSize);

  // Set up stream and data arrays for asynchronous evaluation
  //
  std::vector<cudaStream_t> f_s;
  thrust::device_vector<unsigned int> f_use;
  thrust::device_vector<cuFP_t>       f_mass;

  while (std::distance(first, last)) {
    
    cr->first = first;
    cr->last  = last;
    cr->id    = ++dbg_id;

    // Copy bunch to device
    //
    cC->HostToDev(cr);

    // Sort particles and get coefficient size
    //
    PII lohi = cC->CudaSortByLevel(cr, mlevel, mlevel);

    // Compute grid
    //
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
    if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

#ifdef VERBOSE
    static debug_max_count = 10;
    static debug_cur_count = 0;
    if (debug_cur_count++ < debug_max_count) {
      std::cout << std::endl << "**" << std::endl
		<< "** N      = " << N           << std::endl
		<< "** I low  = " << lohi.first  << std::endl
		<< "** I high = " << lohi.second << std::endl
		<< "** Stride = " << stride      << std::endl
		<< "** Block  = " << BLOCK_SIZE  << std::endl
		<< "** Grid   = " << gridSize    << std::endl
		<< "** Xcen   = " << ctr[0]     << std::endl
		<< "** Ycen   = " << ctr[1]     << std::endl
		<< "** Zcen   = " << ctr[2]     << std::endl
		<< "**" << std::endl;
    }
#endif
  
    if (N) {

      // Adjust cached storage, if necessary
      //
      ar->resize_coefs(ncylorder, mmax, N, gridSize);

      // Shared memory size for the reduction
      //
      int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
    

      thrust::counting_iterator<int> index_begin(0);
      thrust::counting_iterator<int> index_end(gridSize*2*ncylorder);

      // Do the work
      //
				// Compute the coordinate
				// transformation
				// 
      coordKernelCyl<<<gridSize, BLOCK_SIZE, 0, cr->stream>>>
	(toKernel(cr->cuda_particles), toKernel(ar->m_d), toKernel(ar->p_d),
	 toKernel(ar->X_d), toKernel(ar->Y_d), toKernel(ar->iX_d),
	 toKernel(ar->iY_d), stride, lohi, rmax);
      
				// Compute the coefficient
				// contribution for each order
      int osize = ncylorder*2;	// 
      auto beg = ar->df_coef.begin();

      for (int m=0; m<=mmax; m++) {

	coefKernelCyl<<<gridSize, BLOCK_SIZE, 0, cr->stream>>>
	  (toKernel(ar->dN_coef), toKernel(t_d), toKernel(ar->m_d), toKernel(ar->p_d),
	   toKernel(ar->X_d), toKernel(ar->Y_d), toKernel(ar->iX_d), toKernel(ar->iY_d),
	   stride, m, ncylorder, lohi);
      
				// Begin the reduction per grid block
				//
	reduceSum<cuFP_t, BLOCK_SIZE><<<gridSize, BLOCK_SIZE, sMemSize, cr->stream>>>
	  (toKernel(ar->dc_coef), toKernel(ar->dN_coef), osize, N);
      
				// Finish the reduction for this order
				// in parallel
	thrust::reduce_by_key
	  (
	   thrust::cuda::par.on(cr->stream),
	   thrust::make_transform_iterator(index_begin, key_functor(gridSize)),
	   thrust::make_transform_iterator(index_end,   key_functor(gridSize)),
	   ar->dc_coef.begin(), thrust::make_discard_iterator(), ar->dw_coef.begin()
	   );

	thrust::transform(thrust::cuda::par.on(cr->stream),
			  ar->dw_coef.begin(), ar->dw_coef.end(),
			  beg, beg, thrust::plus<cuFP_t>());

	thrust::advance(beg, 2*ncylorder);
      }
    
      // Compute number and total mass of particles used in coefficient
      // determination
      //
      thrust::sort(thrust::cuda::par.on(cr->stream), ar->m_d.begin(), ar->m_d.end());

      // Asynchronously cache result for host side to prevent stream block
      //
      cudaStream_t s1;		// Create a new non blocking stream
      cudaStreamCreateWithFlags(&s1, cudaStreamNonBlocking);
      f_s.push_back(s1);

      size_t fsz = f_s.size();	// Augment the data vectors
      f_use. resize(fsz);
      f_mass.resize(fsz);
				// Call the kernel on a single thread
				// 
      reduceUseCyl<<<1, 1, 0, s1>>>(ar->m_d.begin(), ar->m_d.end(),
				    &f_use[fsz-1], &f_mass[fsz-1]);

      Ntot += N;
    }

    // Advance iterators
    //
    first = last;
    size_t nadv = std::distance(first, end);
    if (nadv < cC->bunchSize) last = end;
    else std::advance(last, cC->bunchSize);

    // Advance stream iterators
    //
    ++cr;			// Coefficient stream
    ++ar;			// Force method storage
  }

  if (Ntot == 0) {
    return;
  }

  // Accumulate the coefficients from the device to the host
  //
  for (auto r : cuRingData) {
    thrust::host_vector<cuFP_t> ret = r.df_coef;
    int offst = 0;
    for (int m=0; m<=mmax; m++) {
      for (size_t j=0; j<ncylorder; j++) {
	host_coefs[Imn(m, 'c', j, ncylorder)] += ret[2*j+offst];
	if (m>0) host_coefs[Imn(m, 's', j, ncylorder)] += ret[2*j+1+offst];
      }
      offst += 2*ncylorder;
    }
  }

  // Get the on-grid count and mass from the threads
  //
  for (auto & s : f_s) {	// Synchronize and dallocate streams
    cudaStreamSynchronize(s);
    cudaStreamDestroy(s);
  }
				// Copy the data from the device
  thrust::host_vector<unsigned int> f_ret1(f_use);
  thrust::host_vector<cuFP_t>       f_ret2(f_mass);

				// Sum counts and mass
  for (auto v : f_ret1) use[0]      += v;
  for (auto v : f_ret2) cylmass0[0] += v;

  // DEBUG, only useful for CUDAtest branch
  //
  if (false) {
    constexpr bool compareC = false;

    if (compareC) {
      std::cout << std::string(2*4+4*20, '-') << std::endl
		<< "---- Cylindrical "      << std::endl
		<< std::string(2*4+4*20, '-') << std::endl;
      std::cout << "M=0 coefficients" << std::endl
		<< std::setprecision(10);

      std::cout << std::setw(4)  << "n"
		<< std::setw(4)  << "i"
		<< std::setw(20) << "GPU"
		<< std::setw(20) << "CPU"
		<< std::setw(20) << "diff"
		<< std::setw(20) << "rel diff"
		<< std::endl;
    } else {
      std::cout << std::string(2*4+20, '-') << std::endl
		<< "---- Cylindrical "      << std::endl
		<< std::string(2*4+20, '-') << std::endl;
      std::cout << "M=0 coefficients" << std::endl
		<< std::setprecision(10);

      std::cout << std::setw(4)  << "n"
		<< std::setw(4)  << "i"
		<< std::setw(20) << "GPU"
		<< std::endl;
    }

    int i = Imn(0, 'c', 0, ncylorder);
    auto cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(0, 'c', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(0, n, 'c');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }

    std::cout << "M=1c coefficients" << std::endl;

    i = Imn(1, 'c', 0, ncylorder);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(1, 'c', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(1, n, 'c');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }

    std::cout << "M=1s coefficients" << std::endl;

    i = Imn(1, 's', 0, ncylorder);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(1, 's', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(1, n, 's');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }

    std::cout << "M=2c coefficients" << std::endl;

    i = Imn(2, 'c', 0, ncylorder);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(2, 'c', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(2, n, 'c');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }
    
    std::cout << "M=2s coefficients" << std::endl;

    i = Imn(2, 's', 0, ncylorder);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(2, 's', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(2, n, 's');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }

    std::cout << std::string(2*4+4*20, '-') << std::endl;
  }


  //
  // TEST comparison of coefficients for debugging
  //
  if (false) {

    struct Element
    {
      double d;
      double f;
      
      int  m;
      int  n;
      
      char cs;
    }
    elem;

    std::multimap<double, Element> compare;

    std::ofstream out("test_cyl.dat");

    // m loop
    for (int m=0; m<=mmax; m++) {
	
      if (m==0) {
	for (int n=0; n<ncylorder; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 'c';
	  elem.d = ortho->get_coef(m, n, 'c');
	  elem.f = host_coefs[Imn(m, 'c', n, ncylorder)];
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);
	  
	  compare.insert(std::make_pair(test, elem));;
	    
	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 'c'
	      << std::setw( 5) << Imn(m, 'c', n, ncylorder)
	      << std::setw(14) << elem.d
	      << std::setw(14) << elem.f
	      << std::endl;
	}

      } else {
	for (int n=0; n<ncylorder; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 'c';
	  elem.d = ortho->get_coef(m, n, 'c');
	  elem.f = host_coefs[Imn(m, 'c', n, ncylorder)];

	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 'c'
	      << std::setw( 5) << Imn(m, 'c', n, ncylorder)
	      << std::setw(14) << elem.d
	      << std::setw(14) << elem.f
	      << std::endl;
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);

	  compare.insert(std::make_pair(test, elem));;
	}

	for (int n=0; n<ncylorder; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 's';
	  elem.d = ortho->get_coef(m, n, 's');
	  elem.f = host_coefs[Imn(m, 's', n, ncylorder)];

	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 's'
	      << std::setw( 5) << Imn(m, 's', n-1, ncylorder)
	      << std::setw(14) << elem.d
	      << std::setw(14) << elem.f
	      << std::endl;
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);
	  
	  compare.insert(std::make_pair(test, elem));;
	}
      }
    }
    
    std::map<double, Element>::iterator best = compare.begin();
    std::map<double, Element>::iterator midl = best;
    std::advance(midl, compare.size()/2);
    std::map<double, Element>::reverse_iterator last = compare.rbegin();
    
    std::cout << std::string(3*2 + 3*20 + 20, '-') << std::endl
	      << "---- Cylinder coefficients" << std::endl
	      << std::string(3*2 + 3*20 + 20, '-') << std::endl;

    std::cout << "Best case: ["
	      << std::setw( 2) << best->second.m << ", "
	      << std::setw( 2) << best->second.n << ", "
	      << std::setw( 2) << best->second.cs << "] = "
	      << std::setw(20) << best->second.d
	      << std::setw(20) << best->second.f
	      << std::setw(20) << fabs(best->second.d - best->second.f)
	      << std::endl;
  
    std::cout << "Mid case:  ["
	      << std::setw( 2) << midl->second.m << ", "
	      << std::setw( 2) << midl->second.n << ", "
	      << std::setw( 2) << midl->second.cs << "] = "
	      << std::setw(20) << midl->second.d
	      << std::setw(20) << midl->second.f
	      << std::setw(20) << fabs(midl->second.d - midl->second.f)
	      << std::endl;
    
    std::cout << "Last case: ["
	      << std::setw( 2) << last->second.m << ", "
	      << std::setw( 2) << last->second.n << ", "
	      << std::setw( 2) << last->second.cs << "] = "
	      << std::setw(20) << last->second.d
	      << std::setw(20) << last->second.f
	      << std::setw(20) << fabs(last->second.d - last->second.f)
	      << std::endl;
  }

}


void Cylinder::determine_acceleration_cuda()
{
  if (initialize_cuda_cyl) {
    initialize_cuda();
    initialize_mapping_constants();
    initialize_cuda_cyl = false;
    // Only copy texture memory once
    t_d = tex;
  }

  std::cout << std::scientific;

  int deviceCount = 0;
  cuda_safe_call(cudaGetDeviceCount(&deviceCount),
		 __FILE__, __LINE__, "could not get device count");

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, deviceCount-1);

  Component::cuRingType cr = *cC->cuRing.get();

  // Assign component center and orientation
  //
  std::vector<cuFP_t> ctr;
  for (auto v : cC->getCenter(Component::Local | Component::Centered)) ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(cylCen, &ctr[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylCen");

  bool orient = (cC->EJ & Orient::AXIS) && !cC->EJdryrun;

  cuda_safe_call(cudaMemcpyToSymbol(cylOrient, &orient,   sizeof(bool),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylOrient");

  if (orient) {
    std::vector<cuFP_t> trans(9);
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++)
	trans[i*3+j] = cC->orient->transformBody()[i][j];
  
    cuda_safe_call(cudaMemcpyToSymbol(cylBody, &trans[0], sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying cylBody");

    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++)
	trans[i*3+j] = cC->orient->transformOrig()[i][j];
  
    cuda_safe_call(cudaMemcpyToSymbol(cylOrig, &trans[0], sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying cylOrig");
  }

  // Copy particles to host vector
  //
  cC->ParticlesToCuda();

  // Loop over bunches
  //
  size_t psize  = cC->host_particles.size();

  Component::hostPartItr begin = cC->host_particles.begin();
  Component::hostPartItr first = begin;
  Component::hostPartItr last  = begin;
  Component::hostPartItr end   = cC->host_particles.end();

  if (psize <= cC->bunchSize) last = end;
  else std::advance(last, cC->bunchSize);

  unsigned Ntot = 0;

  while (std::distance(first, last)) {

    cr->first = first;
    cr->last  = last;
    cr->id    = ++dbg_id;

    // Copy bunch to device
    //
    cC->HostToDev(cr);

    // Sort particles and get coefficient size
    //
    PII lohi = cC->CudaSortByLevel(cr, mlevel, multistep);

    // Compute grid
    //
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
    if (N>0) {

      Ntot += N;

      if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

#ifdef VERBOSE
      static debug_max_count = 10;
      static debug_cur_count = 0;
      if (debug_cur_count++ < debug_max_count) {
	std::cout << std::endl << "**" << std::endl
		  << "** N      = " << N          << std::endl
		  << "** Stride = " << stride     << std::endl
		  << "** Block  = " << BLOCK_SIZE << std::endl
		  << "** Grid   = " << gridSize   << std::endl
		  << "** Xcen   = " << ctr[0]     << std::endl
		  << "** Ycen   = " << ctr[1]     << std::endl
		  << "** Zcen   = " << ctr[2]     << std::endl
		  << "**" << std::endl;
      }
#endif
    
      // Shared memory size for the reduction
      //
      int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
      
      // Maximum radius on grid
      //
      cuFP_t rmax = rcylmax * acyl;
      
      // Do the work
      //
      forceKernelCyl<<<gridSize, BLOCK_SIZE, sMemSize, cr->stream>>>
	(toKernel(cr->cuda_particles), toKernel(dev_coefs), toKernel(t_d),
	 stride, mmax, ncylorder, lohi, rmax, cylmass, use_external);
      
      // Copy particles back to host
      //
      cC->DevToHost(cr);
    }
      
    // Advance iterators
    //
    first = last;
    size_t nadv = std::distance(first, end);
    if (nadv < cC->bunchSize) last = end;
    else std::advance(last, cC->bunchSize);

    // Advance stream iterators
    //
    ++cr;			// Component
				// Force method ring not needed here
  }

  cC->CudaToParticles();

  // DEBUGGING TEST
  if (false) {
    std::cout << std::string(10+7*16, '-') << std::endl;
    std::cout << "---- Acceleration in Cylinder [T=" << tnow
	      << ", N=" << Ntot << ", level=" << mlevel
	      << ", mass=" << cylmass
	      << ", name=" << cC->name << "]" << std::endl;
    std::cout << std::string(10+7*16, '-') << std::endl;
    first = last = begin;
    std::advance(last, 5);
    std::copy(first, last, std::ostream_iterator<cudaParticle>(std::cout, "\n"));
    first = begin;
    last  = end;
    std::advance(first, psize-5);
    std::copy(first, last, std::ostream_iterator<cudaParticle>(std::cout, "\n"));
    std::cout << std::string(10+7*16, '-') << std::endl;
  }
}

void Cylinder::HtoD_coefs()
{
  // Check size
  host_coefs.resize((2*mmax+1)*ncylorder); // Should stay fixed, no reserve

  // Copy from EmpCylSL
  
  // m loop
  //
  for (int m=0; m<=mmax; m++) {
    
    // n loop
    //
    for (int n=0; n<ncylorder; n++) {
      host_coefs[Imn(m, 'c', n, ncylorder)] = ortho->get_coef(m, n, 'c');
      if (m>0) host_coefs[Imn(m, 's', n, ncylorder)] = ortho->get_coef(m, n, 's');
    }
  }

  // Copy to device
  dev_coefs = host_coefs;
}


void Cylinder::DtoH_coefs(int M)
{
  // Copy from host device to EmpCylSL

  // m loop
  //
  for (int m=0; m<=mmax; m++) {
    
    // n loop
    //
    for (int n=0; n<ncylorder; n++) {
      ortho->set_coef(M, m, n, 'c') = host_coefs[Imn(m, 'c', n, ncylorder)];
      if (m>0) ortho->set_coef(M, m, n, 's') = host_coefs[Imn(m, 's', n, ncylorder)];
    }
  }
}

void chdata::destroy_cuda()
{
  for (size_t i=0; i<tex.size(); i++) {
    std::ostringstream sout;
    sout << "trying to free TextureObject [" << i << "]";
    cuda_safe_call(cudaDestroyTextureObject(tex[i]),
		   __FILE__, __LINE__, sout.str());
  }

  for (size_t i=0; i<cuInterpArray.size(); i++) {
    std::ostringstream sout;
    sout << "trying to free cuPitch [" << i << "]";
    cuda_safe_call(cudaFree(cuInterpArray[i]),
		     __FILE__, __LINE__, sout.str());
  }
}
