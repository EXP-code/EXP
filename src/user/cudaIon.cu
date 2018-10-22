// -*- C++ -*-

#include <iostream>
#include <iomanip>
#include <boost/make_shared.hpp>

#include <Ion.H>
#include <Timer.h>

thrust::host_vector<cuIonElement> cuIonElem;

// Global symbols for coordinate transformation
//
__device__ __constant__
cuFP_t ionEminGrid, ionEmaxGrid, ionDeltaEGrid;

__device__ __constant__
int ionEgridNumber, ionRadRecombNumber;

__global__
void testConstantsIon()
{
  printf("** Egrid(min) = %f\n", ionEminGrid);
  printf("** Egrid(max) = %f\n", ionEmaxGrid);
  printf("** Egrid(del) = %f\n", ionDeltaEGrid);
  printf("** Egrid(num) = %d\n", ionEgridNumber);
  printf("** Rgrid(num) = %d\n", ionRadRecombNumber);
}

void chdata::cuda_initialize_textures()
{
  size_t ionSize = IonList.size();

  // Interpolation data array
  //
  cuF0array.resize(ionSize, 0);
  cuFFarray.resize(ionSize, 0);
  cuRCarray.resize(ionSize, 0);
  cuCEarray.resize(ionSize, 0);
  cuCIarray.resize(ionSize, 0);
  cuPIarray.resize(ionSize, 0);

  // Texture object array
  //
  cuIonElem.resize(ionSize);

  // Total photo-ionization rate
  //
  std::vector<cuFP_t> phRate(ionSize, 0.0);

  size_t k = 0;

  for (auto v : IonList) {

    IonPtr I = v.second;
    cuIonElement& E = cuIonElem[k];

    // The free-free array
    if (E.C>1) {
      cudaTextureDesc texDesc;

      memset(&texDesc, 0, sizeof(texDesc));
      texDesc.readMode = cudaReadModeElementType;
      texDesc.filterMode = cudaFilterModePoint;
      texDesc.addressMode[0] = cudaAddressModeClamp;
      texDesc.addressMode[1] = cudaAddressModeClamp;
      texDesc.addressMode[2] = cudaAddressModeClamp;
  
      // Temporary storage
      //
      std::vector<cuFP_t> h_buffer0(I->NfreeFreeGrid, 0.0);

      cuFP_t *d_Interp;

      cuda_safe_call(cudaMalloc((void **)&d_Interp, I->NfreeFreeGrid*CHCUMK*sizeof(cuFP_t)),
		     __FILE__, __LINE__,
		     "Error allocating d_Interp1 for texture construction");
  
      std::vector<cuFP_t> h_buffer1(I->NfreeFreeGrid*CHCUMK, 0.0);

      double delC = 1.0/(CHCUMK-1);

      // Copy cross section values to buffer
      //
      for (int i = 0; i < I->NfreeFreeGrid; i++) {

	h_buffer0[i] = I->freeFreeGrid[i].back();
	
	// Unit normalized cumulative distribution
	//
	size_t tsize = I->freeFreeGrid[i].size();
	std::vector<double> temp(tsize);
	for (int j = 0; j < tsize; j++) {	
	  temp[j] = I->freeFreeGrid[i][j]/h_buffer0[i];
	}

	// End points
	//
	h_buffer1[i                              ] = I->kgrid[0];
	h_buffer1[i + (CHCUMK-1)*I->NfreeFreeGrid] = I->kgrid[tsize-1];

	// Remap to even grid
	//
	for (int j=1; j<CHCUMK-1; j++) {

	  double C = delC*j;

	  // Points to first element that is not < C
	  // but may be equal
	  std::vector<double>::iterator lb = 
	    std::lower_bound(temp.begin(), temp.end(), C);
    
	  // Assign upper end of range to the
	  // found element
	  //
	  std::vector<double>::iterator ub = lb;
	  //
	  // If is the first element, increment
	  // the upper boundary
	  //
	  if (lb == temp.begin()) { if (temp.size()>1) ub++; }
	  //
	  // Otherwise, decrement the lower boundary
	  //
	  else { lb--; }
    
	  // Compute the associated indices
	  //
	  size_t ii = lb - temp.begin();
	  size_t jj = ub - temp.begin();
	  double kk = I->kgrid[ii];
	  
	  // Linear interpolation
	  //
	  if (*ub > *lb) {
	    double d = *ub - *lb;
	    double a = (C - *lb) / d;
	    double b = (*ub - C) / d;
	    /*
	    std::cout << "[a, b] = [" << a << ", " << b << "]"
		      << ", c = " << C << std::endl;
	    */
	    kk  = a * I->kgrid[ii] + b * I->kgrid[jj];
	  }

	  h_buffer1[i + j*I->NfreeFreeGrid] = kk;

	} // END: cumululative array loop

      } // END: energy loop

      // Copy 1-dim data to device
      //
      size_t tsize = I->NfreeFreeGrid*sizeof(cuFP_t);

      // Allocate CUDA array in device memory (a one-dimension 'channel')
      //
#if cuREAL == 4
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<float>();
#else
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<int2>();
#endif
      
      std::cout << "Allocating cuF0array[" << k << "]" << std::endl;
      cuda_safe_call(cudaMallocArray(&cuF0array[k], &channelDesc1, I->NfreeFreeGrid), __FILE__, __LINE__, "malloc cuArray");

      cuda_safe_call(cudaMemcpyToArray(cuF0array[k], 0, 0, &h_buffer0[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "copy texture to array");

      // Specify 1-d texture

      cudaResourceDesc resDesc;

      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array = cuF0array[k];

      cuda_safe_call(cudaCreateTextureObject(&E.ff_0, &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");

      // Copy data to device
      tsize = I->NfreeFreeGrid*CHCUMK*sizeof(cuFP_t);
      cuda_safe_call(cudaMemcpy(d_Interp, &h_buffer1[0], tsize, cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying texture table to device");
    
      // cuda 2d Array Descriptor
      //
#if cuREAL == 4
      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();
#else
      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<int2>();
#endif
      // cuda 2d Array
      //
      cuda_safe_call(cudaMalloc3DArray(&cuFFarray[k], &channelDesc2, make_cudaExtent(I->NfreeFreeGrid, CHCUMK, 1), 0), __FILE__, __LINE__, "Error allocating cuArray for 3d texture");
      
      // Array creation
      //
      cudaMemcpy3DParms copyParams = {0};
  
      copyParams.srcPtr   = make_cudaPitchedPtr(d_Interp, I->NfreeFreeGrid*sizeof(cuFP_t), I->NfreeFreeGrid, CHCUMK);
      copyParams.dstArray = cuFFarray[k];
      copyParams.extent   = make_cudaExtent(I->NfreeFreeGrid, CHCUMK, 1);
      copyParams.kind     = cudaMemcpyDeviceToDevice;
      
      cuda_safe_call(cudaMemcpy3D(&copyParams), __FILE__, __LINE__, "Error in copying 3d pitched array");

      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array  = cuFFarray[k];
    
      cuda_safe_call
	(cudaCreateTextureObject(&E.ff_d, &resDesc, &texDesc, NULL),
	 __FILE__, __LINE__, "Failure in 2d texture creation");
      
      cuda_safe_call(cudaFree(d_Interp), __FILE__, __LINE__, "Failure freeing device memory");
    }

    // Radiative recombination texture (1-d)
    //
    if (E.C>1) {
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
      
      cuda_safe_call(cudaCreateTextureObject(&E.rc_d, &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");
    }

    // The collisional excitation array

    if (E.C <= E.Z) {

      E.ceEmin = I->collideEmin;
      E.ceEmax = I->collideEmax;
      E.ceDelE = I->delCollideE;
      E.NColl  = I->NcollideGrid;

      std::cout << " k=" << k
		<< " Emin=" << E.ceEmin
		<< " Emax=" << E.ceEmax
		<< " delE=" << E.ceDelE
		<< std::endl;

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
      std::cout << "Size(" << I->Z << ", " << I->C << ")="
		<< I->NcollideGrid << std::endl;
      cuda_safe_call(cudaMalloc((void **)&d_Interp, I->NcollideGrid*2*sizeof(cuFP_t)),
		     __FILE__, __LINE__,
		     "Error allocating d_Interp for texture construction");
  
      std::vector<cuFP_t> h_buffer(I->NcollideGrid*2, 0.0);

      // Copy vectors to buffer
      //
      for (int i = 0; i < I->NcollideGrid; i++) {
	h_buffer[i                  ] = I->collideDataGrid[i].back().first;
	h_buffer[i + I->NcollideGrid] = I->collideDataGrid[i].back().second;
      }
      
      // Copy data to device
      cuda_safe_call(cudaMemcpy(d_Interp, &h_buffer[0], I->NcollideGrid*2*sizeof(cuFP_t), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying texture table to device");
    
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
      
      copyParams.srcPtr   = make_cudaPitchedPtr(d_Interp, I->NcollideGrid*sizeof(cuFP_t), I->NcollideGrid, 2);
      copyParams.dstArray = cuCEarray[k];
      copyParams.extent   = make_cudaExtent(I->NcollideGrid, 2, 1);
      copyParams.kind     = cudaMemcpyDeviceToDevice;
      
      cuda_safe_call(cudaMemcpy3D(&copyParams), __FILE__, __LINE__, "Error in copying 3d pitched array");
      
      cudaResourceDesc resDesc;
      
      memset(&resDesc, 0, sizeof(cudaResourceDesc));
      resDesc.resType = cudaResourceTypeArray;
      resDesc.res.array.array  = cuCEarray[k];
      
      cuda_safe_call
	(cudaCreateTextureObject(&E.ce_d, &resDesc, &texDesc, NULL),
	 __FILE__, __LINE__, "Failure in 2d texture creation");
      
      cuda_safe_call(cudaFree(d_Interp), __FILE__, __LINE__, "Failure freeing device memory");
    }

    if (E.C <= E.Z) {

      E.ciEmin = I->ionizeEmin;
      E.ciEmax = I->ionizeEmax;
      E.ciDelE = I->DeltaEGrid;
      E.NIonz  = I->NionizeGrid;

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
      
      thrust::host_vector<cuFP_t> tt(I->NionizeGrid);
      
      std::cout << "Size(" << I->Z << ", " << I->C << ")="
		<< I->NionizeGrid << std::endl;

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
      
      cuda_safe_call(cudaCreateTextureObject(&E.ci_d, &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");

      // Photoionization array
      //

      thrust::host_vector<cuFP_t> piCum.resize(CHCUMK, 0.0);
      piCum[CHCUMK-1] = 1.0;
      
      double delC = 1.0/(CHCUMK-1);
      
      if (not IBinit) IBcreate();
      
      E.piTotl = IBtotl;

      // Copy cross section values to buffer
      //
      for (int i = 1; i < CHCUMK-1; i++) {

	// Location in cumulative cross section grid
	//
	double C = delC*j;

	// Interpolate the cross section array
	//
	
	// Points to first element that is not < rn
	// but may be equal
	std::vector<double>::iterator lb = 
	  std::lower_bound(IBcum.begin(), IBcum.end(), rn);
	
	// Assign upper end of range to the
	// found element
	//
	std::vector<double>::iterator ub = lb;
	//
	// If is the first element, increment
	// the upper boundary
	//
	if (lb == IBcum.begin()) { if (IBcum.size()>1) ub++; }
	//
	// Otherwise, decrement the lower boundary
	//
	else { lb--; }
	
	// Compute the associated indices
	//
	size_t ii = lb - IBcum.begin();
	size_t jj = ub - IBcum.begin();
	double nu = nugrid[ii];
	  
	// Linear interpolation
	//
	if (*ub > *lb) {
	  double d = *ub - *lb;
	  double a = (rn - *lb) / d;
	  double b = (*ub - rn) / d;
	  nu  = a * nugrid[ii] + b * nugrid[jj];
	}
    
	piCum[i] = (nu - 1.0)*ip;
	
      } // END: cumululative array loop

      std::cout << "Allocating pi_0[" << k << "]" << std::endl;

      // Create storage on device
      cuPIarray[k] = piCum;

      // Assign pointer
      E.pi_0 = thrust::raw_pointer_cast(&cuPIarray[k][0]);
    }
    
    // Increment counter
    k++;	
    
  } // END: IonList

}

void chdata::cuda_initialize_grid_constants()
{
  double Emin, Emax, delE;
  int NE, NR;

  for (auto v : IonList) {
    Emin = v.second->EminGrid;
    Emax = v.second->EmaxGrid;
    delE = v.second->DeltaEGrid;

    NE   = v.second->NfreeFreeGrid;

    if (v.first.second>1) {
      NR = v.second->NradRecombGrid;
      break;
    }
  }

  cuFP_t f;

  // Copy constants to device
  //
  cuda_safe_call(cudaMemcpyToSymbol(ionEminGrid, &(f=Emin),
				    sizeof(cuFP_t), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionEminGrid");

  cuda_safe_call(cudaMemcpyToSymbol(ionEmaxGrid, &(f=Emax),
				    sizeof(cuFP_t), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionEmaxGrid");

  cuda_safe_call(cudaMemcpyToSymbol(ionDeltaEGrid, &(f=delE),
				    sizeof(cuFP_t), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionDeltaEGrid");

  cuda_safe_call(cudaMemcpyToSymbol(ionEgridNumber, &NE,
				    sizeof(int), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionEgridNumber");

  cuda_safe_call(cudaMemcpyToSymbol(ionRadRecombNumber, &NR,
				    sizeof(int), size_t(0),
				    cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying ionRadRecombNumber");
}


__device__
void computeFreeFree
(cuFP_t E, cuFP_t rr, cuFP_t& ph, cuFP_t& xc, cuIonElement& elem)
{
  // value of h-bar * c in eV*nm
  //
  constexpr double hbc = 197.327;

  // Enforce minimum and maximum energies
  //
  if (E<ionEminGrid) E = ionEminGrid;
  if (E>ionEmaxGrid) E = ionEmaxGrid;

  size_t indx = std::floor( (E - ionEminGrid)/ionDeltaEGrid );
    
  if (indx >= ionEgridNumber - 1) indx = ionEgridNumber-2;

  double eA = ionEminGrid + ionDeltaEGrid*indx;
  double eB = ionEminGrid + ionDeltaEGrid*(indx+1);
  
  double A = (eB - E)/ionDeltaEGrid;
  double B = (E - eA)/ionDeltaEGrid;
  
  // Location in cumulative cross section grid
  //
  double rn = rr;
  double dC = 1.0/(CHCUMK-1);
  int lb    = rn/dC;
  cuFP_t k[4];

  // Interpolate the cross section array
  //
#if cuREAL == 4
  k[0]  = tex3D<float>(elem.ff_d, indx,   lb  , 0);
  k[1]  = tex3D<float>(elem.ff_d, indx+1, lb  , 0);
  k[2]  = tex3D<float>(elem.ff_d, indx,   lb+1, 0);
  k[3]  = tex3D<float>(elem.ff_d, indx+1, lb+1, 0);
#else
  k[0] = int2_as_double(tex3D<int2>(elem.ff_d, indx,   lb  , 0));
  k[1] = int2_as_double(tex3D<int2>(elem.ff_d, indx+1, lb  , 0));
  k[2] = int2_as_double(tex3D<int2>(elem.ff_d, indx,   lb+1, 0));
  k[3] = int2_as_double(tex3D<int2>(elem.ff_d, indx+1, lb+1, 0));
#endif
  
  // Linear interpolation
  //
  double a = (rn - dC*(lb+0)) / dC;
  double b = (dC*(lb+1) - rn) / dC;

  double K = A*(a*k[0] + b*k[2]) + B*(a*k[1] + b*k[3]);

  // Assign the photon energy
  //
  ph = pow(10, K) * hbc;

  // Use the integrated cross section from the differential grid
  //

  xc = 
#if cuREAL == 4
    A*tex1D<float>(elem.ff_0, indx  ) +
    B*tex1D<float>(elem.ff_0, indx+1) ;
#else
    A*int2_as_double(tex1D<int2>(elem.ff_0, indx  )) +
    B*int2_as_double(tex1D<int2>(elem.ff_0, indx+1)) ;
#endif
}


__global__
void testFreeFree
(dArray<cuFP_t> energy,
 dArray<cuFP_t> randsl,
 dArray<cuFP_t> ph, dArray<cuFP_t> xc,
 cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    computeFreeFree(energy._v[tid], randsl._v[tid], 
		    ph._v[tid], xc._v[tid], elem);
  }

  __syncthreads();
}


__device__
void computeColExcite
(cuFP_t E, cuFP_t& ph, cuFP_t& xc, cuIonElement& elem)
{
  if (E < elem.ceEmin or E > elem.ceEmax) {

    xc = 0.0;
    ph= 0.0;

  } else {

    // Interpolate the values
    //
    int indx = std::floor( (E - elem.ceEmin)/elem.ceDelE );
    
    // Sanity check
    //
    if (indx > elem.NColl-2) indx = elem.NColl - 2;
    if (indx < 0)            indx = 0;
    
    double eA   = elem.ceEmin + elem.ceDelE*indx;
    double eB   = elem.ceEmin + elem.ceDelE*(indx+1);
    
    double A = (eB - E)/elem.ceDelE;
    double B = (E - eA)/elem.ceDelE;
    
#if cuREAL == 4
    xc = 
      A*tex3D<float>(elem.ce_d, indx,   0, 0) +
      B*tex3D<float>(elem.ce_d, indx+1, 0, 0) ;
    ph = 
      A*tex3D<float>(elem.ce_d, indx,   1, 0) +
      B*tex3D<float>(elem.ce_d, indx+1, 1, 0) ;
#else
    xc = 
      A*int2_as_double(tex3D<int2>(elem.ce_d, indx  , 0, 0)) +
      B*int2_as_double(tex3D<int2>(elem.ce_d, indx+1, 0, 0)) ;
    ph= 
      A*int2_as_double(tex3D<int2>(elem.ce_d, indx  , 1, 0)) +
      B*int2_as_double(tex3D<int2>(elem.ce_d, indx+1, 1, 0)) ;
#endif
  }
  // DONE
}

__global__ void testColExcite
(dArray<cuFP_t> energy,
 dArray<cuFP_t> ph, dArray<cuFP_t> xc, cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    computeColExcite(energy._v[tid], ph._v[tid], xc._v[tid], elem);
  }

  __syncthreads();
}

__device__
void computeColIonize
(cuFP_t E, cuFP_t& xc, cuIonElement& elem)
{
  if (E < elem.ciEmin or E > elem.ciEmax) {

    xc = 0.0;

  } else {

    // Interpolate the values
    //
    int indx = std::floor( (E - elem.ciEmin)/elem.ciDelE );

    // Sanity check
    //
    if (indx > elem.NIonz-2) indx = elem.NIonz - 2;
    if (indx < 0)            indx = 0;
    
    double eA   = elem.ciEmin + elem.ciDelE*indx;
    double eB   = elem.ciEmin + elem.ciDelE*(indx+1);
    
    double A = (eB - E)/elem.ciDelE;
    double B = (E - eA)/elem.ciDelE;
    
#if cuREAL == 4
    xc = 
      A*tex1D<float>(elem.ci_d, indx  ) +
      B*tex1D<float>(elem.ci_d, indx+1) ;
#else
    xc = 
      A*int2_as_double(tex1D<int2>(elem.ci_d, indx  )) +
      B*int2_as_double(tex1D<int2>(elem.ci_d, indx+1)) ;
#endif
  }
}


__device__
void computePhotoIonize
(cuFP_t rr, cuFP_t& ph, cuFP_t& xc, cuIonElement& elem)
{
  cuFP_t dC = 1.0/CHCUMK;
  int indx  = rr/ionNuCdel;
  if (indx > CHCUMK-2) indx = CHCUMK - 2;

  // Linear interpolation
  //
  double a = (rr - dC*(indx+0)) / dC;
  double b = (dC*(indx+1) - rr) / dC;

  ph = a*elem.pi_0[indx+0] + b*elem.pi_0[indx+1];
  xc = elem.piTotl;
}


__global__ void testColIonize
(dArray<cuFP_t> energy, dArray<cuFP_t> xc, cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    computeColIonize(energy._v[tid], xc._v[tid], elem);
  }

  __syncthreads();
}

__device__
void computeRadRecomb
(cuFP_t E, cuFP_t& xc, cuIonElement& elem)
{
  if (E < ionEminGrid or E > ionEmaxGrid) {

    xc = 0.0;

  } else {

    // Interpolate the values
    //
    int indx = std::floor( (E - ionEminGrid)/ionDeltaEGrid );

    // Sanity check
    //
    if (indx > ionRadRecombNumber-2) indx = ionRadRecombNumber - 2;
    if (indx < 0)                    indx = 0;
    
    double eA   = ionEminGrid + ionDeltaEGrid*indx;
    double eB   = ionEminGrid + ionDeltaEGrid*(indx+1);
    
    double A = (eB - E)/ionDeltaEGrid;
    double B = (E - eA)/ionDeltaEGrid;
    
#if cuREAL == 4
    xc = 
      A*tex1D<float>(elem.rc_d, indx  ) +
      B*tex1D<float>(elem.rc_d, indx+1) ;
#else
    xc = 
      A*int2_as_double(tex1D<int2>(elem.rc_d, indx  )) +
      B*int2_as_double(tex1D<int2>(elem.rc_d, indx+1)) ;
#endif
  }
  // DONE
}

__global__
void testRadRecomb
(dArray<cuFP_t> energy, dArray<cuFP_t> xc, cuIonElement elem)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {
    computeRadRecomb(energy._v[tid], xc._v[tid], elem);
  }

  __syncthreads();
}


void chdata::testCross(int Nenergy)
{
  // Timers
  //
  Timer serial, cuda;

  // Loop over ions and tabulate statistics
  //
  size_t k = 0;

  thrust::host_vector<cuFP_t> energy_h(Nenergy), randsl_h(Nenergy);

  for (auto v : IonList) {

    IonPtr I = v.second;
    cuIonElement& E = cuIonElem[k];

    // Make an energy grid
    //
    double dE = (I->EmaxGrid - I->EminGrid)/(Nenergy-1) * 0.999;
    for (int i=0; i<Nenergy; i++) {
      energy_h[i] = I->EminGrid + dE*i;
      randsl_h[i] = static_cast<cuFP_t>(rand())/RAND_MAX;
    }

    thrust::device_vector<cuFP_t> energy_d = energy_h;
    thrust::device_vector<cuFP_t> randsl_d = randsl_h;

    // Only free-free for non-neutral species

    thrust::device_vector<cuFP_t> eFF_d(Nenergy), xFF_d(Nenergy);
    thrust::device_vector<cuFP_t> eCE_d(Nenergy), xCE_d(Nenergy);
    thrust::device_vector<cuFP_t> xCI_d(Nenergy), xRC_d(Nenergy);

    unsigned int gridSize  = Nenergy/BLOCK_SIZE;
    if (Nenergy > gridSize*BLOCK_SIZE) gridSize++;

    cuda.start();

    if (E.C>1)
      testFreeFree<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), toKernel(randsl_d),
					     toKernel(eFF_d), toKernel(xFF_d),
					     cuIonElem[k]);

    if (E.C<=E.Z)
      testColExcite<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(eCE_d), toKernel(xCE_d), cuIonElem[k]);
      
    if (E.C<=E.Z)
      testColIonize<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(xCI_d), cuIonElem[k]);
      
    if (E.C>1)
      testRadRecomb<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(xRC_d), cuIonElem[k]);
      
    std::cout << "k=" << k << " delE=" << E.ceDelE << std::endl;

    thrust::host_vector<cuFP_t> eFF_h = eFF_d;
    thrust::host_vector<cuFP_t> xFF_h = xFF_d;
    thrust::host_vector<cuFP_t> eCE_h = eCE_d;
    thrust::host_vector<cuFP_t> xCE_h = xCE_d;
    thrust::host_vector<cuFP_t> xCI_h = xCI_d;
    thrust::host_vector<cuFP_t> xRC_h = xRC_d;
    
    cuda.stop();
    
    std::vector<double> eFF_0(Nenergy, 0), xFF_0(Nenergy, 0);
    std::vector<double> eCE_0(Nenergy, 0), xCE_0(Nenergy, 0);
    std::vector<double> xCI_0(Nenergy, 0), xRC_0(Nenergy, 0);
    
    serial.start();
    
    for (int i=0; i<Nenergy; i++) {
				// Free-free
      auto retFF = I->freeFreeCrossTest(energy_h[i], randsl_h[i], 0);
      if (retFF.first>0.0)
	xFF_0[i]   = (xFF_h[i] - retFF.first )/retFF.first;
      if (retFF.second>0.0)
	eFF_0[i]   = (eFF_h[i] - retFF.second)/retFF.second;

				// Collisional excitation
      auto retCE = I->collExciteCross(energy_h[i], 0).back();
      if (retCE.first>0.0) {
	xCE_0[i]   = (xCE_h[i] - retCE.first )/retCE.first;
	/*
	std::cout << std::setw( 4) << cuZ[k]
		  << std::setw( 4) << cuC[k]
		  << std::setw(14) << energy_h[i]
		  << std::setw(14) << xCE_h[i]
		  << std::setw(14) << retCE.first
		  << std::endl;
	*/
      }
				// Collisional ionization

      auto retCI = I->directIonCross(energy_h[i], 0);
      if (retCI>0.0) {
	xCI_0[i]   = (xCI_h[i] - retCI)/retCI;
	/*
	std::cout << std::setw( 4) << cuZ[k]
		  << std::setw( 4) << cuC[k]
		  << std::setw(14) << energy_h[i]
		  << std::setw(14) << xCI_h[i]
		  << std::setw(14) << retCI
		  << std::endl;
	*/
      }

				// Radiative recombination

      auto retRC = I->radRecombCross(energy_h[i], 0).back();
      if (retRC>0.0) {
	xRC_0[i]   = (xRC_h[i] - retRC)/retRC;
	/*
	std::cout << std::setw( 4) << cuZ[k]
		  << std::setw( 4) << cuC[k]
		  << std::setw(14) << energy_h[i]
		  << std::setw(14) << xRC_h[i]
		  << std::setw(14) << retRC
		  << std::endl;
	*/
      }

      /*
      if (retCE.second>0.0)
	eCE_0[i]   = (eCE_h[i] - retCE.second)/retCE.second;

      if (cuC[k]<=cuZ[k])
	std::cout << std::setw(14) << xCE_h[i]
		  << std::setw(14) << eCE_h[i]
		  << std::endl;
      */
    }

    serial.stop();

    std::sort(xFF_0.begin(), xFF_0.end());
    std::sort(eFF_0.begin(), eFF_0.end());
    std::sort(xCE_0.begin(), xCE_0.end());
    std::sort(eCE_0.begin(), eCE_0.end());
    std::sort(xCI_0.begin(), xCI_0.end());
    std::sort(xRC_0.begin(), xRC_0.end());
    
    std::vector<double> quantiles = {0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99};

    std::cout << "Ion (" << I->Z << ", " << I->C << ")" << std::endl;
    for (auto v : quantiles) {
      int indx = std::min<int>(std::floor(v*Nenergy+0.5), Nenergy-1);
      double FF_xc = 0.0, FF_ph = 0.0, CE_xc = 0.0, CE_ph = 0.0;
      double CI_xc = 0.0, RC_xc = 0.0;
      
      if (E.C>1) {
	FF_xc = xFF_0[indx];
	FF_ph = eFF_0[indx];
	RC_xc = xRC_0[indx];
      }

      if (E.C<=E.Z) {
	CE_xc = xCE_0[indx];
	CE_ph = eCE_0[indx];
	CI_xc = xCI_0[indx];
      }

      std::cout << std::setw(10) << v
		<< " | " << std::setw(14) << FF_xc
		<< " | " << std::setw(14) << FF_ph
		<< " | " << std::setw(14) << CE_xc
		<< " | " << std::setw(14) << CE_ph
		<< " | " << std::setw(14) << CI_xc
		<< " | " << std::setw(14) << RC_xc
		<< std::endl;
    }

    k++;

  } // END: Ion list

  std::cout << std::endl
	    << "Serial time: " << serial() << std::endl
	    << "Cuda time  : " << cuda()   << std::endl;
}
