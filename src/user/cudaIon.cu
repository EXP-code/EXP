// -*- C++ -*-

#include <iostream>
#include <iomanip>
#include <boost/make_shared.hpp>

#include <Ion.H>
#include <Timer.h>

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
  size_t ionSize = cuZ.size();

  // Interpolation data array
  //
  cuF0array.resize(ionSize, 0);
  cuFFarray.resize(ionSize, 0);
  cuRCarray.resize(ionSize, 0);
  cuCEarray.resize(ionSize, 0);
  cuCIarray.resize(ionSize, 0);

  ceEmin   .resize(ionSize, 0);
  ceEmax   .resize(ionSize, 0);
  ceDelE   .resize(ionSize, 0);

  ciEmin   .resize(ionSize, 0);
  ciEmax   .resize(ionSize, 0);
  ciDelE   .resize(ionSize, 0);

  NColl.    resize(ionSize, 0);
  NIonz.    resize(ionSize, 0);

  // Texture object array
  //
  ff_0.resize(ionSize);
  ff_d.resize(ionSize);
  rc_d.resize(ionSize);
  ce_d.resize(ionSize);
  ci_d.resize(ionSize);

  thrust::fill(ff_0.begin(), ff_0.end(), 0);
  thrust::fill(ff_d.begin(), ff_d.end(), 0);
  thrust::fill(rc_d.begin(), rc_d.end(), 0);
  thrust::fill(ce_d.begin(), ce_d.end(), 0);
  thrust::fill(ci_d.begin(), ci_d.end(), 0);

  size_t k = 0;

  for (auto v : IonList) {

    IonPtr I = v.second;
    cuZ[k] = I->Z;
    cuC[k] = I->C;
    
    // The free-free array
    if (cuC[k]>1) {
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

	/*
	for (int j = 0; j < tsize; j++) {	
	  std::cout << std::setw( 2) << I->Z
		    << std::setw( 2) << I->C
		    << std::setw( 4) << i
		    << std::setw( 4) << j
		    << std::setw(18) << I->freeFreeGrid[i][j]
		    << std::setw(18) << temp[j]
		    << std::endl;
	}
	*/


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

      /*
      for (int j=0; j<CHCUMK; j++) {
	std::cout << std::setw(12) << delC*j;
	for (int i = 0; i < std::min<int>(5, I->NfreeFreeGrid); i++) {
	  std::cout << std::setw(12) << h_buffer1[i+j*I->NfreeFreeGrid];
	}
	std::cout << std::endl;
      }
      */

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

      cuda_safe_call(cudaCreateTextureObject(&ff_0[k], &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");

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

      ceEmin[k] = I->collideEmin;
      ceEmax[k] = I->collideEmax;
      ceDelE[k] = I->delCollideE;
      NColl [k] = I->NcollideGrid;

      std::cout << " k=" << k
		<< " Emin=" << ceEmin[k]
		<< " Emax=" << ceEmax[k]
		<< " delE=" << ceDelE[k]
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
	(cudaCreateTextureObject(&ce_d[k], &resDesc, &texDesc, NULL),
	 __FILE__, __LINE__, "Failure in 2d texture creation");
      
      cuda_safe_call(cudaFree(d_Interp), __FILE__, __LINE__, "Failure freeing device memory");
    }

    if (cuC[k] <= cuZ[k]) {

      ciEmin[k] = I->ionizeEmin;
      ciEmax[k] = I->ionizeEmax;
      ciDelE[k] = I->DeltaEGrid;
      NIonz [k] = I->NionizeGrid;

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
      
      cuda_safe_call(cudaCreateTextureObject(&ci_d[k], &resDesc, &texDesc, NULL), __FILE__, __LINE__, "create texture object");
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


__global__ void testFreeFree
(dArray<cuFP_t> energy,
 dArray<cuFP_t> randsl,
 dArray<cuFP_t> ph, dArray<cuFP_t> xc,
 cudaTextureObject_t tex1,
 cudaTextureObject_t tex2)
{
  // value of h-bar * c in eV*nm
  //
  constexpr double hbc = 197.327;

  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {

    cuFP_t E = energy._v[tid];

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
    double rn = randsl._v[tid];
    double dC = 1.0/(CHCUMK-1);
    int lb    = rn/dC;
    cuFP_t k[4];

    // Interpolate the cross section array
    //
#if cuREAL == 4
    k[0]  = tex3D<float>(tex2, indx,   lb  , 0);
    k[1]  = tex3D<float>(tex2, indx+1, lb  , 0);
    k[2]  = tex3D<float>(tex2, indx,   lb+1, 0);
    k[3]  = tex3D<float>(tex2, indx+1, lb+1, 0);
#else
    k[0] = int2_as_double(tex3D<int2>(tex2, indx,   lb  , 0));
    k[1] = int2_as_double(tex3D<int2>(tex2, indx+1, lb  , 0));
    k[2] = int2_as_double(tex3D<int2>(tex2, indx,   lb+1, 0));
    k[3] = int2_as_double(tex3D<int2>(tex2, indx+1, lb+1, 0));
#endif

    // Linear interpolation
    //
    double a = (rn - dC*(lb+0)) / dC;
    double b = (dC*(lb+1) - rn) / dC;

    double K = A*(a*k[0] + b*k[2]) + B*(a*k[1] + b*k[3]);

    // Assign the photon energy
    //
    ph._v[tid] = pow(10, K) * hbc;

    // Use the integrated cross section from the differential grid
    //

    xc._v[tid] = 
#if cuREAL == 4
      A*tex1D<float>(tex1, indx  ) +
      B*tex1D<float>(tex1, indx+1) ;
#else
    A*int2_as_double(tex1D<int2>(tex1, indx  )) +
      B*int2_as_double(tex1D<int2>(tex1, indx+1)) ;
#endif
  }

  __syncthreads();
}

__global__ void testColExcite
(dArray<cuFP_t> energy,
 dArray<cuFP_t> ph, dArray<cuFP_t> xc,
 cudaTextureObject_t tex, cuFP_t Emin, cuFP_t Emax, cuFP_t delE, int Ngrid)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of evals
  //
  const unsigned int N = energy._s;

  if (tid < N) {

    cuFP_t E = energy._v[tid];

    if (E < Emin or E > Emax) {

      xc._v[tid] = 0.0;
      ph._v[tid] = 0.0;

    } else {

      // Interpolate the values
      //
      int indx = std::floor( (E - Emin)/delE );

      // Sanity check
      //
      if (indx > Ngrid-2) indx = Ngrid - 2;
      if (indx < 0)       indx = 0;
      
      double eA   = Emin + delE*indx;
      double eB   = Emin + delE*(indx+1);
    
      double A = (eB - E)/delE;
      double B = (E - eA)/delE;
    
#if cuREAL == 4
      xc._v[tid] = 
	A*tex3D<float>(tex, indx,   0, 0) +
	B*tex1D<float>(tex, indx+1, 0, 0) ;
      ph._v[tid] = 
	A*tex3D<float>(tex, indx,   1, 0) +
	B*tex1D<float>(tex, indx+1, 1, 0) ;
#else
      xc._v[tid] = 
	A*int2_as_double(tex3D<int2>(tex, indx  , 0, 0)) +
	B*int2_as_double(tex3D<int2>(tex, indx+1, 0, 0)) ;
      ph._v[tid] = 
	A*int2_as_double(tex3D<int2>(tex, indx  , 1, 0)) +
	B*int2_as_double(tex3D<int2>(tex, indx+1, 1, 0)) ;
#endif
    }

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

    unsigned int gridSize  = Nenergy/BLOCK_SIZE;
    if (Nenergy > gridSize*BLOCK_SIZE) gridSize++;

    cuda.start();

    if (cuC[k]>1)
      testFreeFree<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), toKernel(randsl_d),
					     toKernel(eFF_d), toKernel(xFF_d),
					     ff_0[k], ff_d[k]);

    if (cuC[k]<=cuZ[k])
      testColExcite<<<gridSize, BLOCK_SIZE>>>(toKernel(energy_d), 
					      toKernel(eCE_d), toKernel(xCE_d), ce_d[k],
					      ceEmin[k], ceEmax[k], ceDelE[k], NColl[k]);
      
    std::cout << "k=" << k << " delE=" << ceDelE[k] << std::endl;

    thrust::host_vector<cuFP_t> eFF_h = eFF_d;
    thrust::host_vector<cuFP_t> xFF_h = xFF_d;
    thrust::host_vector<cuFP_t> eCE_h = eCE_d;
    thrust::host_vector<cuFP_t> xCE_h = xCE_d;
    
    cuda.stop();
    
    std::vector<double> eFF_0(Nenergy, 0), xFF_0(Nenergy, 0);
    std::vector<double> eCE_0(Nenergy, 0), xCE_0(Nenergy, 0);
    
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
	std::cout << std::setw( 4) << cuZ[k]
		  << std::setw( 4) << cuC[k]
		  << std::setw(14) << energy_h[i]
		  << std::setw(14) << xCE_h[i]
		  << std::setw(14) << retCE.first
		  << std::endl;
      }
      if (retCE.second>0.0)
	eCE_0[i]   = (eCE_h[i] - retCE.second)/retCE.second;

      /*
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
    
    std::vector<double> quantiles = {0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99};

    std::cout << "Ion (" << I->Z << ", " << I->C << ")" << std::endl;
    for (auto v : quantiles) {
      int indx = std::min<int>(std::floor(v*Nenergy+0.5), Nenergy-1);
      double FF_xc = 0.0, FF_ph = 0.0, CE_xc = 0.0, CE_ph = 0.0;
      
      if (cuC[k]>1) {
	FF_xc = xFF_0[indx];
	FF_ph = eFF_0[indx];
      }

      if (cuC[k]<=cuZ[k]) {
	CE_xc = xCE_0[indx];
	CE_ph = eCE_0[indx];
      }

      std::cout << std::setw(10) << v
		<< " | " << std::setw(14) << FF_xc
		<< " | " << std::setw(14) << FF_ph
		<< " | " << std::setw(14) << CE_xc
		<< " | " << std::setw(14) << CE_ph
		<< std::endl;
    }

    k++;

  } // END: Ion list

  std::cout << std::endl
	    << "Serial time: " << serial() << std::endl
	    << "Cuda time  : " << cuda()   << std::endl;

}
