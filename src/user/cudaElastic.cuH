#ifndef _CUDAELASTIC_CUH
#define _CUDAELASTIC_CUH

#include <cudaUtil.cuH>

extern thrust::device_vector<cuFP_t> xsc_H, xsc_He, xsc_pH, xsc_pHe;
__constant__ extern cuFP_t cuH_Emin, cuHe_Emin, cuPH_Emin, cuPHe_Emin;

void cudaElasticInit();

__device__ cuFP_t cudaGeometric(int Z);

__device__ cuFP_t cudaElasticInterp(cuFP_t E, cuFP_t Emin, cuFP_t H, dArray<cuFP_t> xsc);

#endif
