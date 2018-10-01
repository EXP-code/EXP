#ifndef _CUDAELASTIC_CUH
#define _CUDAELASTIC_CUH

#include <cudaUtil.cuH>

// Atomic radii in picometers from Clementi, E.; Raimond, D. L.;
// Reinhardt, W. P. (1967). "Atomic Screening Constants from SCF
// Functions. II. Atoms with 37 to 86 Electrons". Journal of Chemical
// Physics 47 (4): 1300-1307.  See also Paper 1, ref. therein.
//
const int numRadii = 87;

__constant__ extern int cudaRadii[numRadii];
extern thrust::device_vector<cuFP_t> xsc_H, xsc_He, xsc_pH, xsc_pHe;

void cudaElasticInit();

#endif
