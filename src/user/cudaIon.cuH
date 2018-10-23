// -*- C++ -*-

#ifndef _CUDAION_CUH
#define _CUDAION_CUH

#include <Ion.H>

__device__
void computeFreeFree
(cuFP_t E, cuFP_t rr, cuFP_t& ph, cuFP_t& xc, cuIonElement& elem);

__device__
void computeColExcite
(cuFP_t E, cuFP_t& ph, cuFP_t& xc, cuIonElement& elem);

__device__
void computeColIonize
(cuFP_t E, cuFP_t& xc, cuIonElement& elem);

__device__
void computeRadRecomb
(cuFP_t E, cuFP_t& xc, cuIonElement& elem);

__device__
void computePhotoIonize
(cuFP_t rr, cuFP_t& ph, cuFP_t& xc, cuIonElement& elem);

#endif


