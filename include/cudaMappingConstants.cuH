// -*- C++ -*-

#ifndef CUDA_MAPPING_CONSTANTS_H
#define CUDA_MAPPING_CONSTANTS_H

#include <cudaUtil.cuH>

struct cudaMappingConstants
{
  cuFP_t rscale;
  cuFP_t hscale;
  cuFP_t xmin;
  cuFP_t xmax;
  cuFP_t ymin;
  cuFP_t ymax;
  cuFP_t dxi;
  cuFP_t dyi;
  int    numr;
  int    numx;
  int    numy;
  int    cmapR;
  int    cmapZ;
};

#endif
