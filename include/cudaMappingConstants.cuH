#ifndef CUDA_MAPPING_CONSTANTS_H
#define CUDA_MAPPING_CONSTANTS_H

struct cudaMappingConstants
{
  float scale;
  float xmin;
  float xmax;
  float dxi;
  int   numr;
  int   cmap;
};

#endif
