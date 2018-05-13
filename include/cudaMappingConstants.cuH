#ifndef CUDA_MAPPING_CONSTANTS_H
#define CUDA_MAPPING_CONSTANTS_H

struct cudaMappingConstants
{
  float rscale;
  float hscale;
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  float dxi;
  float dyi;
  int   numr;
  int   numx;
  int   numy;
  int   cmap;
};

#endif
