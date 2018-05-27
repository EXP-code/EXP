#ifndef CUDA_MAPPING_CONSTANTS_H
#define CUDA_MAPPING_CONSTANTS_H

struct cudaMappingConstants
{
  double rscale;
  double hscale;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double dxi;
  double dyi;
  int    numr;
  int    numx;
  int    numy;
  int    cmap;
};

#endif
