#include <Component.H>
#include <Cylinder.H>
#include <cudaReduce.cuH>

// Define for debugging
//
// #define OFF_GRID_ALERT
// #define BOUNDS_CHECK
// #define VERBOSE

// Global symbols for coordinate transformation
//
__device__ __constant__
cuFP_t cylRscale, cylHscale, cylXmin, cylXmax, cylYmin, cylYmax, cylDxi, cylDyi, cylCen[3], cylBody[9], cylOrig[9];

__device__ __constant__
int   cylNumx, cylNumy, cylCmap;

__device__ __constant__
bool  cylOrient;

__host__ __device__
int Imn(int m, char cs, int n, int nmax)
{
  int ret = 0;

  if (m==0) ret = n;
  else ret = (2*m - 1 + (cs=='s' ? 1 : 0))*nmax + n;

#ifdef BOUNDS_CHECK
  // Verbose sanity check
  if (ret >= (2*m+1)*nmax) {
    printf("Imn oab: %4d %4d %4d [%4d : %4d ]\n", m, n, ret, (2*m+1)*nmax, nmax);
  }
#endif
  return ret;
}

__global__
void testConstantsCyl()
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

				// R coordinate transformation
__device__
cuFP_t cu_r_to_xi_cyl(cuFP_t r)
{
  cuFP_t ret;

  if (cylCmap==1) {
    ret = (r/cylRscale - 1.0)/(r/cylRscale + 1.0);
  } else {
    ret = r;
  }    

  return ret;
}
    
__device__
cuFP_t cu_xi_to_r_cyl(cuFP_t xi)
{
  cuFP_t ret;

  if (cylCmap==1) {
    ret = (1.0 + xi)/(1.0 - xi) * cylRscale;
  } else {
    ret = xi;
  }

  return ret;
}

__device__
cuFP_t cu_d_xi_to_r_cyl(cuFP_t xi)
{
  cuFP_t ret;

  if (cylCmap==1) {
    ret = 0.5*(1.0 - xi)*(1.0 - xi) / cylRscale;
  } else {
    ret = 1.0;
  }

  return ret;
}

				// Z coordinate transformation
__device__
cuFP_t cu_z_to_y_cyl(cuFP_t z)
{ return z/(fabs(z)+FLT_MIN)*asinh(fabs(z/cylHscale)); }

__device__
cuFP_t cu_y_to_z_cyl(cuFP_t y)
{ return cylHscale*sinh(y); }

__device__
cuFP_t cu_d_y_to_z_cyl(cuFP_t y)
{ return cylHscale*cosh(y); }


void Cylinder::initialize_mapping_constants()
{
  // Copy constants to device
  //
  
  cudaMappingConstants f = getCudaMappingConstants();

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
	cfrac      = 1.0 - mfactor;
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

static bool initialize_cuda_cyl = true;

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

  // Sort particles and get coefficient size
  //
  PII lohi = cC->CudaSortByLevel(mlevel, mlevel);

  // Zero out coefficients
  //
  host_coefs.resize((2*mmax+1)*ncylorder);

  // Compute grid
  //
  unsigned int N         = lohi.second - lohi.first;
  unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  unsigned int gridSize  = N/BLOCK_SIZE/stride;

  if (N == 0) {
    use[0] = 0.0;
    cylmass0[0] = 0.0;
    return;
  }

  if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

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

#ifdef VERBOSE
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
#endif

  // Create space for coefficient reduction
  //
  dN_coef.resize(2*ncylorder*N);
  dc_coef.resize(2*ncylorder*gridSize);
  df_coef.resize(2*ncylorder);

  // Space for coordinate arrays
  //
  m_d.resize(N);
  X_d.resize(N);
  Y_d.resize(N);
  p_d.resize(N);
  iX_d.resize(N);
  iY_d.resize(N);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);

  // For debugging (set to false to disable)
  //
  static bool firstime = false;

  if (firstime) {
    testConstantsCyl<<<1, 1>>>();
    cudaDeviceSynchronize();
    firstime = false;
  }

  thrust::counting_iterator<int> index_begin(0);
  thrust::counting_iterator<int> index_end(gridSize*2*ncylorder);

  // Maximum radius on grid
  //
  cuFP_t rmax = rcylmax * acyl * M_SQRT1_2;

  // Do the work
  //
				// Compute the coordinate
				// transformation
				// 
  coordKernelCyl<<<gridSize, BLOCK_SIZE>>>
    (toKernel(cC->cuda_particles), toKernel(m_d), toKernel(p_d),
     toKernel(X_d), toKernel(Y_d), toKernel(iX_d), toKernel(iY_d),
     stride, lohi, rmax);

				// Compute the coefficient
				// contribution for each order
  int osize = ncylorder*2;	// 
  for (int m=0; m<=mmax; m++) {
    coefKernelCyl<<<gridSize, BLOCK_SIZE>>>
      (toKernel(dN_coef), toKernel(t_d), toKernel(m_d), toKernel(p_d),
       toKernel(X_d), toKernel(Y_d), toKernel(iX_d), toKernel(iY_d),
       stride, m, ncylorder, lohi);
    
				// Begin the reduction per grid block
				//
    reduceSum<cuFP_t, BLOCK_SIZE><<<gridSize, BLOCK_SIZE, sMemSize>>>
      (toKernel(dc_coef), toKernel(dN_coef), osize, N);
      
				// Finish the reduction for this order
				// in parallel
    thrust::reduce_by_key
      (
       thrust::make_transform_iterator(index_begin, key_functor(gridSize)),
       thrust::make_transform_iterator(index_end,   key_functor(gridSize)),
       dc_coef.begin(), thrust::make_discard_iterator(), df_coef.begin()
       );
    
    thrust::host_vector<cuFP_t> ret = df_coef;
    for (size_t j=0; j<ncylorder; j++) {
      host_coefs[Imn(m, 'c', j, ncylorder)] = ret[2*j];
      if (m>0) host_coefs[Imn(m, 's', j, ncylorder)] = ret[2*j+1];
    }
  }

  // DEBUG
  //
  if (true) {
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

    int i = Imn(0, 'c', 0, ncylorder);
    auto cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(0, 'c', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(0, n, 'c');
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << i
		<< std::setw(20) << a
		<< std::setw(20) << b
		<< std::setw(20) << a - b
		<< std::setw(20) << (a - b)/fabs(*cmax)
		<< std::endl;
    }

    std::cout << "M=1c coefficients" << std::endl;

    i = Imn(1, 'c', 0, ncylorder);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(1, 'c', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(1, n, 'c');
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << i
		<< std::setw(20) << a
		<< std::setw(20) << b
		<< std::setw(20) << a - b
		<< std::setw(20) << (a - b)/fabs(*cmax)
		<< std::endl;
    }

    std::cout << "M=1s coefficients" << std::endl;

    i = Imn(1, 's', 0, ncylorder);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(1, 's', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(1, n, 's');
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << i
		<< std::setw(20) << a
		<< std::setw(20) << b
		<< std::setw(20) << a - b
		<< std::setw(20) << (a - b)/fabs(*cmax)
		<< std::endl;
    }

    std::cout << "M=2c coefficients" << std::endl;

    i = Imn(2, 'c', 0, ncylorder);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(2, 'c', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(2, n, 'c');
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << i
		<< std::setw(20) << a
		<< std::setw(20) << b
		<< std::setw(20) << a - b
		<< std::setw(20) << (a - b)/fabs(*cmax)
		<< std::endl;
    }
    
    std::cout << "M=2s coefficients" << std::endl;

    i = Imn(2, 's', 0, ncylorder);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+ncylorder, LessAbs<cuFP_t>());

    for (size_t n=0; n<ncylorder; n++) {
      int    i = Imn(2, 's', n, ncylorder);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(2, n, 's');
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << i
		<< std::setw(20) << a
		<< std::setw(20) << b
		<< std::setw(20) << a - b
		<< std::setw(20) << (a - b)/fabs(*cmax)
		<< std::endl;
    }

    std::cout << std::string(2*4+4*20, '-') << std::endl;
  }


  //
  // TEST comparison of coefficients for debugging
  //
  if (true) {

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

  // Compute number and total mass of particles used in coefficient
  // determination
  //
  thrust::sort(m_d.begin(), m_d.end());

  auto m_it   = thrust::upper_bound(m_d.begin(), m_d.end(), 0.0);
  use[0]      = thrust::distance(m_it, m_d.end());
  cylmass0[0] = thrust::reduce  (m_it, m_d.end());
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

  // Sort particles and get coefficient size
  //
  PII lohi = cC->CudaSortByLevel(mlevel, multistep);

  // Compute grid
  //
  unsigned int N         = lohi.second - lohi.first;
  unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  unsigned int gridSize  = N/BLOCK_SIZE/stride;

  if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

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

#ifdef VERBOSE
  std::cout << std::endl << "**" << std::endl
	    << "** N      = " << N          << std::endl
	    << "** Stride = " << stride     << std::endl
	    << "** Block  = " << BLOCK_SIZE << std::endl
	    << "** Grid   = " << gridSize   << std::endl
	    << "** Xcen   = " << ctr[0]     << std::endl
	    << "** Ycen   = " << ctr[1]     << std::endl
	    << "** Zcen   = " << ctr[2]     << std::endl
	    << "**" << std::endl;
#endif

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);

  // Maximum radius on grid
  //
  cuFP_t rmax = rcylmax * acyl;

  // Do the work
  //
  forceKernelCyl<<<gridSize, BLOCK_SIZE, sMemSize>>>
    (toKernel(cC->cuda_particles), toKernel(dev_coefs), toKernel(t_d),
     stride, mmax, ncylorder, lohi, rmax, cylmass, use_external);
}

void Cylinder::HtoD_coefs()
{
  // Check size
  host_coefs.resize((2*mmax+1)*ncylorder);

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

void Cylinder::destroy_cuda()
{
  // std::cout << "texture object array size = " << tex.size() << std::endl;
  for (size_t i=0; i<tex.size(); i++) {
    std::ostringstream sout;
    sout << "trying to free TextureObject [" << i << "]";
    cuda_safe_call(cudaDestroyTextureObject(tex[i]),
		   __FILE__, __LINE__, sout.str());
  }

  // std::cout << "cuInterpArray size = " << cuInterpArray.size() << std::endl;
  for (size_t i=0; i<cuInterpArray.size(); i++) {
    std::ostringstream sout;
    sout << "trying to free cuPitch [" << i << "]";
    cuda_safe_call(cudaFree(cuInterpArray[i]),
		     __FILE__, __LINE__, sout.str());
  }
    
  std::cout << "cuda memory freed" << std::endl;
}

void Cylinder::host_dev_force_compare()
{
  // Copy from device
  cC->host_particles = cC->cuda_particles;
  
  std::streamsize ss = std::cout.precision();
  std::cout.precision(10);

  std::cout << std::string(16+20*8, '-') << std::endl
	    << std::setw(8)  << "Index"  << std::setw(8)  << "Level"
	    << std::setw(20) << "ax [d]" << std::setw(20) << "ay [d]"
	    << std::setw(20) << "az [d]" << std::setw(20) << "ax [h]"
	    << std::setw(20) << "ay [h]" << std::setw(20) << "az [h]"
	    << std::setw(20) << "|Del a|/|a|"
	    << std::setw(20) << "|a|"    << std::endl;
  
  // Compare first and last 5 from the device list
  //
  for (size_t i=0; i<5; i++) 
    {
      auto indx = cC->host_particles[i].indx;
      auto levl = cC->host_particles[i].level;
      
      std::cout << std::setw(8) << indx	<< std::setw(8) << levl;

      for (int k=0; k<3; k++)
	std::cout << std::setw(20) << cC->host_particles[i].acc[k];

      for (int k=0; k<3; k++)
	std::cout << std::setw(20) << cC->Particles()[indx].acc[k];

      double diff = 0.0, norm = 0.0;
      for (int k=0; k<3; k++) {
	double b  = cC->host_particles[i].acc[k];
	double a  = cC->Particles()[indx].acc[k];
	diff += (a - b)*(a - b);
	norm += a*a;
      }
      std::cout << std::setw(20) << sqrt(diff/norm)
		<< std::setw(20) << sqrt(norm) << std::endl;
    }
  
  for (size_t j=0; j<5; j++) 
    {
      size_t i = cC->host_particles.size() - 5 + j;

      auto indx = cC->host_particles[i].indx;
      auto levl = cC->host_particles[i].level;

      std::cout << std::setw(8) << indx	<< std::setw(8) << levl;
      
      for (int k=0; k<3; k++)
	std::cout << std::setw(20) << cC->host_particles[i].acc[k];

      for (int k=0; k<3; k++)
	std::cout << std::setw(20) << cC->Particles()[indx].acc[k];

      double diff = 0.0, norm = 0.0;
      for (int k=0; k<3; k++) {
	double b  = cC->host_particles[i].acc[k];
	double a  = cC->Particles()[indx].acc[k];
	diff += (a - b)*(a - b);
	norm += a*a;
      }
      std::cout << std::setw(20) << sqrt(diff/norm)
		<< std::setw(20) << sqrt(norm) << std::endl;
    }

  std::cout << std::string(16+20*8, '-') << std::endl;
  std::cout.precision(ss);
}
