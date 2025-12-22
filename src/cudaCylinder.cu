// -*- C++ -*-

#include <Component.H>
#include <Cylinder.H>
#include <cudaReduce.cuH>
#include <float.h>
#include "expand.H"

// Define for debugging
//
// #define OFF_GRID_ALERT
// #define BOUNDS_CHECK
// #define VERBOSE_CTR
// #define VERBOSE_DBG

// Global symbols for coordinate transformation
//
__device__ __constant__
cuFP_t cylRscale, cylHscale, cylXmin, cylXmax, cylYmin, cylYmax, cylDxi, cylDyi, cylCen[3], cylBody[9], cylOrig[9];

__device__ __constant__
int cylNumx, cylNumy, cylCmapR, cylCmapZ, cylOrient;

__device__ __constant__
bool cylAcov;

// Index function for sine and cosine coefficients
//
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

// Index function for modulus coefficients
//
__host__ __device__
int Jmn(int m, int n, int nmax)
{
  return m*nmax + n;
}

// Index function for covariance matrix elements
//
__host__ __device__
int Kmn(int m, int n, int o, int nmax)
{
  return (m*nmax + n)*nmax + o;
}

__global__
void testConstantsCyl()
{
  printf("-------------------------\n");
  printf("---Cylindrical constants-\n");
  printf("-------------------------\n");
  printf("   Rscale = %f\n", cylRscale);
  printf("   Hscale = %f\n", cylHscale);
  printf("   Xmin   = %f\n", cylXmin  );
  printf("   Xmax   = %f\n", cylXmax  );
  printf("   Ymin   = %f\n", cylYmin  );
  printf("   Ymax   = %f\n", cylYmax  );
  printf("   Dxi    = %f\n", cylDxi   );
  printf("   Dyi    = %f\n", cylDyi   );
  printf("   Numx   = %d\n", cylNumx  );
  printf("   Numy   = %d\n", cylNumy  );
  printf("   CmapR  = %d\n", cylCmapR );
  printf("   CmapZ  = %d\n", cylCmapZ );
  if (cylAcov)
    printf("   cylAcov: true\n"       );
  else
    printf("   cylAcov: false\n"      );
  printf("-------------------------\n");
}

// R coordinate transformation
//
__device__
cuFP_t cu_r_to_xi_cyl(cuFP_t r)
{
  cuFP_t ret;

  if (cylCmapR==1) {
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

  if (cylCmapR==1) {
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

  if (cylCmapR==1) {
    ret = 0.5*(1.0 - xi)*(1.0 - xi) / cylRscale;
  } else {
    ret = 1.0;
  }

  return ret;
}

// Z coordinate transformation
//
__device__
cuFP_t cu_z_to_y_cyl(cuFP_t z)
{
  cuFP_t ret;

  if (cylCmapZ==1)
    ret = z/(fabs(z)+FLT_MIN)*asinh(fabs(z/cylHscale));
  else if (cylCmapZ==2)
    return z/sqrt(z*z + cylHscale*cylHscale);
  else
    ret = z;

  return ret;
}

__device__
cuFP_t cu_y_to_z_cyl(cuFP_t y)
{
  cuFP_t ret;

  if (cylCmapZ==1)
    ret = cylHscale*sinh(y);
  else if (cylCmapZ==2)
    ret = y * cylHscale/sqrt(1.0 - y*y);
  else
    ret = y;

  return ret;
}


__device__
cuFP_t cu_d_y_to_z_cyl(cuFP_t y)
{
  cuFP_t ret;

  if (cylCmapZ==1)
    ret = cylHscale*cosh(y);
  else if (cylCmapZ==2)
    return cylHscale*pow(1.0-y*y, -1.5);
  else
    ret = 1.0;

  return ret;
}


// Initialize for streams
//
void Cylinder::cuda_initialize()
{
  // Nothing so far
}

// Copy constants to device
//
void Cylinder::initialize_mapping_constants()
{
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

  cuda_safe_call(cudaMemcpyToSymbol(cylCmapR,  &f.cmapR,  sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylCmapR");

  cuda_safe_call(cudaMemcpyToSymbol(cylCmapZ,  &f.cmapZ,  sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylCmapZ");

  cuda_safe_call(cudaMemcpyToSymbol(cylAcov,   &subsamp,  sizeof(bool),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylAcov");
}


__global__ void coordKernelCyl
(dArray<cudaParticle> P, dArray<int> I,
 dArray<cuFP_t> mass, dArray<cuFP_t> phi,
 dArray<cuFP_t> Xfac, dArray<cuFP_t> Yfac,
 dArray<int> IndX, dArray<int> IndY,
 unsigned int stride, PII lohi, cuFP_t rmax)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Particle counter
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {	// Is particle index in range?

#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
    
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
	if (indX>cylNumx-1) indX = cylNumx - 1;
	
	if (indY<0) indY = 0;
	if (indY>cylNumy-1) indY = cylNumy - 1;
	
	Xfac._v[i] = cuFP_t(indX+1) - X;
	IndX._v[i] = indX;

	Yfac._v[i] = cuFP_t(indY+1) - Y;
	IndY._v[i] = indY;

#ifdef OFF_GRID_ALERT
	if (Xfac._v[i]<-0.5 or Xfac._v[i]>1.5)
	  printf("cuda cylinder--X off grid: x=%f R=%f\n", X, R);
	if (Yfac._v[i]<-0.5 or Yfac._v[i]>1.5)
	  printf("cuda cylinder--Y off grid: y=%f z=%f\n", Y, zz);
#endif
#ifdef BOUNDS_CHECK
	if (i>=Xfac._s) printf("cuda cylinder out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=IndX._s) printf("cuda cylinder out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=Yfac._s) printf("cuda cylinder out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=IndY._s) printf("cuda cylinder out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      }
    }
  }
}


__global__ void coefKernelCyl
(dArray<cuFP_t> coef, dArray<cuFP_t> tvar,
 dArray<thrust::complex<cuFP_t>> work,
 dArray<cuFP_t> used, dArray<cudaTextureObject_t> tex,
 dArray<cuFP_t> Mass, dArray<cuFP_t> Phi,
 dArray<cuFP_t> Xfac, dArray<cuFP_t> Yfac,
 dArray<int> indX, dArray<int> indY,
 int stride, int m, unsigned int nmax, PII lohi, bool compute)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int N   = lohi.second - lohi.first;

  const cuFP_t norm = -4.0*M_PI;    // Biorthogonality factor

  for (int n=0; n<stride; n++) {

    // Particle counter
    //
    int i     = tid*stride + n;
    int npart = i + lohi.first;

    if (npart < lohi.second) {	// Check that particle index is in
				// range for consistency with other
				// kernels
      cuFP_t mass = Mass._v[i];
      
#ifdef BOUNDS_CHECK
      if (i>=Mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif      
      if (mass>0.0) {
				// For accumulating mass of used particles
	if (m==0) used._v[i] = mass;

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
	  cuFP_t d00   = tex3D<float>(tex._v[k], indx,   indy  , 0);
	  cuFP_t d10   = tex3D<float>(tex._v[k], indx+1, indy  , 0);
	  cuFP_t d01   = tex3D<float>(tex._v[k], indx,   indy+1, 0);
	  cuFP_t d11   = tex3D<float>(tex._v[k], indx+1, indy+1, 0);

#else
	  cuFP_t d00   = int2_as_double(tex3D<int2>(tex._v[k], indx,   indy  , 0));
	  cuFP_t d10   = int2_as_double(tex3D<int2>(tex._v[k], indx+1, indy  , 0));
	  cuFP_t d01   = int2_as_double(tex3D<int2>(tex._v[k], indx,   indy+1, 0));
	  cuFP_t d11   = int2_as_double(tex3D<int2>(tex._v[k], indx+1, indy+1, 0));
#endif
	  cuFP_t valC  = c00*d00 + c10*d10 + c01*d01 + c11*d11, valS = 0.0;
	  
#ifdef BOUNDS_CHECK
	  if (k>=tex._s)            printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	  if ((2*n+0)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	  coef._v[(2*n+0)*N + i] = valC * cosp * norm * mass;

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
	    valS = c00*d00 + c10*d10 + c01*d01 + c11*d11;
	    coef._v[(2*n+1)*N + i] = valS * sinp * norm * mass;

#ifdef BOUNDS_CHECK
	    if ((2*n+1)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	  }
	  // m==0
	  else {
	    coef._v[(2*n+1)*N + i] = 0.0;
	  }

	  if (compute and tvar._s>0) {
	    cuFP_t val = sqrt(valC*valC*cosp*cosp + valS*valS*sinp*sinp);
	    if (cylAcov) tvar._v[n*N + i   ] = val * norm * mass;
	    else         work._v[i*nmax + n] = thrust::complex<cuFP_t>(valC*cosp + valS*sinp, valC*sinp - valS*cosp) * norm;
	  }
	  
	}
	// END: norder loop

	if (compute and not cylAcov and tvar._s>0) {

	  // Variance
	  int c = 0;
	  for (int r=0; r<nmax; r++) {
	    for (int s=r; s<nmax; s++) {
	      tvar._v[N*c + i] =
		(work._v[i*nmax + r] * thrust::conj(work._v[i*nmax + s])).real() * mass;
	      c++;
	    }
	  }
	  // Mean
	  for (int r=0; r<nmax; r++) {
	    tvar._v[N*c + i] = work._v[i*nmax + r].real() * mass;
	    c++;
	  }
	}

      } else {
	// No contribution from off-grid particles
	for (int n=0; n<nmax; n++) {
	  coef._v[(2*n+0)*N + i] = 0.0;
	  if (m) coef._v[(2*n+1)*N + i] = 0.0;
	}

	if (compute and tvar._s>0) {
	  if (cylAcov) {
	    for (int n=0; n<nmax; n++) {
	      tvar._v[n*N + i] = 0.0;
	    }
	  } else {
	    int c = 0;
	    for (int r=0; r<nmax; r++) {
	      for (int s=r; s<nmax; s++) {
		tvar._v[N*c + i] = 0.0;
		c++;
	      }
	    }
	    for (int r=0; r<nmax; r++) {
	      tvar._v[N*c + i] = 0.0;
	      c++;
	    }
	  }
	}

      } // mass value check

    } // particle index check

  } // stride loop

}

__global__ void
forceKernelCyl(dArray<cudaParticle> P, dArray<int> I,
	       dArray<cuFP_t> coef,
	       dArray<cudaTextureObject_t> tex,
	       int stride, unsigned int mmax, unsigned int mlim,
	       unsigned int nmax, PII lohi,
	       cuFP_t rmax, cuFP_t cylmass)
{
  // Thread ID
  //
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

  // Maximum radius squared
  //
  const cuFP_t rmax2 = rmax*rmax;

  // Algorithm constants
  //
  constexpr cuFP_t ratmin = 0.75;
  constexpr cuFP_t maxerf = 3.0;
  constexpr cuFP_t midpt  = ratmin + 0.5*(1.0 - ratmin);
  constexpr cuFP_t rsmth  = 0.5*(1.0 - ratmin)/maxerf;

  int muse = mmax > mlim ? mlim : mmax;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {	// Check that particle index is in
				// range
      
#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      
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
      
      cuFP_t ratio = sqrt( (R2 + zz*zz)/rmax2 );
      cuFP_t mfactor = 1.0, frac = 1.0, cfrac = 0.0;

      if (ratio >= 1.0) {
	frac  = 0.0;
	cfrac = 1.0;
      } else if (ratio > ratmin) {
	frac  = 0.5*(1.0 - erf( (ratio - midpt)/rsmth )) * mfactor;
	cfrac = 1.0 - frac;
      } else {
	frac  = 1.0;
      }

      // mfactor will apply this a fraction of this component's force
      // when mixture models are implemented (see Cylinder.cc)
      /*
      cfrac *= mfactor;
      frac  *= mfactor;
      */
	
      cuFP_t fr = 0.0;
      cuFP_t fz = 0.0;
      cuFP_t fp = 0.0;
      cuFP_t pp = 0.0;
      cuFP_t pa = 0.0;
      
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
	if (delx0<-0.5 or delx0>1.5) // X value check
	  printf("cuda cylinder--X off grid: x=%f [%d, %d] R=%f\n", delx0, indX, indY, R);

	if (dely0<-0.5 or dely0>1.5) // Y value check
	  printf("cuda cylinder--Y off grid: y=%f [%d, %d] z=%f\n", dely0, indX, indY, zz);
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

	for (int mm=0; mm<=muse; mm++) {

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
	pa     += pp * frac;
      }

      if (ratio > ratmin) {

	cuFP_t r3 = R2 + zz*zz;
	pp = -cylmass/sqrt(r3);	// -M/r
	fr = pp/r3;		// -M/r^3

	acc[0] += xx*fr * cfrac;
	acc[1] += yy*fr * cfrac;
	acc[2] += zz*fr * cfrac;
	pa     += pp    * cfrac;
      }

      if (cylOrient) {
	for (int j=0; j<3; j++) {
	  for (int k=0; k<3; k++) p.acc[j] += cylOrig[3*j+k]*acc[k];
	}
      } else {
	for (int j=0; j<3; j++) p.acc[j] += acc[j];
      }

      p.pot += pa;

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


void Cylinder::cudaStorage::resize_coefs
(int ncylorder, int mmax, int N, int gridSize, int stride,
 int sampT, bool pcavar, bool pcaeof, bool subsamp)
{
  // Reserve space for coefficient reduction
  //
  if (dN_coef.capacity() < 2*ncylorder*N)
    dN_coef.reserve(2*ncylorder*N);
  
  if (dc_coef.capacity() < 2*ncylorder*gridSize)
    dc_coef.reserve(2*ncylorder*gridSize);
  
  if (m_d .capacity() < N) m_d .reserve(N);
  if (u_d .capacity() < N) u_d .reserve(N);
  if (X_d .capacity() < N) X_d .reserve(N);
  if (Y_d .capacity() < N) Y_d .reserve(N);
  if (p_d .capacity() < N) p_d .reserve(N);
  if (iX_d.capacity() < N) iX_d.reserve(N);
  if (iY_d.capacity() < N) iY_d.reserve(N);
  
  
  // Fixed size arrays
  //
  if (pcavar) {
    T_covr.resize(sampT);
    for (int T=0; T<sampT; T++) {
      if (subsamp)
	T_covr[T].resize((mmax+1)*ncylorder);
      else
	T_covr[T].resize((mmax+1)*ncylorder*(ncylorder+3)/2);
    }
  }

  if (pcaeof or pcavar) {
    if (subsamp) {
      dN_tvar.resize(ncylorder*N);
      dc_tvar.resize(ncylorder*gridSize);
      dw_tvar.resize(ncylorder);
    } else {
      int csz = ncylorder*(ncylorder+3)/2;
      dN_tvar.resize(csz*N);
      dW_tvar.resize(ncylorder*gridSize*BLOCK_SIZE*stride);
      dc_tvar.resize(csz*gridSize);
      dw_tvar.resize(csz);
    }
  }
  
  // Set space for current step
  //
  dN_coef.resize(2*ncylorder*N);
  dc_coef.resize(2*ncylorder*gridSize);
  dw_coef.resize(2*ncylorder);	// This will stay fixed

  // Space for coordinate arrays on the current step
  //
  m_d .resize(N);
  u_d .resize(N);
  X_d .resize(N);
  Y_d .resize(N);
  p_d .resize(N);
  iX_d.resize(N);
  iY_d.resize(N);
}

void Cylinder::cuda_zero_coefs()
{
  auto cr = component->cuStream;
  
  // Resize output array
  //
  cuS.df_coef.resize((mmax+1)*2*nmax);
    
  // Zero output array
  //
  thrust::fill(thrust::cuda::par.on(cr->stream),
	       cuS.df_coef.begin(), cuS.df_coef.end(), 0.0);

  // Resize and zero PCA arrays
  //
  if (pcavar or (nint>0)) {
				// (Re)initialize?
    if (cuS.T_covr.size() != sampT) {
      cuS.T_covr.resize(sampT);
      for (int T=0; T<sampT; T++) {
	if (subsamp)
	  cuS.T_covr[T].resize((mmax+1)*nmax);
	else
	  cuS.T_covr[T].resize((mmax+1)*nmax*(nmax+3)/2);
      }
    }
    
    for (int T=0; T<sampT; T++) {
      thrust::fill(thrust::cuda::par.on(cr->stream),
		   cuS.T_covr[T].begin(), cuS.T_covr[T].end(), 0.0);
    }
  }

  if (pcaeof) {
    
    cuS.df_tvar.resize((mmax+1)*nmax*(nmax+1)/2);
    
    thrust::fill(thrust::cuda::par.on(cr->stream),
		 cuS.df_tvar.begin(), cuS.df_tvar.end(), 0.0);
  }
}

void Cylinder::determine_coefficients_cuda()
{
  // Only do this once but copying mapping coefficients and textures
  // must be done every time
  //
  if (initialize_cuda_cyl) {
    initialize_cuda();
    initialize_cuda_cyl = false;
  }

  // Copy coordinate mapping
  //
  initialize_mapping_constants();

  // Copy texture memory
  //
  t_d = tex;

  std::cout << std::scientific;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, component->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  // This will stay fixed for the entire run
  //
				// Sine and cosine components
  host_coefs.resize((2*mmax+1)*nmax);
				// Variance components
  host_covar.resize((mmax+1)*nmax*nmax);

  if (pcavar or (nint>0)) {	// Set sample size
    if (defSampT) sampT = defSampT;
    else          sampT = floor(sqrt(component->CurTotal()));
    host_coefsT.resize(sampT);	// Modulus components
    host_covarT.resize(sampT);	// Upper diagonal
    for (int T=0; T<sampT; T++) {
      host_coefsT[T].resize((mmax+1)*nmax);
      host_covarT[T].resize((mmax+1)*nmax*nmax);
    }
    host_massT.resize(sampT);
    host_numbT.resize(sampT);
  }

  // Set component center and orientation
  //
  std::vector<cuFP_t> ctr;
  for (auto v : component->getCenter(Component::Local | Component::Centered)) ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(cylCen, &ctr[0], sizeof(cuFP_t)*3, size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylCen");

  bool orient = (component->EJ & Orient::AXIS) && !component->EJdryrun;

  int tmp = orient ? 1 : 0;
  cuda_safe_call(cudaMemcpyToSymbol(cylOrient, &tmp, sizeof(int), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylOrient");

  if (orient) {
    std::vector<cuFP_t> trans(9);
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++) trans[i*3+j] = component->orient->transformBody()(i, j);
  
    cuda_safe_call(cudaMemcpyToSymbol(cylBody, &trans[0], sizeof(cuFP_t)*9, size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying cylBody");
  }

  // Get the stream for this component
  //
  auto cs = component->cuStream;

  // VERBOSE diagnostic output on first call
  //
  static bool firstime = true;

  if (firstime and myid==0 and VERBOSE>4) {
    testConstantsCyl<<<1, 1, 0, cs->stream>>>();
    cudaDeviceSynchronize();
    cuda_check_last_error_mpi("cudaDeviceSynchronize", __FILE__, __LINE__, myid);
    firstime = false;
  }
  
  // Zero counter and coefficients
  //
  thrust::fill(host_coefs.begin(), host_coefs.end(), 0.0);

  if (pcavar or (nint>0)) {
    for (int T=0; T<sampT; T++) {
      thrust::fill(host_coefsT[T].begin(), host_coefsT[T].end(), 0.0);
      thrust::fill(host_covarT[T].begin(), host_covarT[T].end(), 0.0);
    }
    thrust::fill(host_massT.begin(), host_massT.end(), 0.0);
    thrust::fill(host_numbT.begin(), host_numbT.end(), 0  );
  }

  // Zero out coefficient storage
  //
  cuda_zero_coefs();

  // Maximum radius on grid; get actual value from EmpCylSL
  //
  cuFP_t rmax = ortho->get_ascale()*ortho->get_rtable();
  // Get sorted particle range for mlevel
  //
  PII lohi = component->CudaGetLevelRange(mlevel, mlevel), cur;

  if (false) {
    for (int n=0; n<numprocs; n++) {
      if (myid==n) std::cout << "[" << myid << "] mlevel=" << mlevel
			     << " coef check (lo, hi) = (" << lohi.first << ", "
			     << lohi.second << ")" << std::endl
			     << std::string(60, '-') << std::endl;
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
  
  unsigned int Ntotal = lohi.second - lohi.first;
  unsigned int Npacks = Ntotal/component->bunchSize + 1;

  // Loop over bunches
  //
  for (int n=0; n<Npacks; n++) {

    // Current bunch
    //
    cur. first = lohi.first + component->bunchSize*n;
    cur.second = lohi.first + component->bunchSize*(n+1);
    cur.second = std::min<unsigned int>(cur.second, lohi.second);

    if (cur.second <= cur.first) break;
    
    // Compute grid
    //
    unsigned int N         = cur.second - cur.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
  
    if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

#ifdef VERBOSE_CTR
    static unsigned debug_max_count = 100;
    static unsigned debug_cur_count = 0;
    if (debug_cur_count++ < debug_max_count) {
      std::cout << std::endl
		<< "** -------------------------" << std::endl
		<< "** cudaCylinder coefficients" << std::endl
		<< "** -------------------------" << std::endl
		<< "** N      = " << N            << std::endl
		<< "** Npacks = " << Npacks       << std::endl
		<< "** I low  = " << cur.first    << std::endl
		<< "** I high = " << cur.second   << std::endl
		<< "** Stride = " << stride       << std::endl
		<< "** Block  = " << BLOCK_SIZE   << std::endl
		<< "** Grid   = " << gridSize     << std::endl
		<< "** Xcen   = " << ctr[0]       << std::endl
		<< "** Ycen   = " << ctr[1]       << std::endl
		<< "** Zcen   = " << ctr[2]       << std::endl
		<< "** Level  = " << mlevel       << std::endl
		<< "** lo     = " << lohi.first   << std::endl
		<< "** hi     = " << lohi.second  << std::endl
		<< "**" << std::endl;
  }
#endif
  
    // Adjust cached storage, if necessary
    //
    cuS.resize_coefs(nmax, mmax, N, gridSize, stride,
		     sampT, pcavar or (nint>0), pcaeof, subsamp);
    
    // Shared memory size for the reduction
    //
    int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
    
    // Compute the coordinate transformation
    // 
    coordKernelCyl<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
      (toKernel(cs->cuda_particles), toKernel(cs->indx1),
       toKernel(cuS.m_d), toKernel(cuS.p_d),
       toKernel(cuS.X_d), toKernel(cuS.Y_d), toKernel(cuS.iX_d),
       toKernel(cuS.iY_d), stride, cur, rmax);
      
    // Compute the coefficient contribution for each order
    //
    int psize = nmax;
    int osize = nmax*2;
    int vsize = nmax*(nmax+1)/2 + nmax;
    auto beg  = cuS.df_coef.begin();
    auto begV = cuS.df_tvar.begin();
    std::vector<thrust::device_vector<cuFP_t>::iterator> bg, bh;

    if (pcavar or (nint > 0)) {
      for (int T=0; T<sampT; T++) {
	bh.push_back(cuS.T_covr[T].begin());
      }
    }

    thrust::fill(cuS.u_d.begin(), cuS.u_d.end(), 0.0);

    for (int m=0; m<=mmax; m++) {

      coefKernelCyl<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
	(toKernel(cuS.dN_coef), toKernel(cuS.dN_tvar), toKernel(cuS.dW_tvar),
	 toKernel(cuS.u_d), toKernel(t_d), toKernel(cuS.m_d), toKernel(cuS.p_d),
	 toKernel(cuS.X_d), toKernel(cuS.Y_d), toKernel(cuS.iX_d), toKernel(cuS.iY_d),
	 stride, m, nmax, cur, compute);
      
      // Begin the reduction by blocks [perhaps this should use a
      // stride?]
      //
      unsigned int gridSize1 = N/BLOCK_SIZE;
      if (N > gridSize1*BLOCK_SIZE) gridSize1++;

      reduceSum<cuFP_t, BLOCK_SIZE>
	<<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
	(toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), osize, N);
      
      // Finish the reduction for this order in parallel
      //
      thrust::counting_iterator<int> index_begin(0);
      thrust::counting_iterator<int> index_end(gridSize1*osize);

      // The key_functor indexes the sum reduced series by array index
      //
      thrust::reduce_by_key
	(
	 thrust::cuda::par.on(cs->stream),
	 thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
	 thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
	 cuS.dc_coef.begin(), thrust::make_discard_iterator(), cuS.dw_coef.begin()
	 );

      thrust::transform(thrust::cuda::par.on(cs->stream),
			cuS.dw_coef.begin(), cuS.dw_coef.end(),
			beg, beg, thrust::plus<cuFP_t>());

      thrust::advance(beg, osize);

      if (compute) {

	// Reuse dN_coef and use dN_tvar to create sampT partitions
	//
	if (pcavar or (nint>0)) {

	  int sN = N/sampT;
	  int nT = sampT;

	  if (sN==0) {	// Fail-safe underrun
	    sN = 1;
	    nT = N;
	  }

	  for (int T=0; T<nT; T++) {
	    int k = sN*T;	// Starting position
	    int s = sN;		// Size of bunch
				// Last bunch
	    if (T==sampT-1) s = N - k;
	    
	    // Begin the reduction per grid block; we need gridsize1
	    // blocks of the current block size
	    //
	    unsigned int gridSize1 = s/BLOCK_SIZE;
	    if (s > gridSize1*BLOCK_SIZE) gridSize1++;
	    
	    if (m==0) {
	      auto mbeg = cuS.u_d.begin();
	      auto mend = mbeg;
	      thrust::advance(mbeg, sN*T);
	      if (T<sampT-1) thrust::advance(mend, sN*(T+1));
	      else mend = cuS.u_d.end();
	      
	      host_massT[T] += thrust::reduce(mbeg, mend);
	      host_numbT[T] += mend - mbeg;
	    }

	    if (subsamp) {

	      // Variance
	      //
	      reduceSumS<cuFP_t, BLOCK_SIZE>
		<<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
		(toKernel(cuS.dc_tvar), toKernel(cuS.dN_tvar), psize, N, k, k+s);
      
	      // Finish the reduction for this order in parallel
	      //
	      thrust::counting_iterator<int> indx2_begin(0);
	      thrust::counting_iterator<int> indx2_end(gridSize1*psize);
	      
	      // The key_functor indexes the sum reduced series by array index
	      //
	      thrust::reduce_by_key
		(
		 thrust::cuda::par.on(cs->stream),
		 thrust::make_transform_iterator(indx2_begin, key_functor(gridSize1)),
		 thrust::make_transform_iterator(indx2_end,   key_functor(gridSize1)),
		 cuS.dc_tvar.begin(), thrust::make_discard_iterator(), cuS.dw_tvar.begin()
		 );
	    
	      thrust::transform(thrust::cuda::par.on(cs->stream),
				cuS.dw_tvar.begin(), cuS.dw_tvar.end(),
				bh[T], bh[T], thrust::plus<cuFP_t>());
	    
	      thrust::advance(bh[T], psize);
	      
	    } else {

	      // Variance
	      //
	      reduceSumS<cuFP_t, BLOCK_SIZE>
		<<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
		(toKernel(cuS.dc_tvar), toKernel(cuS.dN_tvar), vsize, N, k, k+s);
      
	      // Finish the reduction for this order in parallel
	      //
	      thrust::counting_iterator<int> indx2_begin(0);
	      thrust::counting_iterator<int> indx2_end(gridSize1*vsize);
	      
	      // The key_functor indexes the sum reduced series by array index
	      //
	      thrust::reduce_by_key
		(
		 thrust::cuda::par.on(cs->stream),
		 thrust::make_transform_iterator(indx2_begin, key_functor(gridSize1)),
		 thrust::make_transform_iterator(indx2_end,   key_functor(gridSize1)),
		 cuS.dc_tvar.begin(), thrust::make_discard_iterator(), cuS.dw_tvar.begin()
		 );
	    
	      thrust::transform(thrust::cuda::par.on(cs->stream),
				cuS.dw_tvar.begin(), cuS.dw_tvar.end(),
				bh[T], bh[T], thrust::plus<cuFP_t>());
	    
	      thrust::advance(bh[T], vsize);
	    }
	  }
	}
	// END: pcavar block

	// Reduce EOF variance using dN_tvar (which may be reused from
	// pcavar computation)
	//
	if (pcaeof) {
	  
	  reduceSum<cuFP_t, BLOCK_SIZE>
	    <<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
	    (toKernel(cuS.dc_tvar), toKernel(cuS.dN_tvar), vsize, N);
      
	  // Finish the reduction for this order in parallel
	  //
	  thrust::counting_iterator<int> index_begin(0);
	  thrust::counting_iterator<int> index_end(gridSize1*vsize);
	  
	  // The key_functor indexes the sum reduced series by array
	  // index
	  //
	  thrust::reduce_by_key
	    (
	     thrust::cuda::par.on(cs->stream),
	     thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
	     thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
	       cuS.dc_tvar.begin(), thrust::make_discard_iterator(), cuS.dw_tvar.begin()
	     );

	  thrust::transform(thrust::cuda::par.on(cs->stream),
			    cuS.dw_tvar.begin(), cuS.dw_tvar.end(),
			    begV, begV, thrust::plus<cuFP_t>());

	  thrust::advance(begV, vsize);
	}
      }
    }

    // Compute number and total mass of particles used in coefficient
    // determination
    //
    thrust::sort(thrust::cuda::par.on(cs->stream),
		 cuS.m_d.begin(), cuS.m_d.end());
    
    auto exec  = thrust::cuda::par.on(cs->stream);
    auto first = cuS.u_d.begin();
    auto last  = cuS.u_d.end();

				// Sort used masses
				//
    thrust::sort(exec, first, last);

				// Iterator will point to first
				// non-zero element
				//
    thrust::device_vector<cuFP_t>::iterator it;

    // Workaround for: https://github.com/NVIDIA/thrust/pull/1104
    //
    if (thrust_binary_search_workaround) {
      cudaStreamSynchronize(cs->stream);
      cuda_check_last_error_mpi("cudaStreamSynchronize", __FILE__, __LINE__, myid);
      it = thrust::upper_bound(first, last, 0.0);
    } else {
      it = thrust::upper_bound(exec, first, last, 0.0);
    }
      
				// Number of non-zero elements
				//
    used += thrust::distance(it, last);

				// Sum of mass
				// 
    cylmass0[0] += thrust::reduce(exec, it, last);
  }

  if (Ntotal == 0) {
    return;
  }

  // Accumulate the coefficients from the device to the host
  //
  thrust::host_vector<cuFP_t> ret = cuS.df_coef;
  int offst = 0;
  for (int m=0; m<=mmax; m++) {
    for (size_t j=0; j<nmax; j++) {
      host_coefs[Imn(m, 'c', j, nmax)] += ret[2*j+offst];
      if (m>0) host_coefs[Imn(m, 's', j, nmax)] += ret[2*j+1+offst];
    }
    offst += 2*nmax;
  }

  if (compute) {

    // Variance computation
    //
    if (pcavar or (nint>0)) {

      for (int T=0; T<sampT; T++) {

	if (subsamp) {

	  thrust::host_vector<cuFP_t> retV = cuS.T_covr[T];

	  int offst = 0;

	  for (int m=0; m<=mmax; m++) {

	    for (size_t j=0; j<nmax; j++) {
	      host_coefsT[T][Jmn(m, j, nmax)] += retV[j + offst];
	    }

	    offst += nmax;
	  }

	}
	// END: subsample variance
	else {
	  
	  thrust::host_vector<cuFP_t> retM = cuS.T_covr[T];
	  int vffst = 0;
	  for (int m=0; m<=mmax; m++) {

	    // Variance assignment
	    //
	    int c = 0;
	    for (size_t j=0; j<nmax; j++) {
	      for (size_t k=j; k<nmax; k++) {
		host_covarT[T][Kmn(m, j, k, nmax)] += retM[c + vffst];
		if (k!=j)
		  host_covarT[T][Kmn(m, k, j, nmax)] += retM[c + vffst];
		c++;
	      }
	    }
	    
	    // Mean assignment
	    //
	    for (size_t j=0; j<nmax; j++) {
	      host_coefsT[T][Jmn(m, j, nmax)] += retM[c + vffst];
	      c++;
	    }

	    if (myid==0 and c != nmax*(nmax+3)/2)
	      std::cout << "out of bounds: c=" << c << " != "
			<< nmax*(nmax+3)/2 << std::endl;

	    vffst += nmax*(nmax+3)/2;
	  }
	}
	// END: full pop variance
      }
    }
	
    // EOF variance computation
    //
    if (pcaeof) {
      thrust::host_vector<cuFP_t> retV = cuS.df_tvar;
      int csz = nmax*(nmax+1)/2;
      if (retV.size() == (mmax+1)*csz) {
	for (int m=0; m<=mmax; m++) {
	  int c = 0;
	  for (size_t j=0; j<nmax; j++) {
	    for (size_t k=j; k<nmax; k++) {
	      ortho->set_tvar(m, j, k) += retV[csz*m + c];
	      if (j!=k) ortho->set_tvar(m, k, j) += retV[csz*m + c];
	      c++;
	    }
	  }
	}
      }
    }
  }


  // DEBUG, only useful for CUDAtest branch
  //
  if (false) {
    constexpr bool compareC = false;

    if (compareC) {
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
    } else {
      std::cout << std::string(2*4+20, '-') << std::endl
		<< "---- Cylindrical "      << std::endl
		<< std::string(2*4+20, '-') << std::endl;
      std::cout << "M=0 coefficients" << std::endl
		<< std::setprecision(10);

      std::cout << std::setw(4)  << "n"
		<< std::setw(4)  << "i"
		<< std::setw(20) << "GPU"
		<< std::endl;
    }

    int i = Imn(0, 'c', 0, nmax);
    auto cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      int    i = Imn(0, 'c', n, nmax);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(0, n, 'c');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }

    std::cout << "M=1c coefficients" << std::endl;

    i = Imn(1, 'c', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      int    i = Imn(1, 'c', n, nmax);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(1, n, 'c');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }

    std::cout << "M=1s coefficients" << std::endl;

    i = Imn(1, 's', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      int    i = Imn(1, 's', n, nmax);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(1, n, 's');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }

    std::cout << "M=2c coefficients" << std::endl;

    i = Imn(2, 'c', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      int    i = Imn(2, 'c', n, nmax);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(2, n, 'c');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }
    
    std::cout << "M=2s coefficients" << std::endl;

    i = Imn(2, 's', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      int    i = Imn(2, 's', n, nmax);
      cuFP_t a = host_coefs[i];
      cuFP_t b = ortho->get_coef(2, n, 's');
      if (compareC) {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::setw(20) << b
		  << std::setw(20) << a - b
		  << std::setw(20) << (a - b)/fabs(*cmax)
		  << std::endl;
      } else {
	std::cout << std::setw(4)  << n
		  << std::setw(4)  << i
		  << std::setw(20) << a
		  << std::endl;
      }
    }

    std::cout << std::string(2*4+4*20, '-') << std::endl;
  }


  //
  // TEST comparison of coefficients for debugging
  //
  if (false) {

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
	for (int n=0; n<nmax; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 'c';
	  elem.d = ortho->get_coef(m, n, 'c');
	  elem.f = host_coefs[Imn(m, 'c', n, nmax)];
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);
	  
	  compare.insert(std::make_pair(test, elem));;
	    
	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 'c'
	      << std::setw( 5) << Imn(m, 'c', n, nmax)
	      << std::setw(14) << elem.d
	      << std::setw(14) << elem.f
	      << std::endl;
	}

      } else {
	for (int n=0; n<nmax; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 'c';
	  elem.d = ortho->get_coef(m, n, 'c');
	  elem.f = host_coefs[Imn(m, 'c', n, nmax)];

	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 'c'
	      << std::setw( 5) << Imn(m, 'c', n, nmax)
	      << std::setw(14) << elem.d
	      << std::setw(14) << elem.f
	      << std::endl;
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);

	  compare.insert(std::make_pair(test, elem));;
	}

	for (int n=0; n<nmax; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 's';
	  elem.d = ortho->get_coef(m, n, 's');
	  elem.f = host_coefs[Imn(m, 's', n, nmax)];

	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 's'
	      << std::setw( 5) << Imn(m, 's', n-1, nmax)
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

}


void Cylinder::determine_acceleration_cuda()
{
  // Only do this once but copying mapping coefficients and textures
  // must be done every time
  //
  if (initialize_cuda_cyl) {
    initialize_cuda();
    initialize_cuda_cyl = false;
  }

  // Copy coordinate mapping
  //
  initialize_mapping_constants();

  // Copy texture memory
  //
  t_d = tex;

  std::cout << std::scientific;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cC->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  auto cs = cC->cuStream;

  // Assign expansion center and orientation
  //
  std::vector<cuFP_t> ctr;
  for (auto v : component->getCenter(Component::Local | Component::Centered)) ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(cylCen, &ctr[0], sizeof(cuFP_t)*3, size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cylCen");

  bool orient = (component->EJ & Orient::AXIS) && !component->EJdryrun;

  if (orient) {
    std::vector<cuFP_t> trans(9);
    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++)
	trans[i*3+j] = component->orient->transformBody()(i, j);
  
    cuda_safe_call(cudaMemcpyToSymbol(cylBody, &trans[0], sizeof(cuFP_t)*9, size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying cylBody");

    for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++)
	trans[i*3+j] = component->orient->transformOrig()(i, j);
  
    cuda_safe_call(cudaMemcpyToSymbol(cylOrig, &trans[0], sizeof(cuFP_t)*9, size_t(0), cudaMemcpyHostToDevice),
		   __FILE__, __LINE__, "Error copying cylOrig");
  }

  // Get particle index range for levels [mlevel, multistep]
  //
  PII lohi = cC->CudaGetLevelRange(mlevel, multistep);

  // Compute grid
  //
  unsigned int N         = lohi.second - lohi.first;
  unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
  if (N>0) {

    if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

#ifdef VERBOSE_CTR
    static unsigned debug_max_count = 100;
    static unsigned debug_cur_count = 0;
    if (debug_cur_count++ < debug_max_count) {
      std::cout << std::endl
		<< "** -------------------------" << std::endl
		<< "** cudaCylinder acceleration" << std::endl
		<< "** -------------------------" << std::endl
		<< "** N      = " << N            << std::endl
		<< "** Stride = " << stride       << std::endl
		<< "** Block  = " << BLOCK_SIZE   << std::endl
		<< "** Grid   = " << gridSize     << std::endl
		<< "** Xcen   = " << ctr[0]       << std::endl
		<< "** Ycen   = " << ctr[1]       << std::endl
		<< "** Zcen   = " << ctr[2]       << std::endl
		<< "** Level  = " << mlevel       << std::endl
		<< "** lo     = " << lohi.first   << std::endl
		<< "** hi     = " << lohi.second  << std::endl
		<< "**" << std::endl;
    }
#endif
    
    // Shared memory size for the reduction
    //
    int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
      
    // Maximum radius on grid; get value from EmpCylSL
    //
    cuFP_t rmax = ortho->get_ascale()*ortho->get_rtable();
      
    // Do the work
    //
    forceKernelCyl<<<gridSize, BLOCK_SIZE, sMemSize, cs->stream>>>
      (toKernel(cs->cuda_particles), toKernel(cs->indx1),
       toKernel(dev_coefs), toKernel(t_d),
       stride, mmax, mlim, nmax, lohi, rmax, cylmass);
    
  }
}

void Cylinder::HtoD_coefs()
{
  // Check size
  host_coefs.resize((2*mmax+1)*nmax); // Should stay fixed, no reserve

  // Copy from EmpCylSL
  
  // m loop
  //
  for (int m=0; m<=mmax; m++) {
    
    // n loop
    //
    for (int n=0; n<nmax; n++) {
      host_coefs[Imn(m, 'c', n, nmax)] = ortho->get_coef(m, n, 'c');
      if (m>0) host_coefs[Imn(m, 's', n, nmax)] = ortho->get_coef(m, n, 's');
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
    for (int n=0; n<nmax; n++) {
      ortho->set_coef(M, m, n, 'c') = host_coefs[Imn(m, 'c', n, nmax)];
      if (m>0) ortho->set_coef(M, m, n, 's') = host_coefs[Imn(m, 's', n, nmax)];
    }
  }

  if (compute) {

    // T loop
    //
    for (int T=0; T<sampT; T++) {

      // Copy mass and number per sample T
      //
      ortho->set_massT(T) += host_massT[T];
      ortho->set_numbT(T) += host_numbT[T];

      // m loop
      //
      for (int m=0; m<=mmax; m++) {
	
	// n loop
	//
	for (int n=0; n<nmax; n++) {

	  if (pcaeof)
	    ortho->set_coefT(T, m, n) += host_coefsT[T][Jmn(m, n, nmax)];
	  else
	    *ortho->set_VC(T, m, n) += host_coefsT[T][Jmn(m, n, nmax)];

	  // o loop
	  for (int o=0; o<nmax; o++) {
	    if (pcaeof)
	      ortho->set_covrT(T, m, n, o) += host_covarT[T][Kmn(m, n, o, nmax)];
	    else
	      *ortho->set_MV(T, m, n, o) += host_covarT[T][Kmn(m, n, o, nmax)];
	  }
	}
      }
    }
  }

}

void Cylinder::multistep_update_cuda()
{
  if (not self_consistent) return;

  // Sanity check: don't redo a previous update
  //
  static double last_time = -1.0;
  if (last_time < 0.0) last_time = tnow;
  else if (tnow - last_time < 1.0e-18) return;
  last_time = tnow;

  // The plan: for the current active level search above and below for
  // particles for correction to coefficient matrix
  //

  // Sort the device vector by level changes
  //
  auto chg = component->CudaSortLevelChanges();

#ifdef CYL_UPDATE_TABLE
  multistep_add_debug(chg);
#endif

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, component->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);
  auto cs = component->cuStream;

  // Maximum radius on grid
  //
  cuFP_t rmax = rcylmax * acyl;

  // Step through all levels
  //
  for (int olev=mfirst[mdrft]; olev<=multistep; olev++) {

    for (int nlev=0; nlev<=multistep; nlev++) {

      if (olev == nlev) continue;

      // Get range of update block in particle index
      //
      unsigned int Ntotal = chg[olev][nlev].second - chg[olev][nlev].first;

      if (Ntotal==0) continue; // No particles [from, to]=[olev, nlev]

      unsigned int Npacks = Ntotal/component->bunchSize + 1;

      // Zero out coefficient storage
      //
      cuda_zero_coefs();

#ifdef VERBOSE_DBG
      std::cout << "[" << myid << ", " << tnow
		<< "] Adjust cylinder: Ntotal=" << Ntotal << " Npacks=" << Npacks
		<< " for (m, d)=(" << olev << ", " << nlev << ")"
		<< " beg(olev)=" << mfirst[mdrft] << " mdrft=" << mdrft
		<< std::endl;
#endif

      // Loop over bunches
      //
      for (int n=0; n<Npacks; n++) {

	PII cur;
	
	// Current bunch
	//
	cur. first = chg[olev][nlev].first + component->bunchSize*n;
	cur.second = chg[olev][nlev].first + component->bunchSize*(n+1);
	cur.second = std::min<unsigned int>(cur.second, chg[olev][nlev].second);

	if (cur.second <= cur.first) break;
    
	// Compute grid
	//
	unsigned int N         = cur.second - cur.first;
	unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
	unsigned int gridSize  = N/BLOCK_SIZE/stride;
	
	if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

    
	// Adjust cached storage, if necessary
	//
	cuS.resize_coefs(nmax, mmax, N, gridSize, stride,
			 sampT, pcavar or (nint>0), pcaeof, subsamp);
	
	// Shared memory size for the reduction
	//
	int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
    
	// Compute the coordinate transformation
	// 
	coordKernelCyl<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
	  (toKernel(cs->cuda_particles), toKernel(cs->indx2),
	   toKernel(cuS.m_d), toKernel(cuS.p_d),
	   toKernel(cuS.X_d), toKernel(cuS.Y_d), toKernel(cuS.iX_d),
	   toKernel(cuS.iY_d), stride, cur, rmax);
      
	// Compute the coefficient contribution for each order
	//
	int osize = nmax*2;
	auto beg  = cuS.df_coef.begin();

	thrust::fill(cuS.u_d.begin(), cuS.u_d.end(), 0.0);

	for (int m=0; m<=mmax; m++) {

	  coefKernelCyl<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
	    (toKernel(cuS.dN_coef), toKernel(cuS.dN_tvar), toKernel(cuS.dW_tvar),
	     toKernel(cuS.u_d), toKernel(t_d), toKernel(cuS.m_d), toKernel(cuS.p_d),
	     toKernel(cuS.X_d), toKernel(cuS.Y_d), toKernel(cuS.iX_d), toKernel(cuS.iY_d),
	     stride, m, nmax, cur, false);

	  // Begin the reduction per grid block [perhaps this should
	  // use a stride?]
	  //
	  unsigned int gridSize1 = N/BLOCK_SIZE;
	  if (N > gridSize1*BLOCK_SIZE) gridSize1++;

	  reduceSum<cuFP_t, BLOCK_SIZE>
	    <<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
	    (toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), osize, N);
	  
	  // Finish the reduction for this order in parallel
	  //
	  thrust::counting_iterator<int> index_begin(0);
	  thrust::counting_iterator<int> index_end(gridSize1*osize);

	  // The key_functor indexes the sum reduced series by array index
	  //
	  thrust::reduce_by_key
	    (
	     thrust::cuda::par.on(cs->stream),
	     thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
	     thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
	     cuS.dc_coef.begin(), thrust::make_discard_iterator(), cuS.dw_coef.begin()
	     );

	  thrust::transform(thrust::cuda::par.on(cs->stream),
			    cuS.dw_coef.begin(), cuS.dw_coef.end(),
			    beg, beg, thrust::plus<cuFP_t>());
	  
	  thrust::advance(beg, osize);
	}
	// END: m loop
      }
      // END: bunches

      // Accumulate the coefficients from the device to the host
      //
      thrust::host_vector<cuFP_t> ret = cuS.df_coef;

      // Decrement current level and increment new level using the
      // EmpCylSL update matricies
      //
      // m loop
      //
      for (int m=0, offst=0; m<=mmax; m++) {
	// n loop
	//
	for (int n=0; n<nmax; n++) {
	  ortho->differC1[0][olev](m, n) -= ret[2*n+offst];
	  ortho->differC1[0][nlev](m, n) += ret[2*n+offst];
	  if (m>0) {
	    ortho->differS1[0][olev](m, n) -= ret[2*n+1+offst];
	    ortho->differS1[0][nlev](m, n) += ret[2*n+1+offst];
	  }
	}
	// END: n loop

	// Increment update in coefficient array
	offst += 2*nmax;
      }
      // END: m loop
    }
    // DONE: Inner loop
  }
  // DONE: Outer loop
}


void Cylinder::destroy_cuda()
{
  // Deallocate the texture objects
  //
  for (size_t i=0; i<tex.size(); i++) {
    cuda_check_error(cudaDestroyTextureObject(tex[i]),
		     "cudaDestroyTextureObject", __FILE__, __LINE__);
  }

  // Deallocate the texture data arrays
  //
  for (size_t i=0; i<cuInterpArray.size(); i++) {
    cuda_check_error(cudaFreeArray(cuInterpArray[i]),
		     "cudaFreeArray", __FILE__, __LINE__);
  }
}
