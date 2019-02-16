// -*- C++ -*-

#include <tuple>
#include <list>

#include <Component.H>
#include <SphericalBasis.H>
#include <cudaReduce.cuH>

#include <boost/make_shared.hpp>

// Define for debugging
//
// #define OFF_GRID_ALERT
// #define BOUNDS_CHECK
// #define VERBOSE
// #define NAN_CHECK

// Global symbols for coordinate transformation in SphericalBasis
//
__device__ __constant__
cuFP_t sphScale, sphRscale, sphHscale, sphXmin, sphXmax, sphDxi, sphCen[3];

__device__ __constant__
int   sphNumr, sphCmap;

__host__ __device__
int Ilm(int l, int m)
{
  if (l==0) return 0;
  return l*(l+1)/2 + m;
}

__host__ __device__
int Ilmn(int l, int m, char cs, int n, int nmax)
{
  int ret = 0;

  if (l==0) ret = n;
  else if (m==0) ret = l*l*nmax + n;
  else ret = (l*l + 2*m - 1 + (cs=='s' ? 1 : 0))*nmax + n;

#ifdef BOUNDS_CHECK
  if (ret >= (l+1)*(l+1)*nmax) {
    printf("Ilmn oab: %4d %4d %4d [%4d : %4d : %4d]\n", l, m, n, ret, (l+1)*(l+1)*nmax, nmax);
  }
#endif

  return ret;
}

__host__ __device__
void legendre_v(int lmax, cuFP_t x, cuFP_t* p)
{
  cuFP_t fact, somx2, pll, pl1, pl2;

  p[0] = pll = 1.0f;
  if (lmax > 0) {
    somx2 = sqrt( (1.0f - x)*(1.0f + x) );
    fact = 1.0f;
    for (int m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[Ilm(m, m)] = pll;
      fact += 2.0f;
    }
  }
  
  for (int m=0; m<lmax; m++) {
    pl2 = p[Ilm(m, m)];
    p[Ilm(m+1, m)] = pl1 = x*(2*m+1)*pl2;
    for (int l=m+2; l<=lmax; l++) {
      p[Ilm(l, m)] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }
}


__host__ __device__
void legendre_v2(int lmax, cuFP_t x, cuFP_t* p, cuFP_t* dp)
{
  cuFP_t fact, somx2, pll, pl1, pl2;

  p[0] = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (int m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[Ilm(m, m)] = pll;
      fact += 2.0;
    }
  }

  for (int m=0; m<lmax; m++) {
    pl2 = p[Ilm(m, m)];
    p[Ilm(m+1, m)] = pl1 = x*(2*m+1)*pl2;
    for (int l=m+2; l<=lmax; l++) {
      p[Ilm(l, m)] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (1.0-fabs(x) < MINEPS) {
    if (x>0) x =   1.0 - MINEPS;
    else     x = -(1.0 - MINEPS);
  }

  somx2 = 1.0/(x*x - 1.0);
  dp[0] = 0.0;
  for (int l=1; l<=lmax; l++) {
    for (int m=0; m<l; m++)
      dp[Ilm(l, m)] = somx2*(x*l*p[Ilm(l, m)] - (l+m)*p[Ilm(l-1, m)]);
    dp[Ilm(l, l)] = somx2*x*l*p[Ilm(l, l)];
  }
}

template<typename Iterator, typename Pointer>
__global__
void reduceUseSph(Iterator first, Iterator last, Pointer result)
{
  auto it = thrust::lower_bound(thrust::cuda::par, first, last, 0.0);
  *result = thrust::distance(it, last);
}

__global__
void testConstants()
{
  printf("** Scale  = %f\n", sphScale);
  printf("** Rscale = %f\n", sphRscale);
  printf("** Xmin   = %f\n", sphXmin);
  printf("** Xmax   = %f\n", sphXmax);
  printf("** Dxi    = %f\n", sphDxi);
  printf("** Numr   = %d\n", sphNumr);
  printf("** Cmap   = %d\n", sphCmap);
}

__device__
cuFP_t cu_r_to_xi(cuFP_t r)
{
  cuFP_t ret;

  if (sphCmap==1) {
    ret = (r/sphRscale-1.0)/(r/sphRscale+1.0);
  } else if (sphCmap==2) {
    ret = log(r);
  } else {
    ret = r;
  }    

  return ret;
}
    
__device__
cuFP_t cu_xi_to_r(cuFP_t xi)
{
  cuFP_t ret;

  if (sphCmap==1) {
    ret = (1.0+xi)/(1.0 - xi) * sphRscale;
  } else if (sphCmap==2) {
    ret = exp(xi);
  } else {
    ret = xi;
  }

  return ret;
}

__device__
cuFP_t cu_d_xi_to_r(cuFP_t xi)
{
  cuFP_t ret;

  if (sphCmap==1) {
    ret = 0.5*(1.0-xi)*(1.0-xi)/sphRscale;
  } else if (sphCmap==2) {
    ret = exp(-xi);
  } else {
    ret = 1.0;
  }

  return ret;
}

void SphericalBasis::cuda_initialize()
{
  // Initialize for streams
  //
  cuRingData.resize(cuStreams);
  cuRing = boost::make_shared<cuRingType>(cuRingData);
}

void SphericalBasis::initialize_mapping_constants()
{
  // Copy constants to device
  //
  
  cudaMappingConstants f = getCudaMappingConstants();
  cuFP_t z;

  cuda_safe_call(cudaMemcpyToSymbol(sphScale, &(z=scale), sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphScale");

  cuda_safe_call(cudaMemcpyToSymbol(sphRscale, &f.rscale, sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphRscale");

  cuda_safe_call(cudaMemcpyToSymbol(sphXmin,   &f.xmin,   sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphXmin");

  cuda_safe_call(cudaMemcpyToSymbol(sphXmax,   &f.xmax,   sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphXmax");

  cuda_safe_call(cudaMemcpyToSymbol(sphDxi,    &f.dxi,    sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphDxi");

  cuda_safe_call(cudaMemcpyToSymbol(sphNumr,   &f.numr,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphuNumr");

  cuda_safe_call(cudaMemcpyToSymbol(sphCmap,   &f.cmap,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphCmap");
}


__global__ void coordKernel
(dArray<cudaParticle> in, dArray<cuFP_t> mass, dArray<cuFP_t> Afac,
 dArray<cuFP_t> phi, dArray<cuFP_t> Plm, dArray<int> Indx, 
 unsigned int Lmax, unsigned int stride, PII lohi, cuFP_t rmax)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;

  for (int n=0; n<stride; n++) {
    int i = tid*stride + n;
    int npart = i + lohi.first;

    if (npart < lohi.second) {

#ifdef BOUNDS_CHECK
      if (npart>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle p = in._v[npart];
    
      cuFP_t xx = p.pos[0] - sphCen[0];
      cuFP_t yy = p.pos[1] - sphCen[1];
      cuFP_t zz = p.pos[2] - sphCen[2];
      
      cuFP_t r2 = (xx*xx + yy*yy + zz*zz);
      cuFP_t r = sqrt(r2) + FSMALL;
      
#ifdef BOUNDS_CHECK
      if (i>=mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      mass._v[i] = -1.0;
      
      if (r<rmax) {
	
	mass._v[i] = p.mass;
	
	cuFP_t costh = zz/r;
	phi._v[i] = atan2(yy,xx);
	
#ifdef BOUNDS_CHECK
	if (i>=phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	cuFP_t *plm = &Plm._v[psiz*i];
	legendre_v(Lmax, costh, plm);

#ifdef BOUNDS_CHECK
	if (psiz*(i+1)>Plm._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	cuFP_t x  = cu_r_to_xi(r);
	cuFP_t xi = (x - sphXmin)/sphDxi;
	int indx = floor(xi);
	
	if (indx<0) indx = 0;
	if (indx>sphNumr-2) indx = sphNumr - 2;
	  
	Afac._v[i] = cuFP_t(indx+1) - xi;
#ifdef OFF_GRID_ALERT
	if (Afac._v[i]<0.0 or Afac._v[i]>1.0) printf("off grid: x=%f\n", xi);
#endif
	Indx._v[i] = indx;

#ifdef BOUNDS_CHECK
	if (i>=Afac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=Indx._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      }
    }
  }
}


__global__ void coefKernel
(dArray<cuFP_t> coef, dArray<cuFP_t> used, dArray<cudaTextureObject_t> tex,
 dArray<cuFP_t> Mass, dArray<cuFP_t> Afac, dArray<cuFP_t> Phi,
 dArray<cuFP_t> Plm, dArray<int> Indx,  int stride, 
 int l, int m, unsigned Lmax, unsigned int nmax, PII lohi)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;
  const unsigned int N = lohi.second - lohi.first;

  cuFP_t fac0 = 4.0*M_PI;

  for (int istr=0; istr<stride; istr++) {

    int i = tid*stride + istr;

    if (i<N) {

      cuFP_t mass = Mass._v[i];

#ifdef BOUNDS_CHECK
      if (i>=Mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif

      if (mass>0.0) {
				// For accumulating mass of used particles
	if (l==0 and m==0) used._v[i] = mass;

#ifdef BOUNDS_CHECK
	if (i>=Phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	cuFP_t phi  = Phi._v[i];
	cuFP_t cosp = cos(phi*m);
	cuFP_t sinp = sin(phi*m);
	
#ifdef BOUNDS_CHECK
	if (psiz*(i+1)>Plm._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif	
	cuFP_t *plm = &Plm._v[psiz*i];
	
	// Do the interpolation
	//
	cuFP_t a = Afac._v[i];
	cuFP_t b = 1.0 - a;
	int  ind = Indx._v[i];
	
#ifdef BOUNDS_CHECK
	if (i>=Afac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=Indx._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	for (int n=0; n<nmax; n++) {

	  cuFP_t p0 =
#if cuREAL == 4
	    a*tex1D<float>(tex._v[0], ind  ) +
	    b*tex1D<float>(tex._v[0], ind+1) ;
#else
	    a*int2_as_double(tex1D<int2>(tex._v[0], ind  )) +
	    b*int2_as_double(tex1D<int2>(tex._v[0], ind+1)) ;
#endif
	  int k = 1 + l*nmax + n;

#ifdef BOUNDS_CHECK
	  if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	  cuFP_t v = (
#if cuREAL == 4
		     a*tex1D<float>(tex._v[k], ind  ) +
		     b*tex1D<float>(tex._v[k], ind+1)
#else
		     a*int2_as_double(tex1D<int2>(tex._v[k], ind  )) +
		     b*int2_as_double(tex1D<int2>(tex._v[k], ind+1))
#endif
		      ) * p0 * plm[Ilm(l, m)] * Mass._v[i] * fac0;
	  
	  
	  coef._v[(2*n+0)*N + i] = v * cosp;
	  coef._v[(2*n+1)*N + i] = v * sinp;

#ifdef BOUNDS_CHECK
	  if ((2*n+0)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	  if ((2*n+1)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	}
      } else {
	// No contribution from off-grid particles
	for (int n=0; n<nmax; n++) {
	  coef._v[(2*n+0)*N + i] = 0.0;
	  coef._v[(2*n+1)*N + i] = 0.0;
	}
      }
    }
  }
}

__global__ void
forceKernel(dArray<cudaParticle> in, dArray<cuFP_t> coef,
	    dArray<cudaTextureObject_t> tex, dArray<cuFP_t> L1, dArray<cuFP_t> L2,
	    int stride, unsigned Lmax, unsigned int nmax, PII lohi, cuFP_t rmax,
	    bool external)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {
      
#ifdef BOUNDS_CHECK
      if (npart>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle p = in._v[npart];
      
      cuFP_t xx = p.pos[0] - sphCen[0];
      cuFP_t yy = p.pos[1] - sphCen[1];
      cuFP_t zz = p.pos[2] - sphCen[2];
      
      cuFP_t r2 = (xx*xx + yy*yy + zz*zz);
      cuFP_t r  = sqrt(r2) + FSMALL;
      
      cuFP_t costh = zz/r;
      cuFP_t phi   = atan2(yy, xx);
      cuFP_t RR    = xx*xx + yy*yy;
      cuFP_t cosp  = cos(phi);
      cuFP_t sinp  = sin(phi);
      
      cuFP_t *plm1 = &L1._v[psiz*tid];
      cuFP_t *plm2 = &L2._v[psiz*tid];
      
      legendre_v2(Lmax, costh, plm1, plm2);

      int ioff = 0;
      cuFP_t rs = r/sphScale;
      cuFP_t r0 = 0.0;

      if (r>rmax) {
	ioff = 1;
	r0   = r;
	r    = rmax;
	rs   = r/sphScale;
      }

      cuFP_t  x = cu_r_to_xi(rs);
      cuFP_t xi = (x - sphXmin)/sphDxi;
      cuFP_t dx = cu_d_xi_to_r(x)/sphDxi;

      int  ind = floor(xi);
      int  in0 = ind;

      if (in0<0) in0 = 0;
      if (in0>sphNumr-2) in0 = sphNumr - 2;

      if (ind<1) ind = 1;
      if (ind>sphNumr-2) ind = sphNumr - 2;
      
      cuFP_t a0 = (cuFP_t)(in0+1) - xi;
      cuFP_t a1 = (cuFP_t)(ind+1) - xi;
      cuFP_t b0 = 1.0 - a0;
      cuFP_t b1 = 1.0 - a1;

#ifdef OFF_GRID_ALERT
      if (a0<0.0 or a0>1.0)
	printf("forceKernel: off grid [0]: x=%f\n", xi);
      if (a1<0.0 or a1>1.0)
	printf("forceKernel: off grid [1]: x=%f\n", xi);
#endif
      
      // Do the interpolation for the prefactor potential
      //
#if cuREAL == 4
      cuFP_t pm1 = tex1D<float>(tex._v[0], ind-1);
      cuFP_t p00 = tex1D<float>(tex._v[0], ind  );
      cuFP_t pp1 = tex1D<float>(tex._v[0], ind+1);
      cuFP_t p0  = p00;
      cuFP_t p1  = pp1;
      if (in0==0) {
	p1 = p0;
	p0 = tex1D<float>(tex._v[0], 0);
      }
#else
      cuFP_t pm1 = int2_as_double(tex1D<int2>(tex._v[0], ind-1));
      cuFP_t p00 = int2_as_double(tex1D<int2>(tex._v[0], ind  ));
      cuFP_t pp1 = int2_as_double(tex1D<int2>(tex._v[0], ind+1));
      cuFP_t p0  = p00;
      cuFP_t p1  = pp1;
      if (in0==0) {
	p1 = p0;
	p0 = int2_as_double(tex1D<int2>(tex._v[0], 0));
      }
#endif

      // For force accumulation
      //
      cuFP_t potl = 0.0;
      cuFP_t potr = 0.0;
      cuFP_t pott = 0.0;
      cuFP_t potp = 0.0;

      // l loop
      //
      for (int l=0; l<=Lmax; l++) {

	cuFP_t fac1 = (2.0*l + 1.0)/(4.0*M_PI);

	cuFP_t ccos = 1.0;	// For recursion
	cuFP_t ssin = 0.0;

	// m loop
	//
	for (int m=0; m<=l; m++) {
	  
	  int pindx = Ilm(l, m);

	  cuFP_t Plm1 = plm1[pindx];
	  cuFP_t Plm2 = plm2[pindx];
      
#ifdef NAN_CHECK
	  if (std::isnan(Plm1)) {
	    printf("Force isnan for Plm(%d, %d) ioff=%d, costh=%f, z=%f, r=%f\n", l, m, ioff, costh, zz, r);
	  }

	  if (std::isnan(Plm2)) {
	    printf("Force isnan for Plm2(%d, %d) ioff=%d costh=%f, z=%f, r=%f\n", l, m, ioff, costh, zz, r);
	  }
#endif	  

	  cuFP_t pp_c = 0.0;
	  cuFP_t dp_c = 0.0;
	  cuFP_t pp_s = 0.0;
	  cuFP_t dp_s = 0.0;
	  
	  int indxC = Ilmn(l, m, 'c', 0, nmax);
	  int indxS = Ilmn(l, m, 's', 0, nmax);

	  for (size_t n=0; n<nmax; n++) {
	
	    int k = 1 + l*nmax + n;
	
#ifdef BOUNDS_CHECK
	    if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif	
#if cuREAL == 4
	    cuFP_t um1 = tex1D<float>(tex._v[k], ind-1);
	    cuFP_t u00 = tex1D<float>(tex._v[k], ind  );
	    cuFP_t up1 = tex1D<float>(tex._v[k], ind+1);
	    cuFP_t u0  = u00;
	    cuFP_t u1  = up1;
	    if (in0==0) {
	      u1 = u0;
	      u0 = tex1D<float>(tex._v[k], 0);
	    }
#else
	    cuFP_t um1 = int2_as_double(tex1D<int2>(tex._v[k], ind-1));
	    cuFP_t u00 = int2_as_double(tex1D<int2>(tex._v[k], ind  ));
	    cuFP_t up1 = int2_as_double(tex1D<int2>(tex._v[k], ind+1));
	    cuFP_t u0  = u00;
	    cuFP_t u1  = up1;
	    if (in0==0) {
	      u1 = u0;
	      u0 = int2_as_double(tex1D<int2>(tex._v[k], 0));
	    }
#endif
	    cuFP_t v = (a0*u0 + b0*u1)*(a0*p0 + b0*p1);
	    
	    cuFP_t dv =
	      dx * ( (b1 - 0.5)*um1*pm1 - 2.0*b1*u00*p00 + (b1 + 0.5)*up1*pp1 );
	    
#ifdef NAN_CHECK
	    if (std::isnan(v))
	      printf("v tab nan: (%d, %d): a=%f b=%f p0=%f p1=%f u0=%f u1=%f\n", l, m, a1, b1, p00, pp1, u00, up1);

	    if (std::isnan(dv))
	      printf("dv tab nan: (%d, %d): a=%f b=%f pn=%f p0=%f pp=%f un=%f u0=%f up=%f\n", l, m, a1, b1, pm1, p00, pp1, um1, u00, up1);
#endif

	    pp_c -=  v * coef._v[indxC+n];
	    dp_c -= dv * coef._v[indxC+n];
	    if (m>0) {
	      pp_s -=  v * coef._v[indxS+n];
	      dp_s -= dv * coef._v[indxS+n];
	    }

	  } // END: n loop
	  
#ifdef NAN_CHECK
	  if (std::isnan(pp_c)) printf("pp_c eval nan: (%d, %d): r=%f r0=%f\n", l, m, r, r0);
	  if (std::isnan(dp_c)) printf("dp_c eval nan: (%d, %d): r=%f r0=%f\n", l, m, r, r0);
	  if (std::isnan(pp_s)) printf("pp_s eval nan: (%d, %d): r=%f r0=%f\n", l, m, r, r0);
	  if (std::isnan(dp_s)) printf("dp_s eval nan: (%d, %d): r=%f r0=%f\n", l, m, r, r0);
#endif

	  if (m==0) {

	    if (ioff) {
	      pp_c *= pow(rmax/r0, (cuFP_t)(l+1));
	      dp_c  = -pp_c/r0 * (cuFP_t)(l+1);
#ifdef NAN_CHECK
	      if (std::isnan(pp_c)) printf("Force nan [ioff]: l=%d, r=%f r0=%f\n", l, r, r0);
#endif
	    }
	    
	    potl += fac1 * pp_c * Plm1;
	    potr += fac1 * dp_c * Plm1;
	    pott += fac1 * pp_c * Plm2;
	    potp += 0.0;
	    
	  } else {

	    if (ioff) {
	      cuFP_t facp  = pow(rmax/r0,(cuFP_t)(l+1));
	      cuFP_t facdp = -facp/r0 * (l+1);
	      pp_c *= facp;
	      pp_s *= facp;
	      dp_c = pp_c * facdp;
	      dp_s = pp_s * facdp;

#ifdef NAN_CHECK
	      if (std::isnan(pp_c)) {
		printf("Force nan: c l=%d, m=%d, r=%f\n", l, m, r);
	      }
	      if (std::isnan(dp_s)) {
		printf("Force nan: s l=%d, m=%d, r=%f\n", l, m, r);
	      }
	      if (std::isnan(dp_c)) {
		printf("Force nan: dc l=%d, m=%d, r=%f\n", l, m, r);
	      }
	      if (std::isnan(pp_s)) {
		printf("Force nan: ds l=%d, m=%d, r=%f\n", l, m, r);
	      }
#endif
	    }

	    // Factorials
	    //
	    cuFP_t numf = 1.0, denf = 1.0;
	    for (int i=2; i<=l+m; i++) {
	      if (i<=l-m) numf *= i; denf *= i;
	    }
	    
	    cuFP_t fac2 = 2.0 * numf/denf * fac1;
	    
	    potl += fac2 * Plm1 * ( pp_c*ccos + pp_s*ssin);
	    potr += fac2 * Plm1 * ( dp_c*ccos + dp_s*ssin);
	    pott += fac2 * Plm2 * ( pp_c*ccos + pp_s*ssin);
	    potp += fac2 * Plm1 * (-pp_c*ssin + pp_s*ccos)*m;
	  }

	  // Trig recursion to squeeze avoid internal FP fct call
	  //
	  cuFP_t cosM = ccos;
	  cuFP_t sinM = ssin;

	  ccos = cosM * cosp - sinM * sinp;
	  ssin = sinM * cosp + cosM * sinp;

	} // END: m loop

      } // END: l loop
				// Rescale
      potr /= sphScale*sphScale;
      potl /= sphScale;
      pott /= sphScale;
      potp /= sphScale;

      in._v[npart].acc[0] += -(potr*xx/r - pott*xx*zz/(r*r*r));
      in._v[npart].acc[1] += -(potr*yy/r - pott*yy*zz/(r*r*r));
      in._v[npart].acc[2] += -(potr*zz/r + pott*RR/(r*r*r));
      if (RR > FSMALL) {
	in._v[npart].acc[0] +=  potp*yy/RR;
	in._v[npart].acc[1] += -potp*xx/RR;
      }
      if (external)
	in._v[npart].potext += potl;
      else
	in._v[npart].pot    += potl;

#ifdef NAN_CHECK
      // Sanity check
      bool bad = false;
      for (int k=0; k<3; k++) {
	if (std::isnan(in._v[npart].acc[k])) bad = true;
      }

      if (bad)  {
	printf("Force nan value: [%d] x=%f X=%f Y=%f Z=%f r=%f R=%f P=%f dP/dr=%f dP/dx=%f dP/dp=%f\n", in._v[npart].indx, x, xx, yy, zz, r, RR, potl, potr, pott, potp);
	if (a<0.0 or a>1.0)  {
	  printf("Force nan value, no ioff: [%d] x=%f xi=%f dxi=%f a=%f i=%d\n",
		 in._v[npart].indx, x, xi, dx, a, ind);
	}
      }
#endif

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

static bool initialize_cuda_sph = true;
static unsigned dbg_id = 0;

void SphericalBasis::cudaStorage::resize_coefs
(int nmax, int Lmax, int N, int gridSize, int sampT, bool pca)
{
  // Create space for coefficient reduction to prevent continued
  // dynamic allocation
  //
  if (dN_coef.capacity() < 2*nmax*N)
    dN_coef.reserve(2*nmax*N);
  
  if (dc_coef.capacity() < 2*nmax*gridSize)
    dc_coef.reserve(2*nmax*gridSize);
  
  if (plm1_d.capacity() < (Lmax+1)*(Lmax+2)/2*N)
    plm1_d.reserve((Lmax+1)*(Lmax+2)/2*N);
  
  if (r_d.capacity() < N) r_d.reserve(N);
  if (m_d.capacity() < N) m_d.reserve(N);
  if (u_d.capacity() < N) u_d.reserve(N);
  if (a_d.capacity() < N) a_d.reserve(N);
  if (p_d.capacity() < N) p_d.reserve(N);
  if (i_d.capacity() < N) i_d.reserve(N);
  
  // Fixed size arrays
  //
  if (pca) {
    T_coef.resize(sampT);
    for (int T=0; T<sampT; T++) {
      T_coef[T].resize((Lmax+1)*(Lmax+2)*nmax);
    }
  }

  // Set needed space for current step
  //
  dN_coef.resize(2*nmax*N);
  dc_coef.resize(2*nmax*gridSize);

  // This will stay fixed for the entire run
  //
  dw_coef.resize(2*nmax);

  // Space for Legendre coefficients 
  //
  if (plm1_d.capacity() < (Lmax+1)*(Lmax+2)/2*N)
    plm1_d.reserve((Lmax+1)*(Lmax+2)/2*N);
  plm1_d.resize((Lmax+1)*(Lmax+2)/2*N);
  
  // Space for coordinates
  //
  r_d.resize(N);
  m_d.resize(N);
  u_d.resize(N);
  a_d.resize(N);
  p_d.resize(N);
  i_d.resize(N);
}

void SphericalBasis::zero_coefs()
{
  Component::     cuRingType cr = *cC->cuRing.get();
  SphericalBasis::cuRingType ar = *cuRing.get();
  
  for (int n=0; n<cuStreams; n++) {
    // Resize output array
    //
    ar->df_coef.resize((Lmax+1)*(Lmax+2)*nmax);
    
    // Zero output array
    //
    thrust::fill(thrust::cuda::par.on(cr->stream),
		 ar->df_coef.begin(), ar->df_coef.end(), 0.0);

    // Resize and zero PCA arrays
    //
    if (pca) {
      if (ar->T_coef.size() != sampT) {
	ar->T_coef.resize(sampT);
	for (int T=0; T<sampT; T++) {
	  ar->T_coef[T].resize((Lmax+1)*(Lmax+2)*nmax);
	}
	host_massT.resize(sampT);
      }

      for (int T=0; T<sampT; T++)
	thrust::fill(thrust::cuda::par.on(cr->stream),
		     ar->T_coef[T].begin(), ar->T_coef[T].end(), 0.0);

      thrust::fill(host_massT.begin(), host_massT.end(), 0.0);
    }

    // Advance iterators
    //
    ++cr;			// Component stream
    ++ar;			// Force method storage
  }
}

void SphericalBasis::cudaStorage::resize_acc(int Lmax, int Nthread)
{
  // Space for Legendre coefficients 
  //
  if (plm1_d.capacity() < (Lmax+1)*(Lmax+2)/2*Nthread)
    plm1_d.reserve((Lmax+1)*(Lmax+2)/2*Nthread);
  
  if (plm2_d.capacity() < (Lmax+1)*(Lmax+2)/2*Nthread)
    plm2_d.reserve((Lmax+1)*(Lmax+2)/2*Nthread);
  
  plm1_d.resize((Lmax+1)*(Lmax+2)/2*Nthread);
  plm2_d.resize((Lmax+1)*(Lmax+2)/2*Nthread);
}


void SphericalBasis::determine_coefficients_cuda(bool compute)
{
  if (pca) {

    if (sampT == 0) {		// Allocate storage
      sampT = floor(sqrt(cC->nbodies_tot));
      massT    .resize(sampT, 0);
      massT1   .resize(sampT, 0);

      expcoefT .resize(sampT);
      for (auto & t : expcoefT ) t = MatrixP(new Matrix(0, Lmax*(Lmax+2), 1, nmax));
	
      expcoefT1.resize(sampT);
      for (auto & t : expcoefT1) t = MatrixP(new Matrix(0, Lmax*(Lmax+2), 1, nmax));
    }
  }

  if (compute && mlevel==0) {	// Zero arrays
    for (auto & t : expcoefT1) t->zero();
    for (auto & v : massT1)    v = 0;
  }

  if (initialize_cuda_sph) {
    initialize_cuda();
    initialize_mapping_constants();
    initialize_cuda_sph = false;

    // Copy texture objects to device
    //
    t_d = tex;
  }

  std::cout << std::scientific;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, component->cudaDevice);

  Component::     cuRingType cr = *cC->cuRing.get();
  SphericalBasis::cuRingType ar = *cuRing.get();
  
  // This will stay fixed for the entire run
  //
  host_coefs.resize((Lmax+1)*(Lmax+1)*nmax);

  if (pca) {
    host_massT.resize(sampT);
  }

  // Center assignment to symbol data
  //
  std::vector<cuFP_t> ctr;
  for (auto v : cC->getCenter(Component::Local | Component::Centered))
    ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(sphCen, &ctr[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphCen");
  
  // For debugging
  //
  static bool firstime = false;
  
  if (firstime) {
    testConstants<<<1, 1, 0, cr->stream>>>();
    cudaDeviceSynchronize();
    firstime = false;
  }
  
  // Zero counter and coefficients
  //
  unsigned Ntot = 0;
  use[0] = 0.0;
  thrust::fill(host_coefs.begin(), host_coefs.end(), 0.0);

  if (compute) thrust::fill(host_massT.begin(), host_massT.end(), 0.0);

  // Zero out coefficient storage
  //
  zero_coefs();

  // Copy particles to host vector
  //
  cC->ParticlesToCuda();

  // Loop over bunches
  //
  size_t psize  = cC->host_particles.size();

  Component::hostPartItr begin = cC->host_particles.begin();
  Component::hostPartItr first = begin;
  Component::hostPartItr last  = begin;
  Component::hostPartItr end   = cC->host_particles.end();

  if (psize <= cC->bunchSize) last = end;
  else std::advance(last, cC->bunchSize);

  std::vector<cudaStream_t> f_s;
  thrust::device_vector<unsigned int> f_use;

  while (std::distance(first, last)) {
    
    cr->first = first;
    cr->last  = last;
    cr->id    = ++dbg_id;

    // Copy bunch to device
    //
    cC->HostToDev(cr);

    // Sort particles and get coefficient size
    //
    PII lohi = cC->CudaSortByLevel(cr, mlevel, mlevel);
    
    // Compute grid
    //
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;

    if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

#ifdef VERBOSE
    static debug_max_count = 10;
    static debug_cur_count = 0;
    if (debug_cur_count++ < debug_max_count) {
      std::cout << std::endl << "**" << std::endl
		<< "** N      = " << N          << std::endl
		<< "** Stride = " << stride     << std::endl
		<< "** Block  = " << BLOCK_SIZE << std::endl
		<< "** Grid   = " << gridSize   << std::endl
		<< "** Xcen   = " << ctr[0]     << std::endl
		<< "** Ycen   = " << ctr[1]     << std::endl
		<< "** Zcen   = " << ctr[2]     << std::endl
		<< "**" << std::endl;
    }
#endif

    if (N) {
    
      // Resize storage as needed
      //
      ar->resize_coefs(nmax, Lmax, N, gridSize, sampT, pca);
      
      // Shared memory size for the reduction
      //
      int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
      
    
      thrust::counting_iterator<int> index_begin(0);
      thrust::counting_iterator<int> index_end(gridSize*2*nmax);

      // Do the work
      //
				// Compute the coordinate
				// transformation
				// 
      coordKernel<<<gridSize, BLOCK_SIZE, 0, cr->stream>>>
	(toKernel(cr->cuda_particles),
	 toKernel(ar->m_d), toKernel(ar->a_d), toKernel(ar->p_d),
	 toKernel(ar->plm1_d), toKernel(ar->i_d),
	 Lmax, stride, lohi, rmax);
    
				// Compute the coefficient
				// contribution for each order
      int osize = nmax*2;	//
      auto beg = ar->df_coef.begin();
      std::vector<thrust::device_vector<cuFP_t>::iterator> bg;
      for (int T=0; T<sampT; T++) bg.push_back(ar->T_coef[T].begin());

      thrust::fill(ar->u_d.begin(), ar->u_d.end(), 0.0);

      for (int l=0; l<=Lmax; l++) {
	for (int m=0; m<=l; m++) {
				// Compute the contribution to the
				// coefficients from each particle
				//
	  coefKernel<<<gridSize, BLOCK_SIZE, 0, cr->stream>>>
	    (toKernel(ar->dN_coef), toKernel(ar->u_d),
	     toKernel(t_d), toKernel(ar->m_d),
	     toKernel(ar->a_d), toKernel(ar->p_d), toKernel(ar->plm1_d),
	     toKernel(ar->i_d), stride, l, m, Lmax, nmax, lohi);
	  
				// Begin the reduction per grid block
				// [perhaps this should use a stride?]
	  unsigned int gridSize1 = N/BLOCK_SIZE;
	  if (N > gridSize1*BLOCK_SIZE) gridSize1++;
	  reduceSum<cuFP_t, BLOCK_SIZE><<<gridSize1, BLOCK_SIZE, sMemSize, cr->stream>>>
	    (toKernel(ar->dc_coef), toKernel(ar->dN_coef), osize, N);
      
				// Finish the reduction for this order
				// in parallel
	  thrust::reduce_by_key
	    (
	     thrust::cuda::par.on(cr->stream),
	     thrust::make_transform_iterator(index_begin, key_functor(gridSize)),
	     thrust::make_transform_iterator(index_end,   key_functor(gridSize)),
	     ar->dc_coef.begin(), thrust::make_discard_iterator(), ar->dw_coef.begin()
	   );

	  thrust::transform(thrust::cuda::par.on(cr->stream),
			    ar->dw_coef.begin(), ar->dw_coef.end(),
			    beg, beg, thrust::plus<cuFP_t>());
	  
	  thrust::advance(beg, osize);

	  if (compute) {

	    int sN = N/sampT;

	    for (int T=0; T<sampT; T++) {
	      int k = sN*T;	// Starting position
	      int s = sN;
	      if (T==sampT-1) s = N - k;
	    
				// Begin the reduction per grid block
				//
	      /*
	      unsigned int stride1   = s/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
	      unsigned int gridSize1 = s/BLOCK_SIZE/stride1;

	      if (s > gridSize1*BLOCK_SIZE*stride1) gridSize1++;
	      */

	      unsigned int gridSize1 = s/BLOCK_SIZE;
	      if (s > gridSize1*BLOCK_SIZE) gridSize1++;
	      reduceSumS<cuFP_t, BLOCK_SIZE><<<gridSize1, BLOCK_SIZE, sMemSize, cr->stream>>>
		(toKernel(ar->dc_coef), toKernel(ar->dN_coef), osize, N, k, k+s);
      
				// Finish the reduction for this order
				// in parallel
	      thrust::reduce_by_key
		(
		 thrust::cuda::par.on(cr->stream),
		 thrust::make_transform_iterator(index_begin, key_functor(gridSize)),
		 thrust::make_transform_iterator(index_end,   key_functor(gridSize)),
		 ar->dc_coef.begin(), thrust::make_discard_iterator(), ar->dw_coef.begin()
		 );


	      thrust::transform(thrust::cuda::par.on(cr->stream),
				ar->dw_coef.begin(), ar->dw_coef.end(),
				bg[T], bg[T], thrust::plus<cuFP_t>());

	      thrust::advance(bg[T], osize);

	      if (l==0 and m==0) {
		auto mbeg = ar->u_d.begin();
		auto mend = mbeg;
		thrust::advance(mbeg, sN*T);
		if (T<sampT-1) thrust::advance(mend, sN*(T+1));
		else mend = ar->u_d.end();
		
		// TEST
		thrust::host_vector<cuFP_t> tst = ar->u_d;
		// END TEST

		host_massT[T] += thrust::reduce(mbeg, mend);
	      }
	    }
	  }
	}
      }

      // Compute number and total mass of particles used in coefficient
      // determination
      //
      thrust::sort(thrust::cuda::par.on(cr->stream), ar->m_d.begin(), ar->m_d.end());
      
      // Asynchronously cache result for host side to prevent stream block
      //
      cudaStream_t s1;		// Create a new non blocking stream
      cudaStreamCreateWithFlags(&s1, cudaStreamNonBlocking);
      f_s.push_back(s1);
				
      size_t fsz = f_s.size();	// Augment the data vector
      f_use.resize(fsz);
				// Call the kernel on a single thread
				// 
      reduceUseSph<<<1, 1, 0, s1>>>(ar->m_d.begin(), ar->m_d.end(),
				    &f_use[fsz-1]);

      Ntot += N;
    }

    // Advance iterators
    //
    first = last;
    size_t nadv = std::distance(first, end);
    if (nadv <= cC->bunchSize) last = end;
    else std::advance(last, cC->bunchSize);

    // Advance stream iterators (these are designed to increment in
    // lock step)
    //
    ++cr;			// Component stream
    ++ar;			// Force method storage
  }


  // Copy back coefficient data from device and load the host
  //
  for (auto r : cuRingData) {
    thrust::host_vector<cuFP_t> ret = r.df_coef;
    int offst = 0;
    for (int l=0; l<=Lmax; l++) {
      for (int m=0; m<=l; m++) {
	for (size_t j=0; j<nmax; j++) {
	  host_coefs[Ilmn(l, m, 'c', j, nmax)] += ret[2*j+offst];
	  if (m>0) host_coefs[Ilmn(l, m, 's', j, nmax)] += ret[2*j+1+offst];
	}
	offst += nmax*2;
      }
    }
  }
  
  // Get the on-grid count from the threads
  //
  for (auto & s : f_s) {	// Synchronize and dallocate streams
    cudaStreamSynchronize(s);
    cudaStreamDestroy(s);
  }
				// Copy data back from device
  thrust::host_vector<unsigned int> f_ret(f_use);
  for (auto v : f_ret) use[0] += v;

  if (Ntot == 0) {
    return;
  }

  // DEBUG, only useful for CUDAtest branch
  //
  if (false) {

    constexpr bool compareC = false;


    if (compareC) {
      std::cout << std::string(2*4+4*20, '-') << std::endl
		<< "---- Spherical "      << std::endl
		<< std::string(2*4+4*20, '-') << std::endl;
      std::cout << "L=M=0 coefficients" << std::endl
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
		<< "---- Spherical "      << std::endl
		<< std::string(2*4+20, '-') << std::endl;
      std::cout << "L=M=0 coefficients" << std::endl
		<< std::setprecision(10);

      std::cout << std::setw(4)  << "n"
		<< std::setw(4)  << "i"
		<< std::setw(20) << "GPU"
		<< std::endl;
    }

    int i = Ilmn(0, 0, 'c', 0, nmax);
    auto cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      cuFP_t a = host_coefs[i+n];
      cuFP_t b = expcoef[0][n+1];
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

    std::cout << "L=1, M=0 coefficients" << std::endl;

    i = Ilmn(1, 0, 'c', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      cuFP_t a = host_coefs[i+n];
      cuFP_t b = expcoef[1][n+1];
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

    std::cout << "L=1, M=1c coefficients" << std::endl;

    i = Ilmn(1, 1, 'c', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      cuFP_t a = host_coefs[i+n];
      cuFP_t b = expcoef[2][n+1];
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

    std::cout << "L=1, M=1s coefficients" << std::endl;

    i = Ilmn(1, 1, 's', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      cuFP_t a = host_coefs[i+n];
      cuFP_t b = expcoef[3][n+1];
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

    std::cout << "L=2, M=0 coefficients" << std::endl;

    i = Ilmn(2, 0, 'c', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      cuFP_t a = host_coefs[i+n];
      cuFP_t b = expcoef[4][n+1];
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
      
    std::cout << "L=2, M=1c coefficients" << std::endl;

    i = Ilmn(2, 1, 'c', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      cuFP_t a = host_coefs[i+n];
      cuFP_t b = expcoef[5][n+1];
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

    std::cout << "L=2, M=1c coefficients" << std::endl;

    i = Ilmn(2, 1, 's', 0, nmax);
    cmax = std::max_element(host_coefs.begin()+i, host_coefs.begin()+i+nmax, LessAbs<cuFP_t>());

    for (size_t n=0; n<nmax; n++) {
      cuFP_t a = host_coefs[i+n];
      cuFP_t b = expcoef[6][n+1];
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
      
      int  l;
      int  m;
      int  n;
      
      char cs;
    }
    elem;

    std::multimap<double, Element> compare;

    std::ofstream out("test_sph.dat");

    //		l loop
    for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
      //		m loop
      for (int m=0, moffset=0; m<=l; m++) {
	
	if (m==0) {
	  for (int n=1; n<=nmax; n++) {
	    elem.l = l;
	    elem.m = m;
	    elem.n = n;
	    elem.cs = 'c';
	    elem.d = expcoef[loffset+moffset][n];
	    elem.f = host_coefs[Ilmn(l, m, 'c', n-1, nmax)];
	    
	    double test = fabs(elem.d - elem.f);
	    if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);
	    
	    compare.insert(std::make_pair(test, elem));
	    
	    out << std::setw( 5) << l
		<< std::setw( 5) << m
		<< std::setw( 5) << n
		<< std::setw( 5) << 'c'
		<< std::setw( 5) << Ilmn(l, m, 'c', n-1, nmax)
		<< std::setw(14) << elem.d
		<< std::setw(14) << elem.f
		<< std::endl;
	  }
	  
	  moffset++;
	}
	else {
	  for (int n=1; n<=nmax; n++) {
	    elem.l = l;
	    elem.m = m;
	    elem.n = n;
	    elem.cs = 'c';
	    elem.d = expcoef[loffset+moffset][n];
	    elem.f = host_coefs[Ilmn(l, m, 'c', n-1, nmax)];

	    out << std::setw( 5) << l
		<< std::setw( 5) << m
		<< std::setw( 5) << n
		<< std::setw( 5) << 'c'
		<< std::setw( 5) << Ilmn(l, m, 'c', n-1, nmax)
		<< std::setw(14) << elem.d
		<< std::setw(14) << elem.f
		<< std::endl;

	    double test = fabs(elem.d - elem.f);
	    if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);

	    compare.insert(std::make_pair(test, elem));
	  }
	  for (int n=1; n<=nmax; n++) {
	    elem.l = l;
	    elem.m = m;
	    elem.n = n;
	    elem.cs = 's';
	    elem.d = expcoef[loffset+moffset+1][n];
	    elem.f = host_coefs[Ilmn(l, m, 's', n-1, nmax)];

	    out << std::setw( 5) << l
		<< std::setw( 5) << m
		<< std::setw( 5) << n
		<< std::setw( 5) << 's'
		<< std::setw( 5) << Ilmn(l, m, 's', n-1, nmax)
		<< std::setw(14) << elem.d
		<< std::setw(14) << elem.f
		<< std::endl;
	    
	    double test = fabs(elem.d - elem.f);
	    if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);
	    
	    compare.insert(std::make_pair(test, elem));
	  }
	  moffset+=2;
	}
      }
    }
    
    std::map<double, Element>::iterator best = compare.begin();
    std::map<double, Element>::iterator midl = best;
    std::advance(midl, compare.size()/2);
    std::map<double, Element>::reverse_iterator last = compare.rbegin();
    
    std::cout << std::string(4*2 + 3*20 + 22, '-') << std::endl
	      << "---- Spherical coefficients" << std::endl
	      << std::string(4*2 + 3*20 + 22, '-') << std::endl;

    std::cout << "Best case: ["
	      << std::setw( 2) << best->second.l << ", "
	      << std::setw( 2) << best->second.m << ", "
	      << std::setw( 2) << best->second.n << ", "
	      << std::setw( 2) << best->second.cs << "] = "
	      << std::setw(20) << best->second.d
	      << std::setw(20) << best->second.f
	      << std::setw(20) << fabs(best->second.d - best->second.f)
	      << std::endl;
  
    std::cout << "Mid case:  ["
	      << std::setw( 2) << midl->second.l << ", "
	      << std::setw( 2) << midl->second.m << ", "
	      << std::setw( 2) << midl->second.n << ", "
	      << std::setw( 2) << midl->second.cs << "] = "
	      << std::setw(20) << midl->second.d
	      << std::setw(20) << midl->second.f
	      << std::setw(20) << fabs(midl->second.d - midl->second.f)
	      << std::endl;
    
    std::cout << "Last case: ["
	      << std::setw( 2) << last->second.l << ", "
	      << std::setw( 2) << last->second.m << ", "
	      << std::setw( 2) << last->second.n << ", "
	      << std::setw( 2) << last->second.cs << "] = "
	      << std::setw(20) << last->second.d
	      << std::setw(20) << last->second.f
	      << std::setw(20) << fabs(last->second.d - last->second.f)
	      << std::endl;
  }
}

void SphericalBasis::determine_acceleration_cuda()
{
  if (initialize_cuda_sph) {
    initialize_cuda();
    initialize_mapping_constants();
    initialize_cuda_sph = false;

    // Copy texture objects to device
    //
    t_d = tex;
  }

  std::cout << std::scientific;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, component->cudaDevice);

  // Stream structure iterators
  //
  Component::     cuRingType cr = *cC->cuRing.get();
  SphericalBasis::cuRingType ar = *cuRing.get();

  // Assign component center
  //
  std::vector<cuFP_t> ctr;
  for (auto v : cC->getCenter(Component::Local | Component::Centered)) ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(sphCen, &ctr[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphCen");
  
  // Copy particles to host vector
  //
  cC->ParticlesToCuda();

  // Loop over bunches
  //
  size_t psize  = cC->host_particles.size();

  Component::hostPartItr begin = cC->host_particles.begin();
  Component::hostPartItr first = begin;
  Component::hostPartItr last  = begin;
  Component::hostPartItr end   = cC->host_particles.end();

  if (psize <= cC->bunchSize) last = end;
  else std::advance(last, cC->bunchSize);

  unsigned Ntot = 0;

  while (std::distance(first, last)) {

    cr->first = first;
    cr->last  = last;
    cr->id    = ++dbg_id;

    // Copy bunch to device
    //
    cC->HostToDev(cr);

    // Sort particles and get size
    //
    PII lohi = cC->CudaSortByLevel(cr, toplev, multistep);

    // Compute grid
    //
    unsigned int N         = lohi.second - lohi.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
    
    if (N>0) {

      Ntot += N;

      if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

      unsigned int Nthread = gridSize*BLOCK_SIZE;

#ifdef VERBOSE
      static debug_max_count = 10;
      static debug_cur_count = 0;
      if (debug_cur_count++ < debug_max_count) {
	std::cout << std::endl << "**" << std::endl
		  << "** N      = " << N          << std::endl
		  << "** Stride = " << stride     << std::endl
		  << "** Block  = " << BLOCK_SIZE << std::endl
		  << "** Grid   = " << gridSize   << std::endl
		  << "** Xcen   = " << ctr[0]     << std::endl
		  << "** Ycen   = " << ctr[1]     << std::endl
		  << "** Zcen   = " << ctr[2]     << std::endl
		  << "**" << std::endl;
      }
#endif
      
      ar->resize_acc(Lmax, Nthread);

      // Shared memory size for the reduction
      //
      int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
      
      // Do the work
      //
      forceKernel<<<gridSize, BLOCK_SIZE, sMemSize, cr->stream>>>
	(toKernel(cr->cuda_particles), toKernel(dev_coefs), toKernel(t_d),
	 toKernel(ar->plm1_d), toKernel(ar->plm2_d),
	 stride, Lmax, nmax, lohi, rmax, use_external);

      // Copy particles back to host
      //
      cC->DevToHost(cr);
    }

    // Advance iterators
    //
    first = last;
    size_t nadv = std::distance(first, end);
    if (nadv <= cC->bunchSize) last = end;
    else std::advance(last, cC->bunchSize);
    
    // Advance stream iterators
    //
    ++cr;			// Component
    ++ar;			// Force method
  }

  // Finally, do copy from host to component
  //
  cC->CudaToParticles();

  // DEBUGGING TEST
  if (false) {
    std::cout << std::string(10+7*16, '-') << std::endl;
    std::cout << "---- Acceleration in SphericalBasis [T=" << tnow
	      << ", N=" << Ntot << ", level=" << mlevel
	      << ", name=" << cC->name << "]" << std::endl;
    std::cout << std::string(10+7*16, '-') << std::endl;
    first = last = begin;
    std::advance(last, 5);
    std::copy(first, last, std::ostream_iterator<cudaParticle>(std::cout, "\n"));
    first = begin;
    last  = end;
    std::advance(first, psize-5);
    std::copy(first, last, std::ostream_iterator<cudaParticle>(std::cout, "\n"));
    std::cout << std::string(10+7*16, '-') << std::endl;
  }
}

void SphericalBasis::HtoD_coefs(const Matrix& expcoef)
{
  host_coefs.resize((Lmax+1)*(Lmax+1)*nmax); // Should stay fixed, no reserve

  // l loop
  //
  for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
    // m loop
    //
    for (int m=0, moffset=0; m<=l; m++) {
	
      // n loop
      //
      for (int n=1; n<=nmax; n++) {
	host_coefs[Ilmn(l, m, 'c', n-1, nmax)] = expcoef[loffset+moffset][n];
	if (m>0) host_coefs[Ilmn(l, m, 's', n-1, nmax)] = expcoef[loffset+moffset+1][n];
      }

      if (m>0) moffset += 2;
      else     moffset += 1;
    }
  }

  dev_coefs = host_coefs;
}


void SphericalBasis::DtoH_coefs(Matrix& expcoef)
{
  // l loop
  //
  for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {

    // m loop
    //
    for (int m=0, moffset=0; m<=l; m++) {
	
      // n loop
      //
      for (int n=1; n<=nmax; n++) {
	expcoef0[0][loffset+moffset][n] = host_coefs[Ilmn(l, m, 'c', n-1, nmax)];
	if (m>0) expcoef0[0][loffset+moffset+1][n] = host_coefs[Ilmn(l, m, 's', n-1, nmax)];
      }

      if (m>0) moffset += 2;
      else     moffset += 1;
    }
  }

  if (compute) {

    for (int T=0; T<sampT; T++) {
      massT1[T] += host_massT[T];
      muse1 [0] += host_massT[T];
    }

    for (auto r : cuRingData) {
      // T loop
      //
      for (int T=0; T<sampT; T++) {

	thrust::host_vector<cuFP_t> ret = r.T_coef[T];

	int offst = 0;
	int osize = 2.0*nmax;

	// l loop
	//
	for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {

	  // m loop
	  //
	  for (int m=0, moffset=0; m<=l; m++) {
	
	    // n loop
	    //
	    for (int n=1; n<=nmax; n++) {
	      (*expcoefT1[T])[loffset+moffset][n] = ret[2*(n-1) + offst];
	      if (m>0) (*expcoefT1[T])[loffset+moffset+1][n] = ret[2*(n-1) + 1 + offst];
	    }

	    offst += osize;

	    if (m>0) moffset += 2;
	    else     moffset += 1;
	  }
	}
      }
    }
  }
}

void SphericalBasis::destroy_cuda()
{
  for (size_t i=0; i<tex.size(); i++) {
    std::ostringstream sout;
    sout << "trying to free TextureObject [" << i << "]";
    cuda_safe_call(cudaDestroyTextureObject(tex[i]),
		   __FILE__, __LINE__, sout.str());
  }

  for (size_t i=0; i<cuInterpArray.size(); i++) {
    std::ostringstream sout;
    sout << "trying to free cuArray [" << i << "]";
    cuda_safe_call(cudaFreeArray(cuInterpArray[i]),
		     __FILE__, __LINE__, sout.str());
  }
}
