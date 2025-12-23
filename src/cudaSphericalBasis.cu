// -*- C++ -*-

#include <tuple>
#include <list>

#include <Component.H>
#include <SphericalBasis.H>
#include <cudaReduce.cuH>

// Define for debugging
//
// #define OFF_GRID_ALERT
// #define BOUNDS_CHECK
// #define VERBOSE_CTR
// #define NAN_CHECK
// #define VERBOSE_TIMING
// #define VERBOSE_DBG

// Global symbols for coordinate transformation in SphericalBasis
//
__device__ __constant__
cuFP_t sphScale, sphRscale, sphHscale, sphXmin, sphXmax, sphDxi, sphCen[3], sphEPS;

__device__ __constant__
int    sphNumr, sphCmap;

__device__ __constant__
bool   sphAcov, sphNO_L0, sphNO_L1, sphEVEN_L, sphEVEN_M, sphM0only;

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


__device__
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

  if (1.0-fabs(x) < sphEPS) {
    if (x>0) x =   1.0 - sphEPS;
    else     x = -(1.0 - sphEPS);
  }

  somx2 = 1.0/(x*x - 1.0);
  dp[0] = 0.0;
  for (int l=1; l<=lmax; l++) {
    for (int m=0; m<l; m++)
      dp[Ilm(l, m)] = somx2*(x*l*p[Ilm(l, m)] - (l+m)*p[Ilm(l-1, m)]);
    dp[Ilm(l, l)] = somx2*x*l*p[Ilm(l, l)];
  }
}


__global__
void testConstantsSph()
{
  printf("-------------------------\n");
  printf("---Spherical constants---\n");
  printf("-------------------------\n");
  printf("   Scale  = %f\n", sphScale );
  printf("   Rscale = %f\n", sphRscale);
  printf("   Xmin   = %f\n", sphXmin  );
  printf("   Xmax   = %f\n", sphXmax  );
  printf("   Dxi    = %f\n", sphDxi   );
  printf("   Numr   = %d\n", sphNumr  );
  printf("   Cmap   = %d\n", sphCmap  );
  printf("   Feps   = %f\n", sphEPS   );
  printf("-------------------------\n");
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
  // Nothing so far
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

  cuda_safe_call(cudaMemcpyToSymbol(sphCmap,   &f.cmapR,  sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphCmap");

  cuda_safe_call(cudaMemcpyToSymbol(sphAcov,   &subsamp,  sizeof(bool),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphAcov");
  cuda_safe_call(cudaMemcpyToSymbol(sphNO_L0,  &NO_L0,    sizeof(bool),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphNO_L0");
  cuda_safe_call(cudaMemcpyToSymbol(sphNO_L1,  &NO_L1,    sizeof(bool),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphNO_L1");
  cuda_safe_call(cudaMemcpyToSymbol(sphEVEN_L, &EVEN_L,    sizeof(bool),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphEVEN_L");
  cuda_safe_call(cudaMemcpyToSymbol(sphEVEN_M, &EVEN_M,   sizeof(bool),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphEVEN_M");
  cuda_safe_call(cudaMemcpyToSymbol(sphM0only, &M0_only,  sizeof(bool),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphM0only");

#if cuREAL == 4
  cuFP_t feps = 3.0*std::numeric_limits<float>::epsilon();
#else
  cuFP_t feps = 3.0*std::numeric_limits<double>::epsilon();
#endif

  cuda_safe_call(cudaMemcpyToSymbol(sphEPS, &feps,  sizeof(cuFP_t),  size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphEPS");
}


__global__ void coordKernel
(dArray<cudaParticle> P, dArray<int> I, dArray<cuFP_t> mass,
 dArray<cuFP_t> Afac, dArray<cuFP_t> phi, dArray<cuFP_t> Plm, dArray<int> Indx, 
 unsigned int Lmax, unsigned int stride, PII lohi, cuFP_t rmax)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;
    int npart = i + lohi.first;

    if (npart < lohi.second) {

#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
    
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
	if (Afac._v[i]<0.0 or Afac._v[i]>1.0) printf("spherical coord kernel--off grid: x=%f r=%f\n", xi, r);
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
(dArray<cuFP_t> coef, dArray<cuFP_t> tvar, dArray<cuFP_t> work,
 dArray<cuFP_t> used, dArray<cudaTextureObject_t> tex,
 dArray<cuFP_t> Mass, dArray<cuFP_t> Afac, dArray<cuFP_t> Phi,
 dArray<cuFP_t> Plm,  dArray<int> Indx,  int stride, 
 int l, int m, unsigned Lmax, unsigned int nmax, cuFP_t norm,
 PII lohi, bool compute)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;
  const int N     = lohi.second - lohi.first;

  const cuFP_t fac0 = -4.0*M_PI;

  for (int str=0; str<stride; str++) {

    int i     = tid*stride + str;
    int npart = i + lohi.first;

    if (npart < lohi.second) {

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
		      ) * p0 * plm[Ilm(l, m)] * fac0 * norm;
	  
	  coef._v[(2*n+0)*N + i] = v * cosp * mass;
	  coef._v[(2*n+1)*N + i] = v * sinp * mass;

	  // Load work space
	  //
	  if (compute and tvar._s>0) {
	    if (sphAcov) tvar._v[n*N + i   ] = v * mass;
	    else         work._v[i*nmax + n] = thrust::complex<double>(v*cosp, v*sinp);
	  }
	  
#ifdef BOUNDS_CHECK
	  if ((2*n+0)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	  if ((2*n+1)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	}

	if (compute and not sphAcov and tvar._s>0) {

	  int c = 0;
	  // Variance computation
	  for (int r=0; r<nmax; r++) {
	    for (int s=r; s<nmax; s++) {
	      thrust::complex<double> crs = work._v[i*nmax + r] * thrust::conj(work._v[i*nmax + s]);
	      tvar._v[(2*c+0)*N + i] = crs.real() * mass;
	      tvar._v[(2*c+1)*N + i] = crs.imag() * mass;
	      c++;
	    }
	  }
	  // Mean computation
	  for (int r=0; r<nmax; r++) {
	    tvar._v[(2*c+0)*N + i] = work._v[i*nmax + r].real() * mass;
	    tvar._v[(2*c+1)*N + i] = work._v[i*nmax + r].imag() * mass;
	    c++;
	  }
	  if (c != nmax*(nmax+3)/2) printf("out of bounds: wrong c [k]\n");
	}
      } else {
	// No contribution from off-grid particles
	for (int n=0; n<nmax; n++) {
	  coef._v[(2*n+0)*N + i] = 0.0;
	  coef._v[(2*n+1)*N + i] = 0.0;
	}

	if (compute and tvar._s>0) {
	  if (sphAcov) {
	    for (int n=0; n<nmax; n++) {
	      tvar._v[n*N + i] = 0.0;
	    }
	  } else {
	    int c = 0;
	    for (int r=0; r<nmax; r++) {
	      for (int s=r; s<nmax; s++) {
		tvar._v[(2*c+0)*N + i] = 0.0;
		tvar._v[(2*c+1)*N + i] = 0.0;
		c++;
	      }
	    }
	    for (int r=0; r<nmax; r++) {
	      tvar._v[(2*c+0)*N + i] = 0.0;
	      tvar._v[(2*c+1)*N + i] = 0.0;
	      c++;
	    }
	  }
	}
      }
    }
  }

}

__global__ void
forceKernel(dArray<cudaParticle> P, dArray<int> I, dArray<cuFP_t> coef,
	    dArray<cudaTextureObject_t> tex, dArray<cuFP_t> L1, dArray<cuFP_t> L2,
	    int stride, unsigned Lmax, unsigned int nmax, PII lohi, cuFP_t rmax)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {
      
#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      
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

      int ioff  = 0;
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
	printf("spherical forceKernel--off grid [0]: x=%f r=%f rs=%f a=%f b=%f\n", xi, r, rs, a0, b0);
      if (a1<0.0 or a1>1.0)
	printf("spherical forceKernel--off grid [1]: x=%f r=%f rs=%f a=%f b=%f\n", xi, r, rs, a1, b1);
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

	if (sphNO_L0  and l==0        ) continue;
	if (sphNO_L1  and l==1        ) continue;
	if (sphEVEN_L and (l/2)*2 != l) continue;

	cuFP_t fac1 = sqrt( (2.0*l + 1.0)/(4.0*M_PI) );

	cuFP_t ccos = 1.0;	// For recursion
	cuFP_t ssin = 0.0;

	// m loop
	//
	for (int m=0; m<=l; m++) {
	  
	  if (sphEVEN_M and (m/2)*2 != m) continue;
	  if (sphM0only and m != 0      ) continue;

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

	    pp_c +=  v * coef._v[indxC+n];
	    dp_c += dv * coef._v[indxC+n];
	    if (m>0) {
	      pp_s +=  v * coef._v[indxS+n];
	      dp_s += dv * coef._v[indxS+n];
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

	    if (ioff) {		// Factors for external multipole solution
	      cuFP_t facp  = pow(rmax/r0,(cuFP_t)(l+1));
	      cuFP_t facdp = -1.0/r0 * (l+1);
				// Apply the factors
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
	    
	    cuFP_t fac2 = M_SQRT2 * sqrt(numf/denf) * fac1;
	    
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

      p.acc[0] += -(potr*xx/r - pott*xx*zz/(r*r*r));
      p.acc[1] += -(potr*yy/r - pott*yy*zz/(r*r*r));
      p.acc[2] += -(potr*zz/r + pott*RR/(r*r*r));
      if (RR > FSMALL) {
	p.acc[0] +=  potp*yy/RR;
	p.acc[1] += -potp*xx/RR;
      }
      p.pot += potl;

#ifdef NAN_CHECK
      // Sanity check
      bool bad = false;
      for (int k=0; k<3; k++) {
	if (std::isnan(p.acc[k])) bad = true;
      }

      if (bad)  {
	printf("Force nan value: [%d] x=%f X=%f Y=%f Z=%f r=%f R=%f P=%f dP/dr=%f dP/dx=%f dP/dp=%f\n", p.indx, x, xx, yy, zz, r, RR, potl, potr, pott, potp);
	if (a<0.0 or a>1.0)  {
	  printf("Force nan value, no ioff: [%d] x=%f xi=%f dxi=%f a=%f i=%d\n",
		 p.indx, x, xi, dx, a, ind);
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


void SphericalBasis::cudaStorage::resize_coefs
(int nmax, int Lmax, int N, int gridSize, int stride,
 int sampT, bool pcavar, bool pcaeof, bool subsamp)
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
  if (pcavar) {
    T_covr.resize(sampT);
    for (int T=0; T<sampT; T++) {
      if (subsamp)
	T_covr[T].resize((Lmax+1)*(Lmax+2)*nmax);
      else
	T_covr[T].resize(2*(Lmax+1)*(Lmax+2)*(nmax*(nmax+3)/2));
    }
  }

  if (pcaeof or pcavar) {
    if (subsamp) {
      dN_tvar.resize(nmax*N);
      dc_tvar.resize(nmax*gridSize);
      dw_tvar.resize(nmax);
    } else {
      int csz = nmax*(nmax+3)/2;
      dW_tvar.resize(2*nmax*gridSize*BLOCK_SIZE*stride); // Volatile storage
      dN_tvar.resize(2*csz*N);
      dc_tvar.resize(2*csz*gridSize);
      dw_tvar.resize(2*csz);
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

void SphericalBasis::cuda_zero_coefs()
{
  auto cr = component->cuStream;
  
  // Resize output array
  //
  cuS.df_coef.resize((Lmax+1)*(Lmax+2)*nmax);
    
  // Zero output array
  //
  thrust::fill(thrust::cuda::par.on(cr->stream),
	       cuS.df_coef.begin(), cuS.df_coef.end(), 0.0);

  // Resize and zero PCA arrays
  //
  if (pcavar) {
    if (cuS.T_covr.size() != sampT) {
      cuS.T_covr.resize(sampT);
      for (int T=0; T<sampT; T++) {
	if (subsamp)
	  cuS.T_covr[T].resize((Lmax+1)*(Lmax+2)*nmax);
	else
	  cuS.T_covr[T].resize((Lmax+1)*(Lmax+2)*(nmax*(nmax+3)/2));
      }
      host_massT.resize(sampT);
    }
    
    for (int T=0; T<sampT; T++) {
      thrust::fill(thrust::cuda::par.on(cr->stream),
		   cuS.T_covr[T].begin(), cuS.T_covr[T].end(), 0.0);
    }
    
    thrust::fill(host_massT.begin(), host_massT.end(), 0.0);
  }

  if (pcaeof) {
    cuS.df_tvar.resize((Lmax+1)*(Lmax+2)/2*nmax*(nmax+1)/2);
    
    thrust::fill(thrust::cuda::par.on(cr->stream),
		 cuS.df_tvar.begin(), cuS.df_tvar.end(), 0.0);

    if (not pcavar) host_mass_tot = 0.0;
  }
}

void SphericalBasis::cudaStorage::resize_acc(int Lmax, int Nthread)
{
  // Space for Legendre coefficients 
  //
  /*
  if (plm1_d.capacity() < (Lmax+1)*(Lmax+2)/2*Nthread)
    plm1_d.reserve((Lmax+1)*(Lmax+2)/2*Nthread);
  
  if (plm2_d.capacity() < (Lmax+1)*(Lmax+2)/2*Nthread)
    plm2_d.reserve((Lmax+1)*(Lmax+2)/2*Nthread);
  */
  plm1_d.resize((Lmax+1)*(Lmax+2)/2*Nthread);
  plm2_d.resize((Lmax+1)*(Lmax+2)/2*Nthread);
}


void SphericalBasis::determine_coefficients_cuda(bool compute)
{
  // Only do this once but copying mapping coefficients and textures
  // must be done every time
  //
  if (initialize_cuda_sph) {
    initialize_cuda();
    initialize_cuda_sph = false;
  }

  // Copy coordinate mapping constants to device
  //
  initialize_mapping_constants();

  // Copy texture objects to device
  //
  t_d = tex;

  std::cout << std::scientific;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, component->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  auto cr = component->cuStream;

  // This will stay fixed for the entire run
  //
  host_coefs.resize((Lmax+1)*(Lmax+1)*nmax);

  if (pcavar) {
    host_massT.resize(sampT);
  }

  // Center assignment to symbol data
  //
  std::vector<cuFP_t> ctr;
  for (auto v : component->getCenter(Component::Local | Component::Centered))
    ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(sphCen, &ctr[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphCen");
  
  // VERBOSE diagnostic output on first call
  //
  static bool firstime = true;
  
  if (firstime and myid==0 and VERBOSE>4) {
    testConstantsSph<<<1, 1, 0, cr->stream>>>();
    cudaDeviceSynchronize();
    cuda_check_last_error_mpi("cudaDeviceSynchronize", __FILE__, __LINE__, myid);
    firstime = false;
  }
  
  // Zero counter and coefficients
  //
  thrust::fill(host_coefs.begin(), host_coefs.end(), 0.0);

  if (compute) {
    if (pcavar)
      thrust::fill(host_massT.begin(), host_massT.end(), 0.0);
    else
      host_mass_tot = 0.0;
  }

  // Zero out coefficient storage
  //
  cuda_zero_coefs();

  // Get sorted particle range for mlevel
  //
  PII lohi = component->CudaGetLevelRange(mlevel, mlevel), cur;
  
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
    static int debug_max_count = 10;
    static int debug_cur_count = 0;
    if (debug_cur_count++ < debug_max_count) {
      std::cout << std::endl
		<< "** cudaSphericalBasis coefficients" << std::endl
		<< "** N      = " << N          << std::endl
		<< "** Npacks = " << Npacks     << std::endl
		<< "** Stride = " << stride     << std::endl
		<< "** Block  = " << BLOCK_SIZE << std::endl
		<< "** Grid   = " << gridSize   << std::endl
		<< "** Xcen   = " << ctr[0]     << std::endl
		<< "** Ycen   = " << ctr[1]     << std::endl
		<< "** Zcen   = " << ctr[2]     << std::endl
		<< "**" << std::endl;
    }
#endif
    
    // Resize storage as needed
    //
    cuS.resize_coefs(nmax, Lmax, N, gridSize, stride,
		     sampT, pcavar, pcaeof, subsamp);
      
    // Shared memory size for the reduction
    //
    int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
    
    // Compute the coordinate transformation
    // 
    coordKernel<<<gridSize, BLOCK_SIZE, 0, cr->stream>>>
      (toKernel(cr->cuda_particles), toKernel(cr->indx1),
       toKernel(cuS.m_d), toKernel(cuS.a_d), toKernel(cuS.p_d),
       toKernel(cuS.plm1_d), toKernel(cuS.i_d),
       Lmax, stride, cur, rmax);
    
    // Compute the coefficient contribution for each order
    //
    int psize = nmax;
    int osize = nmax*2;
    int xsize = nmax*(nmax+3)/2;
    int vsize = 2*xsize;
    
    auto beg  = cuS.df_coef.begin();
    auto begV = cuS.df_tvar.begin();

    std::vector<thrust::device_vector<cuFP_t>::iterator> bm;
    if (pcavar) {
      for (int T=0; T<sampT; T++) {
	bm.push_back(cuS.T_covr[T].begin());
      }
    }
    thrust::fill(cuS.u_d.begin(), cuS.u_d.end(), 0.0);

    for (int l=0; l<=Lmax; l++) {
      for (int m=0; m<=l; m++) {
	cuFP_t ft = factorial(l, m);

	// Compute the contribution to the coefficients from each
	// particle
	//
	coefKernel<<<gridSize, BLOCK_SIZE, 0, cr->stream>>>
	  (toKernel(cuS.dN_coef), toKernel(cuS.dN_tvar), toKernel(cuS.dW_tvar),
	   toKernel(cuS.u_d), toKernel(t_d), toKernel(cuS.m_d),
	   toKernel(cuS.a_d), toKernel(cuS.p_d), toKernel(cuS.plm1_d),
	   toKernel(cuS.i_d), stride, l, m, Lmax, nmax, ft, cur, compute);
	
	// Begin the reduction per grid block [perhaps this should use
	// a stride?]
	//
	unsigned int gridSize1 = N/BLOCK_SIZE;
	if (N > gridSize1*BLOCK_SIZE) gridSize1++;

	reduceSum<cuFP_t, BLOCK_SIZE>
	  <<<gridSize1, BLOCK_SIZE, sMemSize, cr->stream>>>
	  (toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), osize, N);
      
	// Finish the reduction for this order in parallel
	//
	thrust::counting_iterator<int> index_begin(0);
	thrust::counting_iterator<int> index_end(gridSize1*osize);

	thrust::reduce_by_key
	  (
	   thrust::cuda::par.on(cr->stream),
	   thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
	   thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
	   cuS.dc_coef.begin(), thrust::make_discard_iterator(), cuS.dw_coef.begin()
	   );

	thrust::transform(thrust::cuda::par.on(cr->stream),
			  cuS.dw_coef.begin(), cuS.dw_coef.end(),
			  beg, beg, thrust::plus<cuFP_t>());
	
	thrust::advance(beg, osize);
	
	if (compute) {
	  
	  if (pcavar) {
	    
	    int sN = N/sampT;
	    int nT = sampT;
	    
	    if (sN==0) {	// Fail-safe underrun
	      sN = 1;
	      nT = N;
	    }
	    
	    for (int T=0; T<nT; T++) {
	      int k = sN*T;	// Starting position
	      int s = sN;
	      if (T==sampT-1) s = N - k;
	      
	      // Count and mass accumulation
	      //
	      if (l==0 and m==0) {
		
		// Sample count
		countT1[T] += s;

		auto mbeg = cuS.u_d.begin();
		auto mend = mbeg;
		thrust::advance(mbeg, sN*T);
		if (T<sampT-1) thrust::advance(mend, sN*(T+1));
		else mend = cuS.u_d.end();
		
		// Accumulated mass
		host_massT[T] += thrust::reduce(mbeg, mend);
	      }

	      // Begin the reduction per grid block
	      //
	      /* A reminder to consider implementing strides in reduceSum */
	      /*
		unsigned int stride1   = s/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
		unsigned int gridSize1 = s/BLOCK_SIZE/stride1;
		
		if (s > gridSize1*BLOCK_SIZE*stride1) gridSize1++;
	      */
		
	      unsigned int gridSize1 = s/BLOCK_SIZE;
	      if (s > gridSize1*BLOCK_SIZE) gridSize1++;

	      // Sum reduce into gridsize1 blocks taking advantage of
	      // GPU warp structure
	      //
	      if (subsamp) {
		
		// Mean computation only
		//
		reduceSumS<cuFP_t, BLOCK_SIZE>
		  <<<gridSize1, BLOCK_SIZE, sMemSize, cr->stream>>>
		  (toKernel(cuS.dc_tvar), toKernel(cuS.dN_tvar), psize, N, k, k+s);
		
		// Finish the reduction for this order in parallel
		//
		thrust::counting_iterator<int> index_begin(0);
		thrust::counting_iterator<int> index_end(gridSize1*psize);

		thrust::reduce_by_key
		  (
		   thrust::cuda::par.on(cr->stream),
		   thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
		   thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
		   cuS.dc_tvar.begin(), thrust::make_discard_iterator(), cuS.dw_tvar.begin()
		   );

	      
		thrust::transform(thrust::cuda::par.on(cr->stream),
				  cuS.dw_tvar.begin(), cuS.dw_tvar.end(),
				  bm[T], bm[T], thrust::plus<cuFP_t>());
	      
		thrust::advance(bm[T], psize);

	      } else {
		
		// Variance and mean computation
		//
		reduceSumS<cuFP_t, BLOCK_SIZE>
		  <<<gridSize1, BLOCK_SIZE, sMemSize, cr->stream>>>
		  (toKernel(cuS.dc_tvar), toKernel(cuS.dN_tvar), vsize, N, k, k+s);
		
		// Finish the reduction for this order in parallel
		//
		thrust::counting_iterator<int> indx2_begin(0);
		thrust::counting_iterator<int> indx2_end(gridSize1*vsize);
		
		thrust::reduce_by_key
		  (
		   thrust::cuda::par.on(cr->stream),
		   thrust::make_transform_iterator(indx2_begin, key_functor(gridSize1)),
		   thrust::make_transform_iterator(indx2_end,   key_functor(gridSize1)),
		   cuS.dc_tvar.begin(), thrust::make_discard_iterator(), cuS.dw_tvar.begin()
		   );
		
		thrust::transform(thrust::cuda::par.on(cr->stream),
				  cuS.dw_tvar.begin(), cuS.dw_tvar.end(),
				  bm[T], bm[T], thrust::plus<cuFP_t>());
		
		thrust::advance(bm[T], vsize);
	      }
	    }
	  }
	  // END: pcavar
	    
	  // Reduce EOF variance
	  //
	  if (pcaeof) {
	    
	    reduceSum<cuFP_t, BLOCK_SIZE>
	      <<<gridSize1, BLOCK_SIZE, sMemSize, cr->stream>>>
	      (toKernel(cuS.dc_tvar), toKernel(cuS.dN_tvar), vsize, N);
      
	    // Finish the reduction for this order in parallel
	    //
	    thrust::counting_iterator<int> index_begin(0);
	    thrust::counting_iterator<int> index_end(gridSize1*vsize);
	    
	    thrust::reduce_by_key
	      (
	       thrust::cuda::par.on(cr->stream),
	       thrust::make_transform_iterator(index_begin, key_functor(gridSize1)),
	       thrust::make_transform_iterator(index_end,   key_functor(gridSize1)),
	       cuS.dc_tvar.begin(), thrust::make_discard_iterator(), cuS.dw_tvar.begin()
	       );
	    
	    thrust::transform(thrust::cuda::par.on(cr->stream),
			      cuS.dw_tvar.begin(), cuS.dw_tvar.end(),
			      begV, begV, thrust::plus<cuFP_t>());
	    
	    thrust::advance(begV, vsize);

	    if (not pcavar and l==0 and m==0) {
	      auto mbeg = cuS.u_d.begin();
	      auto mend = cuS.u_d.end();
	      host_mass_tot += thrust::reduce(mbeg, mend);
	    }
	    
	  } // END: pcaeof
	} // END: compute
      } // END: m-loop
    } // END: l-loop

    // Compute number and total mass of particles used in coefficient
    // determination
    //
    thrust::sort(thrust::cuda::par.on(cr->stream), cuS.m_d.begin(), cuS.m_d.end());

    // Call the kernel on a single thread
    // 
    thrust::device_vector<cuFP_t>::iterator it;

    // Workaround for: https://github.com/NVIDIA/thrust/pull/1104
    //
    if (thrust_binary_search_workaround) {
      cudaStreamSynchronize(cr->stream);
      cuda_check_last_error_mpi("cudaStreamSynchronize", __FILE__, __LINE__, myid);
      it = thrust::lower_bound(cuS.m_d.begin(), cuS.m_d.end(), 0.0);
    } else {
      it = thrust::lower_bound(thrust::cuda::par.on(cr->stream),
			       cuS.m_d.begin(), cuS.m_d.end(), 0.0);
    }
    
    use[0] += thrust::distance(it, cuS.m_d.end());
  }

  // Copy back coefficient data from device and load the host
  //
  thrust::host_vector<cuFP_t> ret = cuS.df_coef;
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
  
  if (Ntotal == 0) {
    return;
  }

}

void SphericalBasis::determine_acceleration_cuda()
{
  // Only do this once but copying mapping coefficients and textures
  // must be done every time
  //
  if (initialize_cuda_sph) {
    initialize_cuda();
    initialize_cuda_sph = false;
  }

  // Copy coordinate mapping constants to device
  //
  initialize_mapping_constants();

  // Copy texture objects to device
  //
  t_d = tex;

  std::cout << std::scientific;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cC->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  // Stream structure iterators
  //
  auto cr = cC->cuStream;

  // Assign expansion center
  //
  std::vector<cuFP_t> ctr;
  for (auto v : component->getCenter(Component::Local | Component::Centered))
    ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(sphCen, &ctr[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphCen");
  
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

    unsigned int Nthread = gridSize*BLOCK_SIZE;

#ifdef VERBOSE_CTR
    static int debug_max_count = 10;
    static int debug_cur_count = 0;
    if (debug_cur_count++ < debug_max_count) {
      std::cout << std::endl
		<< "** cudaSphericalBasis acceleration" << std::endl
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
      
    cuS.resize_acc(Lmax, Nthread);

    // Shared memory size for the reduction
    //
    int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
    
    // Do the work
    //
    forceKernel<<<gridSize, BLOCK_SIZE, sMemSize, cr->stream>>>
      (toKernel(cr->cuda_particles), toKernel(cr->indx1),
       toKernel(dev_coefs), toKernel(t_d),
       toKernel(cuS.plm1_d), toKernel(cuS.plm2_d),
       stride, Lmax, nmax, lohi, rmax);
  }
}

void SphericalBasis::HtoD_coefs(const std::vector<VectorP>& expcoef)
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
      for (int n=0; n<nmax; n++) {
	host_coefs[Ilmn(l, m, 'c', n, nmax)] = (*expcoef[loffset+moffset])[n];
	if (m>0) host_coefs[Ilmn(l, m, 's', n, nmax)] = (*expcoef[loffset+moffset+1])[n];
      }

      if (m>0) moffset += 2;
      else     moffset += 1;
    }
  }

  dev_coefs = host_coefs;
}


void SphericalBasis::DtoH_coefs(std::vector<VectorP>& expcoef)
{
  // l loop
  //
  for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {

    // m loop
    //
    for (int m=0, moffset=0; m<=l; m++) {
	
      // n loop
      //
      for (int n=0; n<nmax; n++) {
	(*expcoef0[0][loffset+moffset])[n] = host_coefs[Ilmn(l, m, 'c', n, nmax)];
	if (m>0) (*expcoef0[0][loffset+moffset+1])[n] = host_coefs[Ilmn(l, m, 's', n, nmax)];
      }

      if (m>0) moffset += 2;
      else     moffset += 1;
    }
  }

  if (compute) {

    if (pcavar) {
      for (int T=0; T<sampT; T++) {
	massT1[T] += host_massT[T];
	muse1 [0] += host_massT[T];
      }
    } else {
      muse1[0] += host_mass_tot;
    }


    if (pcavar) {
      
      // T loop
      //
      for (int T=0; T<sampT; T++) {
	
	if (subsamp) {

	  thrust::host_vector<cuFP_t> retV = cuS.T_covr[T];
	  
	  int offst = 0;
	  int psize = nmax;
	  
	  // l loop
	  //
	  for (int l=0, loffset=0; l<=Lmax; loffset+=(l+1), l++) {
	    
	    // m loop
	    //
	    for (int m=0; m<=l; m++) {
	      
	      // n loop
	      //
	      for (int n=0; n<nmax; n++) {
		(*expcoefT1[T][loffset+m])[n] += retV[n + offst];
	      }
	      
	      offst += psize;
	    }
	  }

	} else {

	  thrust::host_vector<cuFP_t> retM = cuS.T_covr[T];
	
	  int vffst = 0;
	  int vsize = 2*nmax*(nmax+3)/2;
	  
	  // l loop
	  //
	  for (int l=0, loffset=0; l<=Lmax; loffset+=(l+1), l++) {
	    
	    // m loop
	    //
	    for (int m=0; m<=l; m++) {

	      // Variance assignment
	      //
	      int c = 0;	// Internal index 
	      for (int n=0; n<nmax; n++) {
	      
		for (int o=n; o<nmax; o++) {
		  std::complex<double> zret(retM[c + vffst], retM[c + vffst + 1]);
		  // Diagonal and upper diagonal
		  (*expcoefM1[T][loffset+m])(n, o) += zret;

		  // Below the diagonal
		  if (o!=n) (*expcoefM1[T][loffset+m])(o, n) += std::conj(zret);
		  c += 2;
		}
	      }

	      // Mean assignment
	      //
	      for (int n=0; n<nmax; n++) {
		std::complex<double> zret(retM[c + vffst], retM[c + vffst + 1]);
		(*expcoefT1[T][loffset+m])[n] += zret;
		c += 2;
	      }
	      
	      if (c != nmax*(nmax+3))
		std::cout << "out of bounds: wrong c [h]" << std::endl;

	      vffst += vsize;
	    }
	  }
	}
      }

      // EOF variance computation
      //
      if (pcaeof) {
	thrust::host_vector<cuFP_t> retV = cuS.df_tvar;
	int csz = nmax*(nmax+1)/2;
	int Ldim = (Lmax + 1)*(Lmax + 2)/2;
	if (retV.size() == Ldim*csz) {
	  for (int l=0; l<Ldim; l++) {
	    int c = 0;
	    for (size_t j=0; j<nmax; j++) {
	      for (size_t k=j; k<nmax; k++) {
		(*tvar[l])(j, k) += retV[csz*l + c];
		if (j!=k) (*tvar[l])(k, j) += retV[csz*l + c];
		c++;
	      }
	    }
	  }
	}
      }
    }
  }

  // Set to false for production
  //   |
  //   v
  if (false) coef_cuda_compare();
}


void SphericalBasis::coef_cuda_compare()
{

  // DEBUG, only useful for CUDAtest branch
  //

  constexpr bool compareC = true; // Compare host and device
  constexpr bool compareB = true;  // Tabulate largest errors

  if (fabs(host_coefs[Ilmn(0, 0, 'c', 0, nmax)])==0.0) return;

  if (compareC) {
    std::cout << std::string(2*4+4*20, '-') << std::endl
	      << "---- Spherical T=" << tnow << " level=" << mlevel
	      << std::endl << std::string(2*4+4*20, '-') << std::endl;
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
	      << "---- Spherical T=" << tnow << " level=" << mlevel
	      << std::endl << std::string(2*4+20, '-') << std::endl;
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
    cuFP_t b = (*expcoef0[0][0])[n];
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
    cuFP_t b = (*expcoef0[0][1])[n];
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
    cuFP_t b = (*expcoef0[0][2])[n];
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
    cuFP_t b = (*expcoef0[0][3])[n];
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
    cuFP_t b = (*expcoef0[0][4])[n];
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
    cuFP_t b = (*expcoef0[0][5])[n];
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
    cuFP_t b = (*expcoef0[0][6])[n];
    if (compareC) {
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << i
		<< std::setw(20) << a
		<< std::setw(20) << b
		<< std::setw(20) << a - b;
      if (fabs(*cmax) > 0.0)
	std::cout << std::setw(20) << (a - b)/fabs(*cmax);
      else
	std::cout << std::setw(20) << 0.0;
      std::cout << std::endl;
    } else {
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << i
		<< std::setw(20) << a
		<< std::endl;
    }
  }
  
  std::cout << std::string(2*4+4*20, '-') << std::endl;


  //
  // TEST comparison of coefficients for debugging
  //
  if (compareB) {

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
	  for (int n=0; n<nmax; n++) {
	    elem.l = l;
	    elem.m = m;
	    elem.n = n;
	    elem.cs = 'c';
	    elem.d = (*expcoef0[0][loffset+moffset])[n];
	    elem.f = host_coefs[Ilmn(l, m, 'c', n, nmax)];
	    
	    double test = fabs(elem.d - elem.f);
	    if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);
	    
	    compare.insert(std::make_pair(test, elem));
	    
	    out << std::setw( 5) << l
		<< std::setw( 5) << m
		<< std::setw( 5) << n
		<< std::setw( 5) << 'c'
		<< std::setw( 5) << Ilmn(l, m, 'c', n, nmax)
		<< std::setw(14) << elem.d
		<< std::setw(14) << elem.f
		<< std::endl;
	  }
	  
	  moffset++;
	}
	else {
	  for (int n=0; n<nmax; n++) {
	    elem.l = l;
	    elem.m = m;
	    elem.n = n;
	    elem.cs = 'c';
	    elem.d = (*expcoef0[0][loffset+moffset])[n];
	    elem.f = host_coefs[Ilmn(l, m, 'c', n, nmax)];

	    out << std::setw( 5) << l
		<< std::setw( 5) << m
		<< std::setw( 5) << n
		<< std::setw( 5) << 'c'
		<< std::setw( 5) << Ilmn(l, m, 'c', n, nmax)
		<< std::setw(14) << elem.d
		<< std::setw(14) << elem.f
		<< std::endl;

	    double test = fabs(elem.d - elem.f);
	    if (fabs(elem.d)>1.0e-12) test /= fabs(elem.d);

	    compare.insert(std::make_pair(test, elem));
	  }
	  for (int n=0; n<nmax; n++) {
	    elem.l = l;
	    elem.m = m;
	    elem.n = n;
	    elem.cs = 's';
	    elem.d = (*expcoef0[0][loffset+moffset+1])[n];
	    elem.f = host_coefs[Ilmn(l, m, 's', n, nmax)];

	    out << std::setw( 5) << l
		<< std::setw( 5) << m
		<< std::setw( 5) << n
		<< std::setw( 5) << 's'
		<< std::setw( 5) << Ilmn(l, m, 's', n, nmax)
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
	      << "---- Spherical coefficients T=" << tnow
	      << " level=" << mlevel << std::endl
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


void SphericalBasis::multistep_update_cuda()
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

  //! Sort the device vector by level changes
#ifdef VERBOSE_TIMING
  auto start0 = std::chrono::high_resolution_clock::now();
  auto start  = std::chrono::high_resolution_clock::now();
#endif

  auto chg = component->CudaSortLevelChanges();

#ifdef SPH_UPDATE_TABLE
  multistep_add_debug(chg);
#endif

#ifdef VERBOSE_TIMING
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::micro> duration = finish - start;
  std::cout << "Time in level sort=" << duration.count()*1.0e-6 << std::endl;
#endif

  // Zero out coefficient storage
  //
  cuda_zero_coefs();

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, component->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  auto cs = component->cuStream;

#ifdef VERBOSE_TIMING
  double coord = 0.0, coefs = 0.0, reduc = 0.0;
#endif
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
		<< "] Adjust sphere: Ntotal=" << Ntotal << " Npacks=" << Npacks
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

	// Resize storage as needed
	//
	cuS.resize_coefs(nmax, Lmax, N, gridSize, stride,
			 sampT, pcavar, pcaeof, subsamp);
	
	// Shared memory size for the reduction
	//
	int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
	
	// Compute the coordinate transformation
	// 
#ifdef VERBOSE_TIMING
	start = std::chrono::high_resolution_clock::now();
#endif
	coordKernel<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
	  (toKernel(cs->cuda_particles), toKernel(cs->indx2),
	   toKernel(cuS.m_d), toKernel(cuS.a_d), toKernel(cuS.p_d),
	   toKernel(cuS.plm1_d), toKernel(cuS.i_d),
	   Lmax, stride, cur, rmax);

#ifdef VERBOSE_TIMING
	finish = std::chrono::high_resolution_clock::now();
	duration = finish - start;
	coord += duration.count()*1.0e-6;
#endif    
	// Compute the coefficient contribution for each order
	//
	int osize = nmax*2;
	auto beg  = cuS.df_coef.begin();
      
	thrust::fill(cuS.u_d.begin(), cuS.u_d.end(), 0.0);

	for (int l=0; l<=Lmax; l++) {
	  for (int m=0; m<=l; m++) {
	    cuFP_t ft = factorial(l, m);

	    // Compute the contribution to the
	    // coefficients from each particle
	    //
#ifdef VERBOSE_TIMING
	    start = std::chrono::high_resolution_clock::now();
#endif
	    coefKernel<<<gridSize, BLOCK_SIZE, 0, cs->stream>>>
	      (toKernel(cuS.dN_coef),
	       toKernel(cuS.dN_tvar), toKernel(cuS.dW_tvar),
	       toKernel(cuS.u_d), toKernel(t_d), toKernel(cuS.m_d),
	       toKernel(cuS.a_d), toKernel(cuS.p_d), toKernel(cuS.plm1_d),
	       toKernel(cuS.i_d), stride, l, m, Lmax, nmax, ft, cur, false);

#ifdef VERBOSE_TIMING
	    finish = std::chrono::high_resolution_clock::now();
	    duration = finish - start;
	    coefs += duration.count()*1.0e-6;
	    start = std::chrono::high_resolution_clock::now();
#endif	  
	    // Begin the reduction per grid block
	    // [perhaps this should use a stride?]
	    //
	    unsigned int gridSize1 = N/BLOCK_SIZE;
	    if (N > gridSize1*BLOCK_SIZE) gridSize1++;

	    reduceSum<cuFP_t, BLOCK_SIZE>
	      <<<gridSize1, BLOCK_SIZE, sMemSize, cs->stream>>>
	      (toKernel(cuS.dc_coef), toKernel(cuS.dN_coef), osize, N);
	    
	    // Finish the reduction for this order
	    // in parallel
	    //
	    thrust::counting_iterator<int> index_begin(0);
	    thrust::counting_iterator<int> index_end(gridSize1*osize);

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

#ifdef VERBOSE_TIMING
	    finish = std::chrono::high_resolution_clock::now();
	    duration = finish - start;
	    reduc += duration.count()*1.0e-6;
#endif	    
	    thrust::advance(beg, osize);
	  }
	  // END: m-loop
	}
	// END: l-loop
      }
      // END: bunches

      // Copy back coefficient data from device and load the host
      //
      thrust::host_vector<cuFP_t> ret = cuS.df_coef;

      // Decrement current level and increment new level using the
      // update matrices
      //
      for (int l=0, loffset=0, offst=0; l<=Lmax; l++) {
	for (int m=0; m<=l; m++) {
	  for (size_t n=0; n<nmax; n++) {
	    differ1[0][olev](loffset, n) -= ret[2*n+offst];
	    differ1[0][nlev](loffset, n) += ret[2*n+offst];
	    if (m>0) {
	      differ1[0][olev](loffset+1, n) -= ret[2*n+1+offst];
	      differ1[0][nlev](loffset+1, n) += ret[2*n+1+offst];
	    }
	  }

	  // Update the offset into the device coefficient array
	  offst += nmax*2;

	  // Update the offset into the host coefficient matrix
	  if (m>0) loffset += 2;
	  else     loffset += 1;
	}
      }
      // END: assign differences
    }
    // END: to new level loop
  }
  // END: from prev level loop

#ifdef VERBOSE_TIMING
  std::cout << "Time in coord=" << coord << std::endl;
  std::cout << "Time in coefs=" << coefs << std::endl;
  std::cout << "Time in reduc=" << reduc << std::endl;
  auto finish0 = std::chrono::high_resolution_clock::now();
  duration = finish0 - start0;
  std::cout << "Total adjust =" << duration.count()*1.0e-6 << std::endl;
#endif
  // DONE
}

void SphericalBasis::destroy_cuda()
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
