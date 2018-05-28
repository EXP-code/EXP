#include <Component.H>
#include <SphericalBasis.H>
#include <cudaReduce.cuH>

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
(dArray<cuFP_t> coef, dArray<cudaTextureObject_t> tex,
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

void SphericalBasis::determine_coefficients_cuda()
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

  // Sort particles and get coefficient size
  //
  PII lohi = cC->CudaSortByLevel(mlevel, mlevel);

  // Zero out coefficients
  //
  host_coefs.resize((Lmax+1)*(Lmax+1)*nmax);

  // Compute grid
  //
  unsigned int N         = lohi.second - lohi.first;
  unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  unsigned int gridSize  = N/BLOCK_SIZE/stride;

  if (N == 0) {
    use[0] = 0.0;
    return;
  }

  if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

  std::vector<cuFP_t> ctr;
  for (auto v : cC->getCenter(Component::Local | Component::Centered)) ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(sphCen, &ctr[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphCen");

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


  // Create space for coefficient reduction
  //
  dN_coef.resize(2*nmax*N);
  dc_coef.resize(2*nmax*gridSize);
  df_coef.resize(2*nmax);

  // Space for Legendre coefficients 
  //
  plm_d.resize((Lmax+1)*(Lmax+2)/2*N);

  // Space for coordinates
  //
  r_d.resize(N);
  m_d.resize(N);
  a_d.resize(N);
  p_d.resize(N);
  i_d.resize(N);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);

  // For debugging
  //
  static bool firstime = false;

  if (firstime) {
    testConstants<<<1, 1>>>();
    cudaDeviceSynchronize();
    firstime = false;
  }

  thrust::counting_iterator<int> index_begin(0);
  thrust::counting_iterator<int> index_end(gridSize*2*nmax);

  // Do the work
  //
				// Compute the coordinate
				// transformation
				// 
  coordKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(cC->cuda_particles),
     toKernel(m_d), toKernel(a_d), toKernel(p_d), toKernel(plm_d),
     toKernel(i_d), Lmax, stride, lohi, rmax);

				// Compute the coefficient
				// contribution for each order
  int osize = nmax*2;		//
  for (int l=0; l<=Lmax; l++) {
    for (int m=0; m<=l; m++) {
				// Compute the contribution to the
				// coefficients from each particle
				//
      coefKernel<<<gridSize, BLOCK_SIZE>>>
	(toKernel(dN_coef), toKernel(t_d), toKernel(m_d),
	 toKernel(a_d), toKernel(p_d), toKernel(plm_d), toKernel(i_d),
	 stride, l, m, Lmax, nmax, lohi);

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

				// Assign reduced output to
				// coefficient array
				//
      thrust::host_vector<cuFP_t> ret = df_coef;
      for (size_t j=0; j<nmax; j++) {
	host_coefs[Ilmn(l, m, 'c', j, nmax)] = ret[2*j];
	if (m>0) host_coefs[Ilmn(l, m, 's', j, nmax)] = ret[2*j+1];
      }
    }
  }

  // Compute number and total mass of particles used in coefficient
  // determination
  //
  thrust::sort(m_d.begin(), m_d.end());

  auto m_it = thrust::upper_bound(m_d.begin(), m_d.end(), 0.0);
  use[0]    = thrust::distance(m_it, m_d.end());
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

  unsigned int Nthread = gridSize*BLOCK_SIZE;

  std::vector<cuFP_t> ctr;
  for (auto v : cC->getCenter(Component::Local | Component::Centered)) ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(sphCen, &ctr[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphCen");

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

  // Space for Legendre coefficients 
  //
  plm1_d.resize((Lmax+1)*(Lmax+2)/2*Nthread);
  plm2_d.resize((Lmax+1)*(Lmax+2)/2*Nthread);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);

  // Do the work
  //
  forceKernel<<<gridSize, BLOCK_SIZE, sMemSize>>>
    (toKernel(cC->cuda_particles), toKernel(dev_coefs), toKernel(t_d),
     toKernel(plm1_d), toKernel(plm2_d), stride, Lmax, nmax, lohi, rmax,
     use_external);
}

void SphericalBasis::HtoD_coefs(const Matrix& expcoef)
{
  host_coefs.resize((Lmax+1)*(Lmax+1)*nmax);

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
	expcoef[loffset+moffset][n] = host_coefs[Ilmn(l, m, 'c', n-1, nmax)];
	if (m>0) expcoef[loffset+moffset+1][n] = host_coefs[Ilmn(l, m, 's', n-1, nmax)];
      }

      if (m>0) moffset += 2;
      else     moffset += 1;
    }
  }
}

void SphericalBasis::destroy_cuda()
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
    sout << "trying to free cuArray [" << i << "]";
    cuda_safe_call(cudaFreeArray(cuInterpArray[i]),
		     __FILE__, __LINE__, sout.str());
  }
    
  // std::cout << "cuda memory freed" << std::endl;
}

void SphericalBasis::host_dev_force_compare()
{
  // Copy from device
  cC->host_particles = cC->cuda_particles;
  
  std::streamsize ss = std::cout.precision();
  std::cout.precision(4);

  std::cout << std::string(16+14*8, '-') << std::endl
	    << std::setw(8)  << "Index"  << std::setw(8)  << "Level"
	    << std::setw(14) << "ax [d]" << std::setw(14) << "ay [d]"
	    << std::setw(14) << "az [d]" << std::setw(14) << "ax [h]"
	    << std::setw(14) << "ay [h]" << std::setw(14) << "az [h]"
	    << std::setw(14) << "|Del a|/|a|"
	    << std::setw(14) << "|a|"    << std::endl;

  // Compare first and last 5 from the device list
  //
  for (size_t i=0; i<5; i++) 
    {
      auto indx = cC->host_particles[i].indx;
      auto levl = cC->host_particles[i].level;
      
      std::cout << std::setw(8) << indx	<< std::setw(8) << levl;

      for (int k=0; k<3; k++)
	std::cout << std::setw(14) << cC->host_particles[i].acc[k];

      for (int k=0; k<3; k++)
	std::cout << std::setw(14) << cC->Particles()[indx].acc[k];

      cuFP_t diff = 0.0, norm = 0.0;
      for (int k=0; k<3; k++) {
	cuFP_t b  = cC->host_particles[i].acc[k];
	cuFP_t a  = cC->Particles()[indx].acc[k];
	diff += (a - b)*(a - b);
	norm += a*a;
      }
      std::cout << std::setw(14) << sqrt(diff/norm)
		<< std::setw(14) << sqrt(norm) << std::endl;
    }
  
  for (size_t j=0; j<5; j++) 
    {
      size_t i = cC->host_particles.size() - 5 + j;

      auto indx = cC->host_particles[i].indx;
      auto levl = cC->host_particles[i].level;

      std::cout << std::setw(8) << indx	<< std::setw(8) << levl;
      
      for (int k=0; k<3; k++)
	std::cout << std::setw(14) << cC->host_particles[i].acc[k];

      for (int k=0; k<3; k++)
	std::cout << std::setw(14) << cC->Particles()[indx].acc[k];

      cuFP_t diff = 0.0, norm = 0.0;
      for (int k=0; k<3; k++) {
	cuFP_t b  = cC->host_particles[i].acc[k];
	cuFP_t a  = cC->Particles()[indx].acc[k];
	diff += (a - b)*(a - b);
	norm += a*a;
      }
      std::cout << std::setw(14) << sqrt(diff/norm)
		<< std::setw(14) << sqrt(norm) << std::endl;
    }

  std::cout << std::string(16+14*8, '-') << std::endl;
  std::cout.precision(ss);
}
