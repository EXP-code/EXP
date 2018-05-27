#include <Component.H>
#include <SphericalBasis.H>
#include <cudaReduce.cuH>

// Define for debugging
//
// #define OFF_GRID_ALERT
// #define BOUNDS_CHECK
// #define VERBOSE
#define NAN_CHECK

// Machine constant for Legendre
//
const double DMINEPS=20.0*DBL_MIN;

// Global symbols for coordinate transformation in SphericalBasis
//
__device__ __constant__
double sphScale, sphRscale, sphHscale, sphXmin, sphXmax, sphDxi, sphCen[3];

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
void legendre_v(int lmax, double x, double* p)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;

  p[0] = pll = 1.0f;
  if (lmax > 0) {
    somx2 = sqrt( (1.0f - x)*(1.0f + x) );
    fact = 1.0f;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[Ilm(m, m)] = pll;
      fact += 2.0f;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p[Ilm(m, m)];
    p[Ilm(m+1, m)] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[Ilm(l, m)] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }
}


__host__ __device__
void legendre_v2(int lmax, double x, double* p, double* dp)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;

  p[0] = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[Ilm(m, m)] = pll;
      fact += 2.0;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p[Ilm(m, m)];
    p[Ilm(m+1, m)] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[Ilm(l, m)] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (1.0-fabs(x) < DMINEPS) {
    if (x>0) x =   1.0 - DMINEPS;
    else     x = -(1.0 - DMINEPS);
  }

  somx2 = 1.0/(x*x - 1.0);
  dp[0] = 0.0;
  for (l=1; l<=lmax; l++) {
    for (m=0; m<l; m++)
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
double cu_r_to_xi(double r)
{
  double ret;

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
double cu_xi_to_r(double xi)
{
  double ret;

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
double cu_d_xi_to_r(double xi)
{
  double ret;

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
  double z;

  cuda_safe_call(cudaMemcpyToSymbol(sphScale, &(z=scale), sizeof(double), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphScale");

  cuda_safe_call(cudaMemcpyToSymbol(sphRscale, &f.rscale, sizeof(double), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphRscale");

  cuda_safe_call(cudaMemcpyToSymbol(sphXmin,   &f.xmin,   sizeof(double), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphXmin");

  cuda_safe_call(cudaMemcpyToSymbol(sphXmax,   &f.xmax,   sizeof(double), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphXmax");

  cuda_safe_call(cudaMemcpyToSymbol(sphDxi,    &f.dxi,    sizeof(double), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphDxi");

  cuda_safe_call(cudaMemcpyToSymbol(sphNumr,   &f.numr,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphuNumr");

  cuda_safe_call(cudaMemcpyToSymbol(sphCmap,   &f.cmap,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying sphCmap");
}

__global__
void testTextureSph(dArray<cudaTextureObject_t> tex, int nmax)
{
  printf("** DEVICE 1d texture compare\n");
  for (int k=0; k<10; k++) {
    for (int i : {0, 1, sphNumr-2, sphNumr-1}) 
      printf("%5d %5d %13.7e\n", k, i, int2_as_double(tex1D<int2>(tex._v[k], i)));
  }
}

__global__ void coordKernel
(dArray<cudaParticle> in, dArray<double> mass, dArray<double> Afac,
 dArray<double> phi, dArray<double> Plm, dArray<int> Indx, 
 unsigned int Lmax, unsigned int stride, PII lohi, double rmax)
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
    
      double xx = p.pos[0] - sphCen[0];
      double yy = p.pos[1] - sphCen[1];
      double zz = p.pos[2] - sphCen[2];
      
      double r2 = (xx*xx + yy*yy + zz*zz);
      double r = sqrt(r2) + FSMALL;
      
#ifdef BOUNDS_CHECK
      if (i>=mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      mass._v[i] = -1.0;
      
      if (r<rmax) {
	
	mass._v[i] = p.mass;
	
	double costh = zz/r;
	phi._v[i] = atan2(yy,xx);
	
#ifdef BOUNDS_CHECK
	if (i>=phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	double *plm = &Plm._v[psiz*i];
	legendre_v(Lmax, costh, plm);

#ifdef BOUNDS_CHECK
	if (psiz*(i+1)>Plm._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	double x  = cu_r_to_xi(r);
	double xi = (x - sphXmin)/sphDxi;
	int indx = floor(xi);
	
	if (indx<0) indx = 0;
	if (indx>sphNumr-2) indx = sphNumr - 2;
	  
	Afac._v[i] = double(indx+1) - xi;
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
(dArray<double> coef, dArray<cudaTextureObject_t> tex,
 dArray<double> Mass, dArray<double> Afac, dArray<double> Phi,
 dArray<double> Plm, dArray<int> Indx,  int stride, 
 int l, int m, unsigned Lmax, unsigned int nmax, PII lohi)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;
  const unsigned int N = lohi.second - lohi.first;

  double fac0 = 4.0*M_PI;

  for (int istr=0; istr<stride; istr++) {

    int i = tid*stride + istr;

    if (i<N) {

      double mass = Mass._v[i];

#ifdef BOUNDS_CHECK
      if (i>=Mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif

      if (mass>0.0) {

#ifdef BOUNDS_CHECK
	if (i>=Phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	double phi  = Phi._v[i];
	double cosp = cos(phi*m);
	double sinp = sin(phi*m);
	
#ifdef BOUNDS_CHECK
	if (psiz*(i+1)>Plm._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif	
	double *plm = &Plm._v[psiz*i];
	
	// Do the interpolation
	//
	double a = Afac._v[i];
	double b = 1.0 - a;
	int ind = Indx._v[i];
	
#ifdef BOUNDS_CHECK
	if (i>=Afac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=Indx._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	for (int n=0; n<nmax; n++) {

	  double p0 =
	    a*int2_as_double(tex1D<int2>(tex._v[0], ind  )) +
	    b*int2_as_double(tex1D<int2>(tex._v[0], ind+1)) ;

	  int k = 1 + l*nmax + n;

#ifdef BOUNDS_CHECK
	  if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	  float v = (
		     a*int2_as_double(tex1D<int2>(tex._v[k], ind  )) +
		     b*int2_as_double(tex1D<int2>(tex._v[k], ind+1))
		     ) * p0 * plm[Ilm(l, m)] * Mass._v[i] * fac0;
	  
	  
	  coef._v[(2*n+0)*N + i] = v * cosp;
	  coef._v[(2*n+1)*N + i] = v * sinp;

#ifdef BOUNDS_CHECK
	  if ((2*n+0)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	  if ((2*n+1)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
	}
      }
    }
  }
}

__global__ void
forceKernel(dArray<cudaParticle> in, dArray<double> coef,
	    dArray<cudaTextureObject_t> tex, dArray<double> L1, dArray<double> L2,
	    int stride, unsigned Lmax, unsigned int nmax, PII lohi, double rmax,
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
      
      double xx = p.pos[0] - sphCen[0];
      double yy = p.pos[1] - sphCen[1];
      double zz = p.pos[2] - sphCen[2];
      
      if (fabs(xx) > 10.0*rmax or
	  fabs(yy) > 10.0*rmax or
	  fabs(zz) > 10.0*rmax) continue;

      double r2 = (xx*xx + yy*yy + zz*zz);
      double r  = sqrt(r2) + FSMALL;
      
      double costh = zz/r;
      double phi   = atan2(yy, xx);
      double RR    = xx*xx + yy*yy;
      
      double *plm1 = &L1._v[psiz*tid];
      double *plm2 = &L2._v[psiz*tid];
      
      const double FLIM = 1.0e-5;
      if (costh >   1.0 - FLIM ) costh =   1.0 - FLIM;
      if (costh < -(1.0 - FLIM)) costh = -(1.0 - FLIM);

      legendre_v2(Lmax, costh, plm1, plm2);

      int ioff = 0;
      double rs = r/sphScale;
      double r0 = 0.0;

      if (r>rmax) {
	ioff = 1;
	r0   = r;
	r    = rmax;
	rs   = r/sphScale;
      }

      double  x = cu_r_to_xi(rs);
      double xi = (x - sphXmin)/sphDxi;
      double dx = cu_d_xi_to_r(x)/sphDxi;
      int  ind = floor(xi);
      
      if (ind<1) ind = 1;
      if (ind>sphNumr-2) ind = sphNumr - 2;
      
      double a = (double)(ind+1) - xi;
#ifdef OFF_GRID_ALERT
      if (a<0.0 or a>1.0) printf("forceKernel: off grid: x=%f\n", xi);
#endif
      double b = 1.0 - a;
      
      // Do the interpolation for the prefactor potential
      //
      double pm1 = int2_as_double(tex1D<int2>(tex._v[0], ind-1));
      double p00 = int2_as_double(tex1D<int2>(tex._v[0], ind  ));
      double pp1 = int2_as_double(tex1D<int2>(tex._v[0], ind+1));

      // For force accumulation
      //
      double potl = 0.0;
      double potr = 0.0;
      double pott = 0.0;
      double potp = 0.0;

      // l loop
      //
      for (int l=0; l<Lmax; l++) {

	double fac1 = (2.0*l + 1.0)/(4.0*M_PI);

	// m loop
	//
	for (int m=0; m<=l; m++) {

	  int pindx = Ilm(l, m);

	  double Plm1 = plm1[pindx];
	  double Plm2 = plm2[pindx];
      
#ifdef NAN_CHECK
	  if (std::isnan(Plm1)) {
	    printf("Force isnan for Plm(%d, %d) ioff=%d, costh=%f, z=%f, r=%f\n", l, m, ioff, costh, zz, r);
	  }

	  if (std::isnan(Plm2)) {
	    printf("Force isnan for Plm2(%d, %d) ioff=%d costh=%f, z=%f, r=%f\n", l, m, ioff, costh, zz, r);
	  }
#endif	  

	  double pp_c = 0.0;
	  double dp_c = 0.0;
	  double pp_s = 0.0;
	  double dp_s = 0.0;
	  
	  int indxC = Ilmn(l, m, 'c', 0, nmax);
	  int indxS = Ilmn(l, m, 's', 0, nmax);

	  double cosp = cos(phi*m);
	  double sinp = sin(phi*m);

	  for (size_t n=0; n<nmax; n++) {
	
	    int k = 1 + l*nmax + n;
	
#ifdef BOUNDS_CHECK
	    if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif	
	    double um1 = int2_as_double(tex1D<int2>(tex._v[k], ind-1));
	    double u00 = int2_as_double(tex1D<int2>(tex._v[k], ind  ));
	    double up1 = int2_as_double(tex1D<int2>(tex._v[k], ind+1));
	    
	    double v = (a*u00 + b*up1)*(a*p00 + b*pp1);
	    
	    double dv =
	      dx * ( (b - 0.5)*um1*pm1 - 2.0*b*u00*p00 + (b + 0.5)*up1*pp1 );
	    
	    if (std::isnan(v))
	      printf("v tab nan: (%d, %d): a=%f b=%f p0=%f p1=%f u0=%f u1=%f\n", l, m, a, b, p00, pp1, u00, up1);

	    if (std::isnan(dv))
	      printf("dv tab nan: (%d, %d): a=%f b=%f pn=%f p0=%f pp=%f un=%f u0=%f up=%f\n", l, m, a, b, pm1, p00, pp1, um1, u00, up1);

	    pp_c +=  v * coef._v[indxC+n];
	    dp_c += dv * coef._v[indxC+n];
	    if (m>0) {
	      pp_s +=  v * coef._v[indxS+n];
	      dp_s += dv * coef._v[indxS+n];
	    }

	  } // END: n loop
	  
	  pp_c *= -1.0;
	  dp_c *= -1.0;
	  pp_s *= -1.0;
	  dp_s *= -1.0;

	  if (std::isnan(pp_c)) printf("pp_c eval nan: (%d, %d): r=%f r0=%f\n", l, m, r, r0);
	  if (std::isnan(dp_c)) printf("dp_c eval nan: (%d, %d): r=%f r0=%f\n", l, m, r, r0);
	  if (std::isnan(pp_s)) printf("pp_s eval nan: (%d, %d): r=%f r0=%f\n", l, m, r, r0);
	  if (std::isnan(dp_s)) printf("dp_s eval nan: (%d, %d): r=%f r0=%f\n", l, m, r, r0);

	  if (m==0) {

	    if (ioff) {
	      pp_c *= pow(rmax/r0, (double)(l+1));
	      dp_c  = -pp_c/r0 * (double)(l+1);
#ifdef NAN_CHECK
	      if (std::isnan(pp_c)) printf("Force nan [ioff]: l=%d, r=%f r0=%f\n", l, r, r0);
#endif
	    }
	    
	    potl += fac1 * pp_c * Plm1;
	    potr += fac1 * dp_c * Plm1;
	    pott += fac1 * pp_c * Plm2;
	    potp += 0.0;
	    
	  } else {

	    double cosm = cos(phi*m);
	    double sinm = sin(phi*m);
	    
	    if (ioff) {
	      double facp  = pow(rmax/r0,(double)(l+1));
	      double facdp = -facp/r0 * (l+1);
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
	    double numf = 1.0, denf = 1.0;
	    for (int i=2; i<=l+m; i++) {
	      if (i<=l-m) numf *= i; denf *= i;
	    }
	    
	    double fac2 = 2.0 * numf/denf * fac1;
	    
	    potl += fac2 * Plm1 * ( pp_c*cosm + pp_s*sinm);
	    potr += fac2 * Plm1 * ( dp_c*cosm + dp_s*sinm);
	    pott += fac2 * Plm2 * ( pp_c*cosm + pp_s*sinm);
	    potp += fac2 * Plm1 * (-pp_c*sinm + pp_s*cosm)*m;
	  }

	} // END: m loop

      } // END: l loop
				// Rescale
      potr /= sphScale*sphScale;
      potl /= sphScale;
      pott /= sphScale;
      potp /= sphScale;

      in._v[npart].acc[0] += -(potr*xx/r - pott*xx*zz/(r*r*r));
      in._v[npart].acc[1] += -(potr*yy/r - pott*yy*zz/(r*r*r));
      in._v[npart].acc[2] += -(potr*zz/r - pott*RR/(r*r*r));
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
  thrust::fill(host_coefs.begin(), host_coefs.end(), 0.0);

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

  std::vector<double> ctr;
  for (auto v : cC->getCenter(Component::Local | Component::Centered)) ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(sphCen, &ctr[0], sizeof(double)*3,
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
  thrust::device_vector<double> dN_coef(2*nmax*N);
  thrust::device_vector<double> dc_coef(2*nmax*gridSize);
  thrust::device_vector<double> df_coef(2*nmax);

  // Texture objects
  //
  thrust::device_vector<cudaTextureObject_t> t_d = tex;

  // Space for Legendre coefficients 
  //
  thrust::device_vector<double> plm_d((Lmax+1)*(Lmax+2)/2*N);
  thrust::device_vector<double> r_d(N), m_d(N), a_d(N), p_d(N);
  thrust::device_vector<int>   i_d(N);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(double);

  // For debugging
  //
  static bool firstime = true;

  if (firstime) {
    testConstants<<<1, 1>>>();
    cudaDeviceSynchronize();
    /*
    testTextureSph<<<1, 1>>>(toKernel(t_d), nmax);
    cudaDeviceSynchronize();
    */
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
      coefKernel<<<gridSize, BLOCK_SIZE>>>
	(toKernel(dN_coef), toKernel(t_d), toKernel(m_d),
	 toKernel(a_d), toKernel(p_d), toKernel(plm_d), toKernel(i_d),
	 stride, l, m, Lmax, nmax, lohi);

				// Begin the reduction per grid block
      				// 
      reduceSum<double, BLOCK_SIZE><<<gridSize, BLOCK_SIZE, sMemSize>>>
	(toKernel(dc_coef), toKernel(dN_coef), osize, N);
      
				// Finish the reduction for this order
				// in parallel
      thrust::reduce_by_key
	(
	 thrust::make_transform_iterator(index_begin, key_functor(gridSize)),
	 thrust::make_transform_iterator(index_end,   key_functor(gridSize)),
	 dc_coef.begin(), thrust::make_discard_iterator(), df_coef.begin()
	 );

      thrust::host_vector<double> ret = df_coef;
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

  std::vector<double> ctr;
  for (auto v : cC->getCenter(Component::Local | Component::Centered)) ctr.push_back(v);

  cuda_safe_call(cudaMemcpyToSymbol(sphCen, &ctr[0], sizeof(double)*3,
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

  // Texture objects
  //
  thrust::device_vector<cudaTextureObject_t> t_d = tex;

  // Space for Legendre coefficients 
  //
  thrust::device_vector<double> plm1_d((Lmax+1)*(Lmax+2)/2*Nthread);
  thrust::device_vector<double> plm2_d((Lmax+1)*(Lmax+2)/2*Nthread);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(double);

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

      double diff = 0.0, norm = 0.0;
      for (int k=0; k<3; k++) {
	double b  = cC->host_particles[i].acc[k];
	double a  = cC->Particles()[indx].acc[k];
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

      double diff = 0.0, norm = 0.0;
      for (int k=0; k<3; k++) {
	double b  = cC->host_particles[i].acc[k];
	double a  = cC->Particles()[indx].acc[k];
	diff += (a - b)*(a - b);
	norm += a*a;
      }
      std::cout << std::setw(14) << sqrt(diff/norm)
		<< std::setw(14) << sqrt(norm) << std::endl;
    }

  std::cout << std::string(16+14*8, '-') << std::endl;
  std::cout.precision(ss);
}
