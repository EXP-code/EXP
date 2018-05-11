#include <Component.H>
#include <SphericalBasis.H>

#define FSMALL 1.0e-16
#define MINEPS 1.0e-10

// Template structure to pass to kernel to avoid all of those pesky
// raw pointer casts
//
template <typename T>
struct dArray
{
  T*     _v;			// Pointer to underlying array data
  size_t _s;			// Number of elements
};
 
// Template function to convert device_vector to structure for passing
// from host to device in a cuda kernel
//
template <typename T>
dArray<T> toKernel(thrust::device_vector<T>& dev)
{
  dArray<T> dA;
  dA._v = thrust::raw_pointer_cast(&dev[0]);
  dA._s  = dev.size();
 
  return dA;
}

typedef std::pair<unsigned int, unsigned int> PII;

// Loop unroll
//
template <typename T, unsigned int blockSize>
__device__
void warpReduce(volatile T *sdata, unsigned int tid)
{
  if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
  if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
  if (blockSize >= 16) sdata[tid] += sdata[tid +  8];
  if (blockSize >=  8) sdata[tid] += sdata[tid +  4];
  if (blockSize >=  4) sdata[tid] += sdata[tid +  2];
  if (blockSize >=  2) sdata[tid] += sdata[tid +  1];
}

template <typename T, unsigned int blockSize>
__global__ void reduceSum(dArray<float> out, dArray<float> in,
			  unsigned int dim, unsigned int n)
{
  extern __shared__ T sdata[];

  unsigned int tid      = threadIdx.x;
  unsigned int gridSize = blockSize*gridDim.x*2;
    
  // printf("in reduceSum: b=%d/%d id=%d\n", blockIdx.x, gridDim.x, tid);

  for (unsigned j=0; j<dim; j++) {

    sdata[tid] = 0;
    
    unsigned int i = blockIdx.x*blockSize*2 + tid;

    while (i < n) {
      // Sanity check
      if (i+n*j>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
      if (i+blockSize<n)
	if (i+n*j+blockSize>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

      sdata[tid] +=
	in._v[i + n*j] + (i+blockSize<n ? in._v[i + blockSize + n*j] : T(0));
      i += gridSize;
    }
  
    __syncthreads();
  
    if (blockSize >= 1024) {
      if (tid < 512) {
	sdata[tid] += sdata[tid + 512];
	__syncthreads();
      }
    }
    
    if (blockSize >= 512) {
      if (tid < 256) {
	sdata[tid] += sdata[tid + 256];
      }
      __syncthreads();
    }

    if (blockSize >= 256) {
      if (tid < 128) {
	sdata[tid] += sdata[tid + 128];
      }
      __syncthreads();
    }

    if (blockSize >= 128) {
      if (tid < 64) {
	sdata[tid] += sdata[tid + 64];
      }
      __syncthreads();
    }

    if (tid < 32) {
      warpReduce<T, blockSize>(&sdata[0], tid);
    }

    if (tid == 0) {
      if (blockIdx.x + j*gridDim.x>=out._s)
	printf("reduceSum: out of bounds, b=%d/%d j=%d\n", blockIdx.x, gridDim.x, j);
     out._v[blockIdx.x + j*gridDim.x] = sdata[tid];
     /*
     printf("reduceSum: i=%d b=%d/%d j=%d\n",
	    blockIdx.x + j*gridDim.x, blockIdx.x, gridDim.x, j);
     */
    }
    __syncthreads();
  }
}

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

  if (ret >= (l+1)*(l+1)*nmax) {
    printf("Ilmn oab: %4d %4d %4d [%4d : %4d : %4d]\n", l, m, n, ret, (l+1)*(l+1)*nmax, nmax);
  }

  return ret;
}

__host__ __device__
void legendre_v(int lmax, float x, float* p)
{
  float fact, somx2, pll, pl1, pl2;
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
void legendre_v2(int lmax, float x, float* p, float* dp)
{
  float fact, somx2, pll, pl1, pl2;
  int m, l;

  p[0] = pll = 1.0f;
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

  if (1.0-fabs(x) < MINEPS) {
    if (x>0) x =   1.0 - MINEPS;
    else     x = -(1.0 - MINEPS);
  }

  somx2 = 1.0/(x*x - 1.0);
  dp[0] = 0.0;
  for (l=1; l<=lmax; l++) {
    for (m=0; m<l; m++)
      dp[Ilm(l, m)] = somx2*(x*l*p[Ilm(l, m)] - (l+m)*p[Ilm(l-1, m)]);
    dp[Ilm(l, l)] = somx2*x*l*p[Ilm(l, l)];
  }
}

__device__ __constant__ float cuRscale, cuXmin, cuXmax, cuDxi;
__device__ __constant__ int   cuNumr, cuCmap;

__global__
void testConstants()
{
  printf("** Rscale = %f\n", cuRscale);
  printf("** Xmin   = %f\n", cuXmin);
  printf("** Xmax   = %f\n", cuXmax);
  printf("** Dxi    = %f\n", cuDxi);
  printf("** Numr   = %d\n", cuNumr);
  printf("** Cmap   = %d\n", cuCmap);
}

__device__
float cu_r_to_xi(float r)
{
  float ret;

  if (cuCmap==1) {
    ret = (r/cuRscale-1.0)/(r/cuRscale+1.0);
  } else if (cuCmap==2) {
    ret = log(r);
  } else {
    ret = r;
  }    

  return ret;
}
    
__device__
float cu_xi_to_r(float xi)
{
  float ret;

  if (cuCmap==1) {
    ret = (1.0+xi)/(1.0 - xi) * cuRscale;
  } else if (cuCmap==2) {
    ret = exp(xi);
  } else {
    ret = xi;
  }

  return ret;
}

__device__
float cu_d_xi_to_r(float xi)
{
  float ret;

  if (cuCmap==1) {
    ret = 0.5*(1.0-xi)*(1.0-xi)/cuRscale;
  } else if (cuCmap==2) {
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

  cuda_safe_call(cudaMemcpyToSymbol(cuRscale, &f.scale, sizeof(float), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuRscale");

  cuda_safe_call(cudaMemcpyToSymbol(cuXmin, &f.xmin, sizeof(float), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuXmin");

  cuda_safe_call(cudaMemcpyToSymbol(cuXmax, &f.xmax, sizeof(float), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuXmax");

  cuda_safe_call(cudaMemcpyToSymbol(cuDxi, &f.dxi, sizeof(float), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuDxi");

  cuda_safe_call(cudaMemcpyToSymbol(cuNumr, &f.numr, sizeof(int), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuNumr");

  cuda_safe_call(cudaMemcpyToSymbol(cuCmap, &f.cmap, sizeof(int), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuCmap");
}

__global__
void testTexture(dArray<cudaTextureObject_t> tex, int nmax)
{
  printf("**DEVICE Texture compare\n");
  for (int l : {0, 1, 2}) {
    for (int j=0; j<10; j++) {
      int k = 1 + l*nmax;
      for (int i : {3980, 3990, 3995, 3999}) 
	printf("%5d %5d %5d %13.7e\n", l, j, i, tex1D<float>(tex._v[k+j], i));
    }
  }
}

__global__ void coordKernel
(dArray<cudaParticle> in, dArray<float> mass, dArray<float> Afac,
 dArray<float> phi, dArray<float> Plm, dArray<int> Indx, 
 unsigned int Lmax, unsigned int stride, PII lohi, float rmax)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;

  /*
    vector<double> ctr;
    if (mix) mix->getCenter(ctr);
  */
  float ctr[3] {0.0f, 0.0f, 0.0f};

  for (int n=0; n<stride; n++) {
    int i = tid*stride + n;
    int npart = i + lohi.first;

    if (npart < lohi.second) {

      if (npart>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

      cudaParticle p = in._v[npart];
    
      float xx = p.pos[0] - ctr[0];
      float yy = p.pos[1] - ctr[1];
      float zz = p.pos[2] - ctr[2];
      
      float r2 = (xx*xx + yy*yy + zz*zz);
      float r = sqrt(r2) + FSMALL;
      
      if (i>=mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

      mass._v[i] = -1.0;
      
      if (r<rmax) {
	
	mass._v[i] = p.mass;
	
	float costh = zz/r;
	phi._v[i] = atan2(yy,xx);
	
	if (i>=phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	float *plm = &Plm._v[psiz*i];
	legendre_v(Lmax, costh, plm);

	if (psiz*(i+1)>Plm._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	float x  = cu_r_to_xi(r);
	float xi = (x - cuXmin)/cuDxi;
	int indx = floor(xi);
	
	if (indx<0) indx = 0;
	if (indx>cuNumr-2) indx = cuNumr - 2;
	  
	
	Afac._v[i] = float(indx+1) - xi;
	if (Afac._v[i]<0.0 or Afac._v[i]>1.0)
	  printf("off grid: x=%f\n", xi);
	Indx._v[i] = indx;

	if (i>=Afac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=Indx._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
      }
    }
  }
}


__global__ void coefKernel
(dArray<float> coef, dArray<cudaTextureObject_t> tex,
 dArray<float> M, dArray<float> A, dArray<float> P,
 dArray<float> L, dArray<int> I,  int stride, 
 int l, int m, unsigned Lmax, unsigned int nmax, PII lohi)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;
  const unsigned int N = lohi.second - lohi.first;

  float fac0 = 4.0*M_PI;

  for (int istr=0; istr<stride; istr++) {

    int i = tid*stride + istr;

    if (i<N) {

      float mass = M._v[i];

      if (i>=M._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

      if (mass>0.0) {

	if (i>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	float phi  = P._v[i];
	float cosp = cos(phi*m);
	float sinp = sin(phi*m);
	
	if (psiz*(i+1)>L._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	
	float *plm = &L._v[psiz*i];
	
	for (int n=0; n<nmax; n++) {
	  // Do the interpolation
	  //
	  float a = A._v[i];
	  float b = 1.0 - a;
	  int ind = I._v[i];
	  
	  if (i>=A._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	  if (i>=I._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	  float p0 =
	    a*tex1D<float>(tex._v[0], ind  ) +
	    b*tex1D<float>(tex._v[0], ind+1) ;


	  int k = 1 + l*nmax + n;

	  if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	  float v = (
		     a*tex1D<float>(tex._v[k], ind  ) +
		     b*tex1D<float>(tex._v[k], ind+1)
		     ) * p0 * plm[Ilm(l, m)] * M._v[i] * fac0;
	  
	  
	  coef._v[(2*n+0)*N + i] = v * cosp;
	  coef._v[(2*n+1)*N + i] = v * sinp;

	  if ((2*n+0)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	  if ((2*n+1)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	}
      }
    }
  }
}

__global__ void
forceKernel(dArray<cudaParticle> in, dArray<float> coef,
	    dArray<cudaTextureObject_t> tex, dArray<float> L1, dArray<float> L2,
	    int stride, unsigned Lmax, unsigned int nmax, PII lohi, float rmax)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;
  // const unsigned int N = lohi.second - lohi.first;

  // const float fac0 = 4.0*M_PI;
  // const float dfac = 0.25/M_PI;

  /*
    vector<double> ctr;
    if (mix) mix->getCenter(ctr);
  */
  float ctr[3] {0.0f, 0.0f, 0.0f};

  for (int n=0; n<stride; n++) {
    int i = tid*stride + n;
    int npart = i + lohi.first;

    if (npart < lohi.second) {
      
      if (npart>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

      cudaParticle p = in._v[npart];
      
      float xx = p.pos[0] - ctr[0];
      float yy = p.pos[1] - ctr[1];
      float zz = p.pos[2] - ctr[2];
      
      float r2 = (xx*xx + yy*yy + zz*zz);
      float r  = sqrt(r2) + FSMALL;
      
      int ioff = 0;
      // float rs = 0.0;
      float r0;

      if (r>rmax) {
	ioff = 1;
	r0 = r;
	r = rmax;
	// rs = r/cuRscale;
      }

      float costh = zz/r;
      float phi   = atan2(yy,xx);
      float RR    = xx*xx + yy*yy;
      
      float *plm  = &L1._v[psiz*i];
      float *plm2 = &L2._v[psiz*i];
      legendre_v2(Lmax, costh, plm, plm2);

      float  x = cu_r_to_xi(r);
      float xi = (x - cuXmin)/cuDxi;
      int  ind = floor(xi);
      
      if (ind<1) ind = 1;
      if (ind>cuNumr-2) ind = cuNumr - 2;
      
      float a = float(ind+1) - xi;
      if (a<0.0 or a>1.0) printf("off grid: x=%f\n", xi);
      float b = 1.0 - a;
      
      // Do the interpolation
      //
      float pm1 = tex1D<float>(tex._v[0], ind-1);
      float p00 = tex1D<float>(tex._v[0], ind  );
      float pp1 = tex1D<float>(tex._v[0], ind+1);

      // For force accumulation
      //
      float potl = 0.0;
      float potr = 0.0;
      float pott = 0.0;
      float potp = 0.0;

      // l loop
      //
      for (int l=0; l<Lmax; l++) {

	float fac1 = (2.0*l + 1.0)/(4.0*M_PI);

	// m loop
	//
	for (int m=0; m<=l; m++) {

	  int pindx = Ilm(l, m);

	  float Plm  =  plm[pindx];
	  float Plm2 = plm2[pindx];
      
	  float pp_c = 0.0;
	  float dp_c = 0.0;
	  float pp_s = 0.0;
	  float dp_s = 0.0;
	  
	  int indxC = Ilmn(l, m, 0, 'c', nmax);
	  int indxS = Ilmn(l, m, 0, 's', nmax);

	  float cosp = cos(phi*m);
	  float sinp = sin(phi*m);

	  for (size_t n=0; n<nmax; n++) {
	
	    int k = 1 + l*nmax + n;
	
	    if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	
	    float um1 = tex1D<float>(tex._v[k], ind-1);
	    float u00 = tex1D<float>(tex._v[k], ind  );
	    float up1 = tex1D<float>(tex._v[k], ind+1);
	    
	    float v = a*u00 + b*up1;
	    
	    float dv =
	      ( (b - 0.5)*um1*pm1 - 2.0*b*u00*p00 + (b + 0.5)*up1*pp1 );
	    
	    pp_c +=  v * coef._v[indxC+n];
	    dp_c += dv * coef._v[indxC+n];
	    if (m>0) {
	      pp_s +=  v * coef._v[indxS+n];
	      dp_s += dv * coef._v[indxS+n];
	    }

	  } // END: n loop
	  
	  if (m==0) {

	    float fac2 = fac1 * Plm;
	
	    if (ioff) {
	      pp_c *= pow(rmax/r0, (float)(l+1));
	      dp_c  = -pp_c/r0 * (float)(l+1);
	    }
	    
	    potl += fac2 * pp_c;
	    potr += fac2 * dp_c;
	    pott += fac1 * pp_c * plm2[Ilm(l, m)];
	    potp += 0.0;
	    
	  } else {

	    float cosm = cos(phi*m);
	    float sinm = sin(phi*m);
	    
	    if (ioff) {
	      float facp = pow(rmax/r0,(double)(l+1));
	      float facdp = -1.0/r0 * (l+1);
	      pp_c *= facp;
	      pp_s *= facp;
	      dp_c = pp_c * facdp;
	      dp_s = pp_s * facdp;
	    }

	    float numf = 1.0, denf = 1.0;
	    for (int i=1; i<=l-m; i++) numf *= i;
	    for (int i=1; i<=l+m; i++) denf *= i;
	    
	    float fac2 = 2.0 * numf/denf * fac1;
	    
	    potl += fac2 * Plm  * ( pp_c*cosm + pp_s*sinm);
	    potr += fac2 * Plm  * ( dp_c*cosm + dp_s*sinm);
	    pott += fac2 * Plm2 * ( pp_c*cosm + pp_s*sinm);
	    potp += fac2 * Plm  * (-pp_c*sinm + pp_s*cosm)*m;
	  }

	} // END: m loop

      } // END: l loop

      in._v[npart].acc[0] = -(potr*xx/r - pott*xx*zz/(r*r*r));
      in._v[npart].acc[1] = -(potr*yy/r - pott*yy*zz/(r*r*r));
      in._v[npart].acc[2] = -(potr*zz/r - pott*RR/(r2*r));
      in._v[npart].pot    = potl;

    } // Particle index block

  } // END: stride loop

}

void SphericalBasis::determine_coefficients_cuda(const Matrix& expcoef)
{
  std::cout << std::scientific;

  // Sort particles and get coefficient size
  //
  PII lohi = cC->CudaSortByLevel(mlevel, multistep);

  unsigned int N         = lohi.second - lohi.first;
  unsigned int stride    = 1;
  unsigned int gridSize  = N/BLOCK_SIZE/stride;

  /*
  if (gridSize>128) {
    stride = N/BLOCK_SIZE/128 + 1;
    gridSize = N/BLOCK_SIZE/stride;
  }
  */

  if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

  std::cout << "**" << std::endl
	    << "** N      = " << N          << std::endl
	    << "** Stride = " << stride     << std::endl
	    << "** Block  = " << BLOCK_SIZE << std::endl
	    << "** Grid   = " << gridSize   << std::endl
	    << "**" << std::endl;

  // Create space for coefficient reduction
  //
  thrust::device_vector<float> dN_coef(2*nmax*N);
  thrust::device_vector<float> dc_coef(2*nmax*gridSize);

  // Texture objects
  //
  thrust::device_vector<cudaTextureObject_t> t_d = tex;

  // Space for per particle values
  //
  thrust::device_vector<float> work_d(gridSize*BLOCK_SIZE);

  // Space for Legendre coefficients 
  //
  thrust::device_vector<float> plm_d((Lmax+1)*(Lmax+2)/2*N);
  thrust::device_vector<float> r_d(N), m_d(N), a_d(N), p_d(N);
  thrust::device_vector<int>   i_d(N);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(float);

  // For debugging
  //
  if (false) {
    testConstants<<<1, 1>>>();
    
    static bool firstime = true;
    testTexture<<<1, 1>>>(toKernel(t_d), nmax);
    firstime == false;
  }

  std::vector<float> coefs((Lmax+1)*(Lmax+1)*nmax);

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
  for (int l=0; l<=Lmax; l++) {
    for (int m=0; m<=l; m++) {
      coefKernel<<<gridSize, BLOCK_SIZE>>>
	(toKernel(dN_coef), toKernel(t_d), toKernel(m_d),
	 toKernel(a_d), toKernel(p_d), toKernel(plm_d), toKernel(i_d),
	 stride, l, m, Lmax, nmax, lohi);

      cudaDeviceSynchronize();

				// Begin the reduction per grid block
      int osize = nmax*2;	// 
      reduceSum<float, BLOCK_SIZE><<<gridSize, BLOCK_SIZE, sMemSize>>>
	(toKernel(dc_coef), toKernel(dN_coef), osize, N);
      
      cudaDeviceSynchronize();
				// Finish the reduction
				//
      for (size_t j=0; j<nmax; j++) {
	coefs[Ilmn(l, m, 'c', j, nmax)] =
	  thrust::reduce(dc_coef.begin() + gridSize*(2*j+0),
			 dc_coef.begin() + gridSize*(2*j+1));
	if (m>0)
	  coefs[Ilmn(l, m, 's', j, nmax)] =
	    thrust::reduce(dc_coef.begin() + gridSize*(2*j+1),
			   dc_coef.begin() + gridSize*(2*j+2));
      }
    }
  }

  // DEBUG
  //
  if (false) {
    std::cout << "L=M=0 coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Ilmn(0, 0, 'c', n, nmax)]
		<< std::setw(16) << expcoef[0][n+1]
		<< std::endl;
    }

    std::cout << "L=1, M=0 coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Ilmn(1, 0, 'c', n, nmax)]
		<< std::setw(16) << expcoef[1][n+1]
		<< std::endl;
    }

    std::cout << "L=1, M=1c coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Ilmn(1, 1, 'c', n, nmax)]
		<< std::setw(16) << expcoef[2][n+1]
		<< std::endl;
    }

    std::cout << "L=1, M=1s coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Ilmn(1, 1, 's', n, nmax)]
		<< std::setw(16) << expcoef[3][n+1]
		<< std::endl;
    }
  }

  //
  // TEST comparison of coefficients for debugging
  //
  if (false) {

    struct Element
    {
      double d;
      float  f;
      
      int  l;
      int  m;
      int  n;
      
      char cs;
    }
    elem;

    std::map<double, Element> compare;

    std::ofstream out("test.dat");

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
	    elem.f = coefs[Ilmn(l, m, 'c', n-1, nmax)];
	    
	    double test = fabs(elem.d - elem.f);
	    if (fabs(elem.d)>1.0e-4) test /= fabs(elem.d);
	    
	    compare[test] = elem;
	    
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
	    elem.f = coefs[Ilmn(l, m, 'c', n-1, nmax)];

	    out << std::setw( 5) << l
		<< std::setw( 5) << m
		<< std::setw( 5) << n
		<< std::setw( 5) << 'c'
		<< std::setw( 5) << Ilmn(l, m, 'c', n-1, nmax)
		<< std::setw(14) << elem.d
		<< std::setw(14) << elem.f
		<< std::endl;

	    double test = fabs(elem.d - elem.f);
	    if (fabs(elem.d)>1.0e-4) test /= fabs(elem.d);

	    compare[test] = elem;
	  }
	  for (int n=1; n<=nmax; n++) {
	    elem.l = l;
	    elem.m = m;
	    elem.n = n;
	    elem.cs = 's';
	    elem.d = expcoef[loffset+moffset+1][n];
	    elem.f = coefs[Ilmn(l, m, 's', n-1, nmax)];

	    out << std::setw( 5) << l
		<< std::setw( 5) << m
		<< std::setw( 5) << n
		<< std::setw( 5) << 's'
		<< std::setw( 5) << Ilmn(l, m, 's', n-1, nmax)
		<< std::setw(14) << elem.d
		<< std::setw(14) << elem.f
		<< std::endl;
	    
	    double test = fabs(elem.d - elem.f);
	    if (fabs(elem.d)>1.0e-4) test /= fabs(elem.d);
	    
	    compare[test] = elem;
	  }
	  moffset+=2;
	}
      }
    }
    
    std::map<double, Element>::iterator best = compare.begin();
    std::map<double, Element>::iterator midl = best;
    std::advance(midl, compare.size()/2);
    std::map<double, Element>::reverse_iterator last = compare.rbegin();
    
    std::cout << "Best case: ["
	      << std::setw( 2) << best->second.l << ", "
	      << std::setw( 2) << best->second.m << ", "
	      << std::setw( 2) << best->second.n << ", "
	      << std::setw( 2) << best->second.cs << "] = "
	      << std::setw(15) << best->second.d
	      << std::setw(15) << best->second.f
	      << std::setw(15) << fabs(best->second.d - best->second.f)
	      << std::endl;
  
    std::cout << "Mid case:  ["
	      << std::setw( 2) << midl->second.l << ", "
	      << std::setw( 2) << midl->second.m << ", "
	      << std::setw( 2) << midl->second.n << ", "
	      << std::setw( 2) << midl->second.cs << "] = "
	      << std::setw(15) << midl->second.d
	      << std::setw(15) << midl->second.f
	      << std::setw(15) << fabs(midl->second.d - midl->second.f)
	      << std::endl;
    
    std::cout << "Last case: ["
	      << std::setw( 2) << last->second.l << ", "
	      << std::setw( 2) << last->second.m << ", "
	      << std::setw( 2) << last->second.n << ", "
	      << std::setw( 2) << last->second.cs << "] = "
	      << std::setw(15) << last->second.d
	      << std::setw(15) << last->second.f
	      << std::setw(15) << fabs(last->second.d - last->second.f)
	      << std::endl;
  }

}


void SphericalBasis::determine_acceleration_cuda()
{
  std::cout << std::scientific;

  // Sort particles and do all particles at or above mlevel
  //
  PII lohi = cC->CudaSortByLevel(mlevel, multistep);

  unsigned int N         = lohi.second - lohi.first;
  unsigned int stride    = 1;
  unsigned int gridSize  = N/BLOCK_SIZE/stride;

  /*
  if (gridSize>128) {
    stride = N/BLOCK_SIZE/128 + 1;
    gridSize = N/BLOCK_SIZE/stride;
  }
  */

  if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

  std::cout << "**" << std::endl
	    << "** N      = " << N          << std::endl
	    << "** Stride = " << stride     << std::endl
	    << "** Block  = " << BLOCK_SIZE << std::endl
	    << "** Grid   = " << gridSize   << std::endl
	    << "**" << std::endl;

  // Texture objects
  //
  thrust::device_vector<cudaTextureObject_t> t_d = tex;

  // Space for Legendre coefficients 
  //
  thrust::device_vector<float> plm1_d((Lmax+1)*(Lmax+2)/2*N);
  thrust::device_vector<float> plm2_d((Lmax+1)*(Lmax+2)/2*N);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(float);

  // For debugging
  //
  std::vector<float> coefs((Lmax+1)*(Lmax+1)*nmax);

  // Do the work
  //
  forceKernel<<<gridSize, BLOCK_SIZE, sMemSize>>>
    (toKernel(cC->cuda_particles), toKernel(dev_coefs), toKernel(t_d),
     toKernel(plm1_d), toKernel(plm2_d), stride, Lmax, nmax, lohi, rmax);
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
  host_coefs = dev_coefs;

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

