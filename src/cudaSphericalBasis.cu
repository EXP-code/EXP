#include <Component.H>
#include <SphericalBasis.H>

#define FSMALL 1.0e-16


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
__device__ void reduceSum(T *out, T *work, T *sdata, unsigned int n)
{
  unsigned int tid      = threadIdx.x;
  unsigned int gridSize = (blockSize*2)*gridDim.x;
    
  unsigned int i = blockIdx.x*(blockSize*2) + tid;

  sdata[tid] = 0.0;

  while (i < n) {
    sdata[tid] += work[i] + (i+blockSize<n ? work[i+blockSize] : T(0));
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

  if (tid == 0 and blockIdx.x < gridDim.x/2) {
    *out = sdata[tid];
  }

  __syncthreads();
}


__host__ __device__
int Ilm(int l, int m)
{
  if (l==0) return 0;
  return (l+1)*(l+2)/2 + m;
}

__host__ __device__
int Ilmn(int l, int m, char cs, int n, int nmax)
{
  if (l==0) return n;
  if (m==0) return l*l*nmax + n;
  return (l*l + 1 + 2*(m-1) + (cs=='s' ? 1 : 0))*nmax + n;
}

__host__ __device__
void legendre_array(int lmax, float x, float* p)
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
    pl2 = p[m*(lmax+1) + m];
    p[Ilm(m+1, m)] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[Ilm(l, m)] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }
}

extern __constant__ float cuRscale, cuXmin, cuXmax, cuDxi;
extern __constant__ int   cuNumr, cuCmap;

__device__
float cu_r_to_xi(float r)
{
  float ret;

  printf("cuRscale=%f\n", cuRscale);

  if (cuCmap==1) {
    ret =  (r/cuRscale-1.0)/(r/cuRscale+1.0);
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
    ret =(1.0+xi)/(1.0 - xi) * cuRscale;
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
  
  cuda_mapping_constants f = get_cuda_mapping_constants();

  cuda_safe_call(cudaMemcpyToSymbol(cuRscale, &f.scale, sizeof(float), size_t(0), cudaMemcpyHostToDevice), "Error copying cuRscale");

  cuda_safe_call(cudaMemcpyToSymbol(cuXmin, &f.min, sizeof(float), size_t(0), cudaMemcpyHostToDevice), "Error copying cuXmin");

  cuda_safe_call(cudaMemcpyToSymbol(cuXmax, &f.xmax, sizeof(float), size_t(0), cudaMemcpyHostToDevice), "Error copying cuXmax");

  cuda_safe_call(cudaMemcpyToSymbol(cuDxi, &f.dxi, sizeof(float), size_t(0), cudaMemcpyHostToDevice), "Error copying cuDxi");

  cuda_safe_call(cudaMemcpyToSymbol(cuNumr, &f.numr, sizeof(int), size_t(0), cudaMemcpyHostToDevice), "Error copying cuNumr");

  cuda_safe_call(cudaMemcpyToSymbol(cuCmap, &f.cmap, sizeof(int), size_t(0), cudaMemcpyHostToDevice), "Error copying cuCmap");
}


__global__ void coefficient_kernel
(float *coef, float*work, cudaTextureObject_t *tex, float *plm, cudaParticle* in,
 int Lmax, int nmax, unsigned int nbeg, unsigned int nend, float rmax)
{
  extern __shared__ float sdata[];
		
  const int        tid = blockDim.x * blockIdx.x + threadIdx.x;

  float cosp0, sinp0, cosp1, sinp1, cosp2, sinp2, cosp, sinp;
  float fac0=4.0*M_PI;

  float adb = 1.0; // = component->Adiabatic();

  /*
    vector<double> ctr;
    if (mix) mix->getCenter(ctr);
  */
  float ctr[3] {0.0f, 0.0f, 0.0f};

  cudaParticle p = in[tid + nbeg];
    
  float mass =  p.mass * adb;

  // Adjust mass for subset
  // if (subset) mass /= ssfrac;
    
  float xx = p.pos[0] - ctr[0];
  float yy = p.pos[1] - ctr[1];
  float zz = p.pos[2] - ctr[2];


  float r2 = (xx*xx + yy*yy + zz*zz);
  float r = sqrt(r2) + FSMALL;
  
  if (r<rmax) {
    
    float costh = zz/r;
    float phi = atan2(yy,xx);
    
    cosp0 = cos(phi);
    sinp0 = sin(phi);
    
    legendre_array(Lmax, costh, plm);
    
    float x  = cu_r_to_xi(r);
    float xi = (x - cuXmin)/cuDxi;
    int indx = floor(xi);
    
    if (indx<0) indx = 0;
    if (indx>cuNumr-2) indx = cuNumr - 2;
    
    float a = float(indx+1) - xi;
    float b = 1.0 - a;
    
    //		l loop
    for (int l=0; l<=Lmax; l++) {
      //		m loop
      for (int m=0; m<=l; m++) {
	
	if (m==0) {
	  
	  cosp = 1.0;
	  sinp = 0.0;
	  
	  for (int n=0; n<nmax; n++) {
	    // Do the interpolation
	    
	    int k = l*nmax + n;
	    
	    float v =
	      a*tex1D<float>(tex[k], indx  ) +
	      b*tex1D<float>(tex[k], indx+1) ;
	    
	    work[n-nbeg] = v * plm[Ilmn(l, m, 'c', n, nmax)] * mass * fac0;

	    __syncthreads();
	    
	    reduceSum<float, BLOCK_SIZE>(&coef[Ilmn(l, m, 'c', n, nmax)], work, sdata, nend - nbeg);
	  }
	  
	  cosp1 = cosp;
	  sinp1 = sinp;
	  cosp  = cosp0;
	  sinp  = sinp0;
	}
	else {
	  for (int n=1; n<=nmax; n++) {
	    
	    int k = l*nmax + n;

	    float v =
	      a*tex1D<float>(tex[k], indx  ) +
	      b*tex1D<float>(tex[k], indx+1) ;
	    
	    v *= plm[Ilm(l, m)] * mass * fac0;
	    
	    work[n-nbeg] = v * cosp;
	    
	    __syncthreads();
	    
	    reduceSum<float, BLOCK_SIZE>(&coef[Ilmn(l, m, 'c', n, nmax)], work, sdata, nend - nbeg);
	    
	    work[n-nbeg] = v * sinp;
	    
	    reduceSum<float, BLOCK_SIZE>(&coef[Ilmn(l, m, 's', n, nmax)], work, sdata, nend - nbeg);
	    
	    __syncthreads();
	    
	  }
	  
	  cosp2 = cosp1;
	  sinp2 = sinp1;
	  cosp1 = cosp;
	  sinp1 = sinp;
	  
	  cosp  = cosp0*cosp1 - cosp2;
	  sinp  = sinp0*sinp1 - sinp2;
	}
      }
    }
  }

}


void SphericalBasis::determine_coefficients_cuda(const Matrix& expcoef)
{
  // Create space for coefficients
  size_t csize = (Lmax+1)*(Lmax+1)*nmax;
  thrust::host_vector<float> h_coef(csize, 0.0f);
  thrust::device_vector<float> d_coef = h_coef;

  // Sort particles and get coefficient size
  std::pair<unsigned, unsigned> lohi = cC->CudaSortByLevel(mlevel, multistep);

  unsigned int N         = lohi.second - lohi.first;
  unsigned int gridSize  = N/BLOCK_SIZE + (N % BLOCK_SIZE > 0 ? 1 : 0);

  thrust::device_vector<cudaTextureObject_t> t_d = tex;
  thrust::device_vector<float> work_d((Lmax+1)*(Lmax+1));
  thrust::device_vector<float> p_d((Lmax+1)*(Lmax+2)/2);


  float               *C = thrust::raw_pointer_cast(&d_coef[0]);
  float               *W = thrust::raw_pointer_cast(&work_d[0]);
  cudaTextureObject_t *T = thrust::raw_pointer_cast(&t_d[0]);
  float               *P = thrust::raw_pointer_cast(&p_d[0]);
  cudaParticle        *I = thrust::raw_pointer_cast(&cC->cuda_particles[0]);

  int sMemSize = csize * sizeof(float);

  coefficient_kernel<<<gridSize, BLOCK_SIZE, sMemSize>>>
    (C, W, T, P, I, Lmax, nmax, lohi.first, lohi.second, rmax);
  
  // Copy from device to host
  //
  h_coef = d_coef;
  
  //
  // TEST comparison of coefficients
  //
  double dmax = 0.0;
  int  l0 = 0;
  int  m0 = 0;
  int  n0 = 0;
  char cs = 'c';

  // int loffset, moffset;

  //		l loop
  for (int l=0, loffset=0; l<=Lmax; loffset+=(2*l+1), l++) {
    //		m loop
    for (int m=0, moffset=0; m<=l; m++) {

      if (m==0) {
	for (int n=1; n<=nmax; n++) {
	  double test = fabs(expcoef[loffset+moffset][n]  - h_coef[Ilmn(l, m, 'c', n-1, nmax)]);
	  if (test > dmax) {
	    dmax = test;
	    l0 = l;
	    m0 = m;
	    n0 = n;
	    cs = 'c';
	  }
	}

	moffset++;
      }
      else {
	for (int n=1; n<=nmax; n++) {

	  double test = fabs(expcoef[loffset+moffset][n]  - h_coef[Ilmn(l, m, 'c', n-1, nmax)]);
	  if (test > dmax) {
	    dmax = test;
	    l0 = l;
	    m0 = m;
	    n0 = n;
	    cs = 'c';
	  }

	  test = fabs(expcoef[loffset+moffset+1][n]  - h_coef[Ilmn(l, m, 's', n-1, nmax)]);
	  if (test > dmax) {
	    dmax = test;
	    l0 = l;
	    m0 = m;
	    n0 = n;
	    cs = 's';
	  }
	}
	moffset+=2;
      }
    }
  }

  std::cout << "CUDA coefficient comparision" << std::endl
	    << std::string(28, '-')           << std::endl
	    << "Worst case: [" << l0 << ", " << m0 << ", " << n0 << ", " << cs << "] = " << dmax << std::endl
	    << std::string(28, '-')           << std::endl;

}

