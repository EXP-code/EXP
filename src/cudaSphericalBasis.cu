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
    // printf("Ans=%e\n", *out);
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
    pl2 = p[Ilm(m, m)];
    p[Ilm(m+1, m)] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[Ilm(l, m)] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
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

  cuda_safe_call(cudaMemcpyToSymbol(cuRscale, &f.scale, sizeof(float), size_t(0), cudaMemcpyHostToDevice), "Error copying cuRscale");

  cuda_safe_call(cudaMemcpyToSymbol(cuXmin, &f.xmin, sizeof(float), size_t(0), cudaMemcpyHostToDevice), "Error copying cuXmin");

  cuda_safe_call(cudaMemcpyToSymbol(cuXmax, &f.xmax, sizeof(float), size_t(0), cudaMemcpyHostToDevice), "Error copying cuXmax");

  cuda_safe_call(cudaMemcpyToSymbol(cuDxi, &f.dxi, sizeof(float), size_t(0), cudaMemcpyHostToDevice), "Error copying cuDxi");

  cuda_safe_call(cudaMemcpyToSymbol(cuNumr, &f.numr, sizeof(int), size_t(0), cudaMemcpyHostToDevice), "Error copying cuNumr");

  cuda_safe_call(cudaMemcpyToSymbol(cuCmap, &f.cmap, sizeof(int), size_t(0), cudaMemcpyHostToDevice), "Error copying cuCmap");
}


__global__ void coefKernel
(float *coef, float*work, cudaTextureObject_t *tex, float *plm, cudaParticle* in,
 int Lmax, int nmax, unsigned int nbeg, unsigned int nend, float rmax)
{
  extern __shared__ float sdata[];
		
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  bool exist = true;
  bool rgood = true;

  if (tid+nbeg >= nend) {
    exist = rgood = false;
  }

  float adb = 1.0; // = component->Adiabatic();

  float cosp0, sinp0, cosp1, sinp1, cosp2, sinp2, cosp, sinp, v, a, b, xi;
  float fac0 = 4.0*M_PI, mass = 0.0;
  int indx;

  /*
    vector<double> ctr;
    if (mix) mix->getCenter(ctr);
  */
  float ctr[3] {0.0f, 0.0f, 0.0f};

  if (exist) {

    cudaParticle p = in[tid + nbeg];
    
    mass =  p.mass * adb;

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
    
      // float x, xi;
      float x;

      x  = cu_r_to_xi(r);
      xi = (x - cuXmin)/cuDxi;
      indx = floor(xi);
      
      if (indx<0) indx = 0;
      if (indx>cuNumr-2) indx = cuNumr - 2;
      
      a = float(indx+1) - xi;
      b = 1.0 - a;
    } else {
      rgood = false;
    }
  }

    
  //		l loop
  for (int l=0; l<=Lmax; l++) {
    //          m loop
    for (int m=0; m<=l; m++) {
      
      if (m==0) {
	
	cosp = 1.0;
	sinp = 0.0;
	
	for (int n=0; n<nmax; n++) {
	  if (rgood) {
	    // Do the interpolation
	    //
	    int k = l*nmax + n;
	    
	    v =
	      a*tex1D<float>(tex[k], indx  ) +
	      b*tex1D<float>(tex[k], indx+1) ;
	      
	    printf("[%3d %3d %3d] %f %f\n", l, m, n, plm[Ilm(l, m)], v);

	    work[tid] = v * plm[Ilm(l, m)] * mass * fac0;

	  } else {
	    // Ignore value
	    //
	    if (exist) {
	      work[tid] = 0.0; // out of bounds
	      if (tid >= nend-nbeg) {
		printf("Error: %d %d %d\n", tid, nbeg, nend);
	      }
	    }
	  }

	  __syncthreads();
	    
	  float z;
	  reduceSum<float, BLOCK_SIZE>(&z, work, sdata, nend - nbeg);
	  coef[Ilmn(l, m, 'c', n, nmax)] = z;
	  // printf("[%d, %d] z=%e\n", l, m, z);
	}
	
	cosp1 = cosp;
	sinp1 = sinp;
	cosp  = cosp0;
	sinp  = sinp0;
      }
      else {
	for (int n=0; n<nmax; n++) {
	  
	  if (rgood) {
	    // Do the interpolation
	    //
	    int k = l*nmax + n;
	    
	    v =
	      a*tex1D<float>(tex[k], indx  ) +
	      b*tex1D<float>(tex[k], indx+1) ;
	    
	    v *= plm[Ilm(l, m)] * mass * fac0;
	    
	    printf("[%3d %3d %3d] %f %f\n", l, m, n, plm[Ilm(l, m)], v);

	    work[tid] = v * cosp;
	  } else {
	    // Ignore value
	    //
	    if (exist) {
	      work[tid] = 0.0;
	      if (tid >= nend-nbeg) {
		printf("Error: %d %d %d\n", tid, nbeg, nend);
	      }
	    }
	  }
	    
	  __syncthreads();
	  
	  float z;
	  reduceSum<float, BLOCK_SIZE>(&z, work, sdata, nend - nbeg);
	  coef[Ilmn(l, m, 'c', n, nmax)] = z;
	  // printf("[%d, %d] z=%e\n", l, m, z);
	      
	  if (rgood)
	    work[tid] = v * sinp;
	  else {
	    if (exist) {
	      work[tid] = 0.0;
	      if (tid >= nend-nbeg) {
		printf("Error: %d %d %d\n", tid, nbeg, nend);
	      }
	    }
	  }
	  
	  reduceSum<float, BLOCK_SIZE>(&z, work, sdata, nend - nbeg);
	  coef[Ilmn(l, m, 's', n, nmax)] = z;
	  // printf("[%d, %d] z=%e\n", l, m, z);
	    
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

  assert(gridSize*BLOCK_SIZE >= N);

  thrust::device_vector<cudaTextureObject_t> t_d = tex;
  thrust::device_vector<float> work_d(lohi.second - lohi.first+1);
  thrust::device_vector<float> p_d((Lmax+1)*(Lmax+2)/2);


  float               *C = thrust::raw_pointer_cast(&d_coef[0]);
  float               *W = thrust::raw_pointer_cast(&work_d[0]);
  cudaTextureObject_t *T = thrust::raw_pointer_cast(&t_d[0]);
  float               *P = thrust::raw_pointer_cast(&p_d[0]);
  cudaParticle        *I = thrust::raw_pointer_cast(&cC->cuda_particles[0]);

  int sMemSize = csize * sizeof(float);

  testConstants<<<1, 1>>>();

  coefKernel<<<gridSize, BLOCK_SIZE, sMemSize>>>
    (C, W, T, P, I, Lmax, nmax, lohi.first, lohi.second, rmax);
  
  // Copy from device to host
  //
  h_coef = d_coef;
  
  //
  // TEST comparison of coefficients
  //

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
	  elem.f = h_coef[Ilmn(l, m, 'c', n-1, nmax)];

	  compare[fabs(elem.d - elem.f)] = elem;
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
	  elem.f = h_coef[Ilmn(l, m, 'c', n-1, nmax)];

	  compare[fabs(elem.d - elem.f)] = elem;

	  elem.l = l;
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 's';
	  elem.d = expcoef[loffset+moffset+1][n];
	  elem.f = h_coef[Ilmn(l, m, 's', n-1, nmax)];

	  compare[fabs(elem.d - elem.f)] = elem;
	}
	moffset+=2;
      }
    }
  }

  std::map<double, Element>::iterator best = compare.begin();
  std::map<double, Element>::iterator midl = best;
  std::advance(midl, compare.size()/2);
  std::map<double, Element>::reverse_iterator last = compare.rbegin();

  std::cout << std::scientific;

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

