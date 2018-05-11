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
__global__ void reduceSum(float *out, float *in,
			  unsigned int dim, unsigned int n)
{
  extern __shared__ T sdata[];

  unsigned int tid      = threadIdx.x;
  unsigned int gridSize = blockSize*gridDim.x*2;
    
  for (unsigned j=0; j<dim; j++) {

    sdata[tid] = 0;
    
    unsigned int i = blockIdx.x*blockSize*2 + tid;

    while (i < n) {
      sdata[tid] +=
	in[i + n*j] + (i+blockSize<n ? in[i + blockSize + n*j] : T(0));
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
      out[blockIdx.x + j*gridDim.x] = sdata[tid];
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

__global__
void testTexture(cudaTextureObject_t *tex, int nmax)
{
  printf("**DEVICE Texture compare\n");
  for (int l : {0, 1, 2}) {
    for (int j=0; j<10; j++) {
      int k = 1 + l*nmax;
      for (int i : {3980, 3990, 3995, 3999}) 
	printf("%5d %5d %5d %13.7e\n", l, j, i, tex1D<float>(tex[k+j], i));
    }
  }
}

__global__ void coordKernel
(cudaParticle* in, float *M, float *A, float *P, float *L, int *I,
 unsigned int Lmax, unsigned int stride, unsigned int nbeg, unsigned int nend, float rmax)
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
    int npart = i + nbeg;

    if (npart < nend) {

      cudaParticle p = in[npart];
    
      float xx = p.pos[0] - ctr[0];
      float yy = p.pos[1] - ctr[1];
      float zz = p.pos[2] - ctr[2];
      
      float r2 = (xx*xx + yy*yy + zz*zz);
      float r = sqrt(r2) + FSMALL;
      
      M[i] = -1.0;
      
      if (r<rmax) {
	
	M[i] = p.mass;
	
	float costh = zz/r;
	P[i] = atan2(yy,xx);
	
	float *plm = &L[psiz*i];
	legendre_array(Lmax, costh, plm);
	
	float x  = cu_r_to_xi(r);
	float xi = (x - cuXmin)/cuDxi;
	int indx = floor(xi);
	
	if (indx<0) indx = 0;
	if (indx>cuNumr-2) indx = cuNumr - 2;
	
	A[i] = float(indx+1) - xi;
	if (A[i]<0.0 or A[i]>1.0) printf("off grid: x=%f\n", xi);
	I[i] = indx;
      }
    }
  }
}


__global__ void coefKernel
(float *coef, cudaTextureObject_t *tex,
 float *M, float *A, float *P, float *L, int *I,  int stride, 
 int l, int m, unsigned Lmax, unsigned int nmax, unsigned int nbeg, unsigned int nend)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const int psiz  = (Lmax+1)*(Lmax+2)/2;

  float fac0 = 4.0*M_PI;

  for (int istr=0; istr<stride; istr++) {

    int i = tid*stride + istr;

    if (i<nend - nbeg) {

      float mass = M[i];

      if (mass>0.0) {

	float phi  = P[i];
	float cosp = cos(phi*m);
	float sinp = sin(phi*m);
	
	float *plm = &L[psiz*i];
	
	for (int n=0; n<nmax; n++) {
	  // Do the interpolation
	  //
	  float a = A[i];
	  float b = 1.0 - a;
	  int ind = I[i];
	  
	  float p0 =
	    a*tex1D<float>(tex[0], ind  ) +
	    b*tex1D<float>(tex[0], ind+1) ;


	  int k = 1 + l*nmax + n;

	  float v = (
		     a*tex1D<float>(tex[k], ind  ) +
		     b*tex1D<float>(tex[k], ind+1)
		     ) * p0 * plm[Ilm(l, m)] * M[i] * fac0;
	  
	  
	  coef[(2*n+0)*(nend - nbeg) + i] = v * cosp;
	  coef[(2*n+1)*(nend - nbeg) + i] = v * sinp;
	}
      }
    }
  }
}

void SphericalBasis::determine_coefficients_cuda(const Matrix& expcoef)
{
  std::cout << std::scientific;

  // Sort particles and get coefficient size
  std::pair<unsigned, unsigned> lohi = cC->CudaSortByLevel(mlevel, multistep);

  unsigned int N         = lohi.second - lohi.first;
  unsigned int stride     = 1;
  unsigned int gridSize  = N/BLOCK_SIZE/stride;

  if (gridSize>128) {
    stride = N/BLOCK_SIZE/128 + 1;
    gridSize = N/BLOCK_SIZE/stride;
  }

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

  // Make pointers for the device
  //
  cudaParticle        *S = thrust::raw_pointer_cast(&cC->cuda_particles[0]);

  float               *C = thrust::raw_pointer_cast(&dN_coef[0]);
  float               *D = thrust::raw_pointer_cast(&dc_coef[0]);
  float               *M = thrust::raw_pointer_cast(&m_d[0]);
  float               *A = thrust::raw_pointer_cast(&a_d[0]);
  float               *P = thrust::raw_pointer_cast(&p_d[0]);
  float               *L = thrust::raw_pointer_cast(&plm_d[0]);
  int                 *I = thrust::raw_pointer_cast(&i_d[0]);

  cudaTextureObject_t *T = thrust::raw_pointer_cast(&t_d[0]);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(float);


  // For debugging
  //
  if (false) {
    testConstants<<<1, 1>>>();
    
    static bool firstime = true;
    testTexture<<<1, 1>>>(T, nmax);
    firstime == false;
  }

  std::vector<float> coefs((Lmax+1)*(Lmax+1)*nmax);

  // Do the work
  //
  coordKernel<<<gridSize, BLOCK_SIZE>>>(S, M, A, P, L, I, Lmax, stride,
					lohi.first, lohi.second, rmax);

  cudaDeviceSynchronize();


  for (int l=0; l<=Lmax; l++) {
    for (int m=0; m<=l; m++) {
      coefKernel<<<gridSize, BLOCK_SIZE>>>(C, T, M, A, P, L, I, stride,
					   l, m, Lmax, nmax,
					   lohi.first, lohi.second);
      cudaDeviceSynchronize();

      int osize = nmax*2;
      reduceSum<float, BLOCK_SIZE><<<gridSize, BLOCK_SIZE, sMemSize>>>(D, C, osize, N);

      cudaDeviceSynchronize();

      // Finish the reduction
      //
      for (size_t j=0; j<nmax; j++) {
	coefs[Ilmn(l, m, 'c', j, nmax)] =
	  thrust::reduce(dc_coef.begin() + gridSize*(2*j+0),
			 dc_coef.begin() + gridSize*(2*j+1));
	cudaDeviceSynchronize();
	if (m>0)
	  coefs[Ilmn(l, m, 's', j, nmax)] =
	    thrust::reduce(dc_coef.begin() + gridSize*(2*j+1),
			   dc_coef.begin() + gridSize*(2*j+2));
	cudaDeviceSynchronize();
      }
    }
  }

  // DEBUG
  //
  if (true) {
    std::cout << "L=M=0 coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Ilmn(0, 0, 'c', n, nmax)]
		<< std::setw(16) << expcoef[0][n+1]
		<< std::endl;
    }
  }

  //
  // TEST comparison of coefficients
  //
  bool compare = false;

  if (compare) {

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

