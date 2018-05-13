#include <Component.H>
#include <Cylinder.H>

#include <cudaReduce.cuH>

__host__ __device__
int Imn(int m, char cs, int n, int nmax)
{
  int ret = 0;

  if (m==0) ret = n;
  else ret = (2*m - 1 + (cs=='s' ? 1 : 0))*nmax + n;

  if (ret >= (2*m+1)*nmax) {
    printf("Imn oab: %4d %4d %4d [%4d : %4d ]\n", m, n, ret, (2*m+1)*nmax, nmax);
  }

  return ret;
}

__global__
void testConstantsCyl()
{
  printf("** Rscale = %f\n", cuRscale);
  printf("** Hscale = %f\n", cuHscale);
  printf("** Xmin   = %f\n", cuXmin);
  printf("** Xmax   = %f\n", cuXmax);
  printf("** Ymin   = %f\n", cuYmin);
  printf("** Ymax   = %f\n", cuYmax);
  printf("** Dxi    = %f\n", cuDxi);
  printf("** Dyi    = %f\n", cuDxi);
  printf("** Numx   = %d\n", cuNumx);
  printf("** Numy   = %d\n", cuNumy);
  printf("** Cmap   = %d\n", cuCmap);
}

__device__
float cu_r_to_xi_cyl(float r)
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
float cu_xi_to_r_cyl(float xi)
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
float cu_d_xi_to_r_cyl(float xi)
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

				// Z coordinate transformation

__device__
float cu_z_to_y_cyl(float z)
{ return z/(fabs(z)+FLT_MIN)*asinh(fabs(z/cuHscale)); }

__device__
float cu_y_to_z_cyl(double y)
{ return cuHscale*sinh(y); }

__device__
float cu_d_y_to_z_cyl(float y)
{ return cuHscale*cosh(y); }


void Cylinder::initialize_mapping_constants()
{
  // Copy constants to device
  //
  
  cudaMappingConstants f = getCudaMappingConstants();

  cuda_safe_call(cudaMemcpyToSymbol(cuRscale, &f.rscale, sizeof(float), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuRscale");

  cuda_safe_call(cudaMemcpyToSymbol(cuHscale, &f.hscale, sizeof(float), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuHscale");

  cuda_safe_call(cudaMemcpyToSymbol(cuXmin, &f.xmin, sizeof(float), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuXmin");

  cuda_safe_call(cudaMemcpyToSymbol(cuXmax, &f.xmax, sizeof(float), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuXmax");

  cuda_safe_call(cudaMemcpyToSymbol(cuDxi, &f.dxi, sizeof(float), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuDxi");

  cuda_safe_call(cudaMemcpyToSymbol(cuNumr, &f.numr, sizeof(int), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuNumr");

  cuda_safe_call(cudaMemcpyToSymbol(cuCmap, &f.cmap, sizeof(int), size_t(0), cudaMemcpyHostToDevice), __FILE__, __LINE__, "Error copying cuCmap");
}


__global__
void testTextureCyl(dArray<cudaTextureObject_t> tex, int nmax)
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

__global__ void coordKernelCyl
(dArray<cudaParticle> in, dArray<float> mass, dArray<float> phi,
 dArray<float> Xfac, dArray<float> Yfac,
 dArray<int> IndX, dArray<int> IndY,
 unsigned int stride, PII lohi, float rmax)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

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
	
	phi._v[i] = atan2(yy, xx);

	if (i>=phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	float X  = (cu_r_to_xi_cyl(r) - cuXmin)/cuDxi;
	float Y  = (cu_z_to_y_cyl(zz) - cuYmin)/cuDyi;

	int indX = floor(X);
	int indY = floor(Y);
	
	if (indX<0) indX = 0;
	if (indX>cuNumx-2) indX = cuNumr - 2;
	
	if (indY<0) indY = 0;
	if (indY>cuNumy-2) indY = cuNumr - 2;
	
	Xfac._v[i] = float(indX+1) - X;
	if (Xfac._v[i]<0.0 or Xfac._v[i]>1.0)
	  printf("X off grid: x=%f\n", X);
	IndX._v[i] = indX;

	Yfac._v[i] = float(indY+1) - Y;
	if (Yfac._v[i]<0.0 or Yfac._v[i]>1.0)
	  printf("Y off grid: y=%f\n", Y);
	IndY._v[i] = indY;

	if (i>=Xfac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=IndX._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=Yfac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=IndY._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
      }
    }
  }
}


__global__ void coefKernelCyl
(dArray<float> coef, dArray<cudaTextureObject_t> tex,
 dArray<float> Mass, dArray<float> Phi, dArray<float> Xfac, dArray<float> Yfac,
 dArray<int> indX, dArray<int> indY, int stride, int m, unsigned int nmax, PII lohi)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;
  const unsigned int N = lohi.second - lohi.first;
  const float norm = -4.0*M_PI;

  for (int istr=0; istr<stride; istr++) {

    int i = tid*stride + istr;

    if (i<N) {

      float mass = Mass._v[i];
      
      if (i>=Mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
      
      if (mass>0.0) {

	if (i>=Phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	float phi  = Phi._v[i];
	float cosp = cos(phi*m);
	float sinp = sin(phi*m);
	
	// Do the interpolation
	//
	float delx0 = Xfac._v[i];
	float dely0 = Yfac._v[i];
	float delx1 = 1.0 - delx0;
	float dely1 = 1.0 - dely0;

	if (i>=Xfac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=indX._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=Yfac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=indY._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	float c00 = delx0*dely0;
	float c10 = delx1*dely0;
	float c01 = delx0*dely1;
	float c11 = delx1*dely1;

	int indx = indX._v[i];
	int indy = indY._v[i];

	for (int n=0; n<nmax; n++) {

	  if (m==0) {

	    int k = 3*n;

	    if (n>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	    // Fetch the values from the texture

	    const float d00  = tex2D<float>(tex._v[k], indx,   indy  );
	    const float d10  = tex2D<float>(tex._v[k], indx+1, indy  );
	    const float d01  = tex2D<float>(tex._v[k], indx,   indy+1);
	    const float d11  = tex2D<float>(tex._v[k], indx+1, indy+1);

	    if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	    coef._v[(2*n+0)*N + i] = (c00*d00 + c10*d10 + c01*d01 + c11*d11) * norm * mass;
	    
	    if ((2*n+0)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	  } else {

	    int k = 3*(2*m - 1)*nmax + 6*n;

	    const float d00  = tex2D<float>(tex._v[k  ], indx,   indy  );
	    const float d10  = tex2D<float>(tex._v[k  ], indx+1, indy  );
	    const float d01  = tex2D<float>(tex._v[k  ], indx,   indy+1);
	    const float d11  = tex2D<float>(tex._v[k  ], indx+1, indy+1);

	    if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	    const float e00  = tex2D<float>(tex._v[k+1], indx,   indy  );
	    const float e10  = tex2D<float>(tex._v[k+1], indx+1, indy  );
	    const float e01  = tex2D<float>(tex._v[k+1], indx,   indy+1);
	    const float e11  = tex2D<float>(tex._v[k+1], indx+1, indy+1);

	    if (k+1>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	    coef._v[(2*n+0)*N + i] = (c00*d00 + c10*d10 + c01*d01 + c11*d11) * cosp * norm * mass;
	    coef._v[(2*n+1)*N + i] = (c00*d00 + c10*d10 + c01*d01 + c11*d11) * sinp * norm * mass;

	    if ((2*n+0)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	    if ((2*n+1)*N+i>=coef._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	  }
	} // norder loop
      } // mass value check
    } // particle index check
  } // stride loop
}

__global__ void
forceKernelCyl(dArray<cudaParticle> in, dArray<float> coef,
	       dArray<cudaTextureObject_t> tex,
	       int stride, unsigned int mmax, unsigned int nmax, PII lohi, float cylmass)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

  /*
    vector<double> ctr;
    if (mix) mix->getCenter(ctr);
  */
  float ctr[3] {0.0f, 0.0f, 0.0f};

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {
      
      if (npart>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

      cudaParticle p = in._v[npart];
      
      float xx = p.pos[0] - ctr[0];
      float yy = p.pos[1] - ctr[1];
      float zz = p.pos[2] - ctr[2];
      
      float phi   = atan2(yy, xx);
      float R2    = xx*xx + yy*yy + FSMALL;
      float  R    = sqrt(R2);
      
      const double ratmin = 0.75;
      const double maxerf = 3.0;
      const double midpt = ratmin + 0.5*(1.0 - ratmin);
      const double rsmth = 0.5*(1.0 - ratmin)/maxerf;

      float ratio = sqrt( (R2 + zz*zz)/R2 );
      float mfactor = 1.0, frac = 1.0, cfrac = 0.0;

      if (ratio >= 1.0) {
	cfrac      = 1.0 - mfactor;
	in._v[npart].acc[0] = 0.0;
	in._v[npart].acc[1] = 0.0;
	in._v[npart].acc[2] = 0.0;

      } else if (ratio > ratmin) {
	frac  = 0.5*(1.0 - erf( (ratio - midpt)/rsmth )) * mfactor;
	cfrac = 1.0 - frac;
      } else {
	frac  = mfactor;
      }

      float fr = 0.0;
      float fz = 0.0;
      float fp = 0.0;
      float pp = 0.0;
      
      if (ratio < 1.0) {

	float X  = (cu_r_to_xi_cyl(R) - cuXmin)/cuDxi;
	float Y  = (cu_z_to_y_cyl(zz) - cuYmin)/cuDyi;

	int indX = floor(X);
	int indY = floor(Y);
	
	float delx0 = float(indX+1) - X;
	if (delx0<0.0 or delx0>1.0)
	  printf("X off grid: x=%f\n", delx0);

	float dely0 = float(indY+1) - Y;
	if (dely0<0.0 or dely0>1.0)
	  printf("Y off grid: y=%f\n", dely0);

	float delx1 = 1.0 - delx0;
	float dely1 = 1.0 - delx1;
      
	float c00 = delx0*dely0;
	float c10 = delx1*dely0;
	float c01 = delx0*dely1;
	float c11 = delx1*dely1;

	for (int mm=0; mm<=mmax; mm++) {

	  float ccos = cos(phi*mm);
	  float ssin = sin(phi*mm);

	  for (int n=0; n<nmax; n++) {
      
	    int ic = n;
	    if (mm) ic +=  (2*mm - 1)*nmax;

	    float fac = coef._v[ic] * ccos;
      
	    int k = 3*n;
	    if (mm) k = 3*(2*mm - 1)*nmax + 6*n;

	    pp += fac *
	      (
	       tex2D<float>(tex._v[k  ], indX,   indY  ) * c00 +
	       tex2D<float>(tex._v[k  ], indX+1, indY  ) * c10 +
	       tex2D<float>(tex._v[k  ], indX,   indY+1) * c01 +
	       tex2D<float>(tex._v[k  ], indX+1, indY+1) * c11 
	       );
	    
	    fr += fac *
	      (
	       tex2D<float>(tex._v[k+1], indX,   indY  ) * c00 +
	       tex2D<float>(tex._v[k+1], indX+1, indY  ) * c10 +
	       tex2D<float>(tex._v[k+1], indX,   indY+1) * c01 +
	       tex2D<float>(tex._v[k+1], indX+1, indY+1) * c11 
	       );
      
	    fz += fac *
	      (
	       tex2D<float>(tex._v[k+2], indX,   indY  ) * c00 +
	       tex2D<float>(tex._v[k+2], indX+1, indY  ) * c10 +
	       tex2D<float>(tex._v[k+2], indX,   indY+1) * c01 +
	       tex2D<float>(tex._v[k+2], indX+1, indY+1) * c11 
	       );
	    
	    fac = coef._v[ic] * ssin;
	    
	    fp += fac * mm *
	      (
	       tex2D<float>(tex._v[k  ], indX,   indY  ) * c00 +
	       tex2D<float>(tex._v[k  ], indX+1, indY  ) * c10 +
	       tex2D<float>(tex._v[k  ], indX,   indY+1) * c01 +
	       tex2D<float>(tex._v[k  ], indX+1, indY+1) * c11 
	       );
      
      
	    if (mm) {
	
	      ic +=  nmax;

	      fac = coef._v[ic] * ssin;
	      
	      pp += fac *
		(
		 tex2D<float>(tex._v[k+3], indX,   indY  ) * c00 +
		 tex2D<float>(tex._v[k+3], indX+1, indY  ) * c10 +
		 tex2D<float>(tex._v[k+3], indX,   indY+1) * c01 +
		 tex2D<float>(tex._v[k+3], indX+1, indY+1) * c11 
		 );
	      
	      fr += fac *
		(
		 tex2D<float>(tex._v[k+4], indX,   indY  ) * c00 +
		 tex2D<float>(tex._v[k+4], indX+1, indY  ) * c10 +
		 tex2D<float>(tex._v[k+4], indX,   indY+1) * c01 +
		 tex2D<float>(tex._v[k+4], indX+1, indY+1) * c11 
		 );
	      
	      fz += fac *
		(
		 tex2D<float>(tex._v[k+5], indX,   indY  ) * c00 +
		 tex2D<float>(tex._v[k+5], indX+1, indY  ) * c10 +
		 tex2D<float>(tex._v[k+5], indX,   indY+1) * c01 +
		 tex2D<float>(tex._v[k+5], indX+1, indY+1) * c11 
		 );
	      
	      fac = -coef._v[ic] * ccos;
	
	      fp += fac * mm *
		(
		 tex2D<float>(tex._v[k+3], indX,   indY  ) * c00 +
		 tex2D<float>(tex._v[k+3], indX+1, indY  ) * c10 +
		 tex2D<float>(tex._v[k+3], indX,   indY+1) * c01 +
		 tex2D<float>(tex._v[k+3], indX+1, indY+1) * c11 
		 );
	      
	    }
	  }
	  
	}

	in._v[npart].acc[0] = ( fr*xx/R - fp*yy/R2 ) * frac;
	in._v[npart].acc[1] = ( fr*yy/R + fp*xx/R2 ) * frac;
	in._v[npart].acc[2] = fz * frac;
      }

      if (ratio > ratmin) {

	float r3 = R2 + zz*zz;
	pp = -cylmass/sqrt(r3);	// -M/r
	fr = pp/r3;		// -M/r^3

	in._v[npart].acc[0] += xx*fr * cfrac;
	in._v[npart].acc[1] += yy*fr * cfrac;
	in._v[npart].acc[2] += zz*fr * cfrac;
      }

      in._v[npart].pot += pp;

    } // Particle index block

  } // END: stride loop

}


void Cylinder::determine_coefficients_cuda(const Matrix& expcoef)
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

  // unsigned int Nthread = gridSize*BLOCK_SIZE;

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
  thrust::device_vector<float> df_coef(2*nmax);

  // Texture objects
  //
  thrust::device_vector<cudaTextureObject_t> t_d = tex;

  // Space for Legendre coefficients 
  //
  thrust::device_vector<float> m_d(N), X_d(N), Y_d(N), p_d(N);
  thrust::device_vector<int>   iX_d(N), iY_d(N);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(float);

  // For debugging
  //
  if (false) {
    testConstantsCyl<<<1, 1>>>();
    
    static bool firstime = true;
    testTextureCyl<<<1, 1>>>(toKernel(t_d), nmax);
    firstime == false;
  }

  std::vector<float> coefs((2*mmax+1)*nmax);

  thrust::counting_iterator<int> index_begin(0);
  thrust::counting_iterator<int> index_end(gridSize*2*nmax);

  // Do the work
  //
				// Compute the coordinate
				// transformation
				// 
  coordKernelCyl<<<gridSize, BLOCK_SIZE>>>
    (toKernel(cC->cuda_particles), toKernel(m_d), toKernel(p_d),
     toKernel(X_d), toKernel(Y_d), toKernel(iX_d), toKernel(iY_d),
     stride, lohi, rcylmax);

				// Compute the coefficient
				// contribution for each order
  for (int m=0; m<=mmax; m++) {
    coefKernelCyl<<<gridSize, BLOCK_SIZE>>>
      (toKernel(dN_coef), toKernel(t_d), toKernel(m_d), toKernel(p_d),
       toKernel(X_d), toKernel(Y_d), toKernel(iX_d), toKernel(iY_d),
       stride, m, nmax, lohi);
    
				// Begin the reduction per grid block
    int osize = nmax*2;	// 
    reduceSum<float, BLOCK_SIZE><<<gridSize, BLOCK_SIZE, sMemSize>>>
      (toKernel(dc_coef), toKernel(dN_coef), osize, N);
      
				// Finish the reduction for this order
				// in parallel
    thrust::reduce_by_key
      (
       thrust::make_transform_iterator(index_begin, key_functor(gridSize)),
       thrust::make_transform_iterator(index_end,   key_functor(gridSize)),
       dc_coef.begin(), thrust::make_discard_iterator(), df_coef.begin()
       );
    
    thrust::host_vector<float> ret = df_coef;
    for (size_t j=0; j<nmax; j++) {
      coefs[Imn(m, 'c', j, nmax)] = ret[2*j];
      if (m>0) coefs[Imn(m, 's', j, nmax)] = ret[2*j+1];
    }
  }

  // DEBUG
  //
  if (false) {
    std::cout << "M=0 coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Imn(0, 'c', n, nmax)]
		<< std::setw(16) << ortho->get_coef(0, n, 'c')
		<< std::endl;
    }

    std::cout << "M=1c coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Imn(1, 'c', n, nmax)]
		<< std::setw(16) << ortho->get_coef(1, n, 'c')
		<< std::endl;
    }

    std::cout << "M=1s coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Imn(1, 's', n, nmax)]
		<< std::setw(16) << ortho->get_coef(1, n, 's')
		<< std::endl;
    }

    std::cout << "M=2c coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Imn(2, 'c', n, nmax)]
		<< std::setw(16) << ortho->get_coef(2, n, 'c')
		<< std::endl;
    }
    
    std::cout << "M=2s coefficients" << std::endl;
    for (size_t n=0; n<nmax; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(16) << coefs[Imn(2, 'c', n, nmax)]
		<< std::setw(16) << ortho->get_coef(2, n, 's')
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
      
      int  m;
      int  n;
      
      char cs;
    }
    elem;

    std::map<double, Element> compare;

    std::ofstream out("test.dat");

    // m loop
    for (int m=0; m<=mmax; m++) {
	
      if (m==0) {
	for (int n=0; n<nmax; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 'c';
	  elem.d = ortho->get_coef(m, n, 'c');
	  elem.f = coefs[Imn(m, 'c', n, nmax)];
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-4) test /= fabs(elem.d);
	  
	  compare[test] = elem;
	    
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
	  elem.f = coefs[Imn(m, 'c', n, nmax)];

	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 'c'
	      << std::setw( 5) << Imn(m, 'c', n, nmax)
	      << std::setw(14) << elem.d
	      << std::setw(14) << elem.f
	      << std::endl;
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-4) test /= fabs(elem.d);

	  compare[test] = elem;
	}

	for (int n=0; n<nmax; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 's';
	  elem.d = ortho->get_coef(m, n, 'c');
	  elem.f = coefs[Imn(m, 's', n, nmax)];

	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 's'
	      << std::setw( 5) << Imn(m, 's', n-1, nmax)
	      << std::setw(14) << elem.d
	      << std::setw(14) << elem.f
	      << std::endl;
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-4) test /= fabs(elem.d);
	  
	  compare[test] = elem;
	}
      }
    }
    
    std::map<double, Element>::iterator best = compare.begin();
    std::map<double, Element>::iterator midl = best;
    std::advance(midl, compare.size()/2);
    std::map<double, Element>::reverse_iterator last = compare.rbegin();
    
    std::cout << "Best case: ["
	      << std::setw( 2) << best->second.m << ", "
	      << std::setw( 2) << best->second.n << ", "
	      << std::setw( 2) << best->second.cs << "] = "
	      << std::setw(15) << best->second.d
	      << std::setw(15) << best->second.f
	      << std::setw(15) << fabs(best->second.d - best->second.f)
	      << std::endl;
  
    std::cout << "Mid case:  ["
	      << std::setw( 2) << midl->second.m << ", "
	      << std::setw( 2) << midl->second.n << ", "
	      << std::setw( 2) << midl->second.cs << "] = "
	      << std::setw(15) << midl->second.d
	      << std::setw(15) << midl->second.f
	      << std::setw(15) << fabs(midl->second.d - midl->second.f)
	      << std::endl;
    
    std::cout << "Last case: ["
	      << std::setw( 2) << last->second.m << ", "
	      << std::setw( 2) << last->second.n << ", "
	      << std::setw( 2) << last->second.cs << "] = "
	      << std::setw(15) << last->second.d
	      << std::setw(15) << last->second.f
	      << std::setw(15) << fabs(last->second.d - last->second.f)
	      << std::endl;
  }
}


void Cylinder::determine_acceleration_cuda()
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

  // unsigned int Nthread = gridSize*BLOCK_SIZE;

  std::cout << "**" << std::endl
	    << "** N      = " << N          << std::endl
	    << "** Stride = " << stride     << std::endl
	    << "** Block  = " << BLOCK_SIZE << std::endl
	    << "** Grid   = " << gridSize   << std::endl
	    << "**" << std::endl;

  // Texture objects
  //
  thrust::device_vector<cudaTextureObject_t> t_d = tex;

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(float);

  // Do the work
  //
  forceKernelCyl<<<gridSize, BLOCK_SIZE, sMemSize>>>
    (toKernel(cC->cuda_particles), toKernel(dev_coefs), toKernel(t_d),
     stride, mmax, nmax, lohi, cylmass);
}

void Cylinder::HtoD_coefs(const Matrix& expcoef)
{
  host_coefs.resize((2*mmax+1)*nmax);

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

  dev_coefs = host_coefs;
}


void Cylinder::DtoH_coefs(Matrix& expcoef)
{
  host_coefs = dev_coefs;

  // m loop
  //
  for (int m=0; m<=mmax; m++) {
    
    // n loop
    //
    for (int n=0; n<nmax; n++) {
      ortho->get_coef(m, n, 'c') = host_coefs[Imn(m, 'c', n, nmax)];
      if (m>0)
	ortho->get_coef(m, n, 's') = host_coefs[Imn(m, 's', n, nmax)];
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
    
  // std::cout << "cuda memory freed" << std::endl;
}

void Cylinder::host_dev_force_compare()
{
  // Copy from device
  cC->host_particles = cC->cuda_particles;
  
  std::streamsize ss = std::cout.precision();
  std::cout.precision(4);

  std::cout << std::string(16+14*7, '-') << std::endl
	    << std::setw(8)  << "Index"  << std::setw(8)  << "Level"
	    << std::setw(14) << "ax [d]" << std::setw(14) << "ay [d]"
	    << std::setw(14) << "az [d]" << std::setw(14) << "ax [h]"
	    << std::setw(14) << "ay [h]" << std::setw(14) << "az [h]"
	    << std::setw(14) << "|Del a|/|a|"  << std::endl;

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
      std::cout << std::setw(14) << sqrt(diff/norm) << std::endl;
    }
  
  for (size_t j=0; j<5; j++) 
    {
      size_t i = cC->host_particles.size() - 6 + j;

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
      std::cout << std::setw(14) << sqrt(diff/norm) << std::endl;
    }

  std::cout << std::string(16+14*7, '-') << std::endl;
  std::cout.precision(ss);
}

    
