#include <Component.H>
#include <Cylinder.H>

#include <cudaReduce.cuH>

// Global symbols for coordinate transformation
//
__device__ __constant__
float cylRscale, cylHscale, cylXmin, cylXmax, cylYmin, cylYmax, cylDxi, cylDyi;

__device__ __constant__
int   cylNumx, cylNumy, cylCmap;

__host__ __device__
int Imn(int m, char cs, int n, int nmax)
{
  int ret = 0;

  if (m==0) ret = n;
  else ret = (2*m - 1 + (cs=='s' ? 1 : 0))*nmax + n;

  // Verbose sanity check
  if (ret >= (2*m+1)*nmax) {
    printf("Imn oab: %4d %4d %4d [%4d : %4d ]\n", m, n, ret, (2*m+1)*nmax, nmax);
  }

  return ret;
}

__global__
void testConstantsCyl()
{
  printf("** Rscale = %f\n", cylRscale);
  printf("** Hscale = %f\n", cylHscale);
  printf("** Xmin   = %f\n", cylXmin);
  printf("** Xmax   = %f\n", cylXmax);
  printf("** Ymin   = %f\n", cylYmin);
  printf("** Ymax   = %f\n", cylYmax);
  printf("** Dxi    = %f\n", cylDxi);
  printf("** Dyi    = %f\n", cylDyi);
  printf("** Numx   = %d\n", cylNumx);
  printf("** Numy   = %d\n", cylNumy);
  printf("** Cmap   = %d\n", cylCmap);
}

				// R coordinate transformation
__device__
float cu_r_to_xi_cyl(float r)
{
  float ret;

  if (cylCmap==1) {
    ret = (r/cylRscale - 1.0)/(r/cylRscale + 1.0);
  } else {
    ret = r;
  }    

  return ret;
}
    
__device__
float cu_xi_to_r_cyl(float xi)
{
  float ret;

  if (cylCmap==1) {
    ret = (1.0 + xi)/(1.0 - xi) * cylRscale;
  } else {
    ret = xi;
  }

  return ret;
}

__device__
float cu_d_xi_to_r_cyl(float xi)
{
  float ret;

  if (cylCmap==1) {
    ret = 0.5*(1.0 - xi)*(1.0 - xi) / cylRscale;
  } else {
    ret = 1.0;
  }

  return ret;
}

				// Z coordinate transformation
__device__
float cu_z_to_y_cyl(float z)
{ return z/(fabs(z)+FLT_MIN)*asinh(fabs(z/cylHscale)); }

__device__
float cu_y_to_z_cyl(double y)
{ return cylHscale*sinh(y); }

__device__
float cu_d_y_to_z_cyl(float y)
{ return cylHscale*cosh(y); }


void Cylinder::initialize_mapping_constants()
{
  // Copy constants to device
  //
  
  cudaMappingConstants f = getCudaMappingConstants();

  cuda_safe_call(cudaMemcpyToSymbol(cylRscale, &f.rscale, sizeof(float), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuRscale");

  cuda_safe_call(cudaMemcpyToSymbol(cylHscale, &f.hscale, sizeof(float), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuHscale");

  cuda_safe_call(cudaMemcpyToSymbol(cylXmin,   &f.xmin,   sizeof(float), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuXmin");

  cuda_safe_call(cudaMemcpyToSymbol(cylXmax,   &f.xmax,   sizeof(float), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuXmax");

  cuda_safe_call(cudaMemcpyToSymbol(cylDxi,    &f.dxi,    sizeof(float), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuDxi");

  cuda_safe_call(cudaMemcpyToSymbol(cylNumx,   &f.numx,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuNumx");

  cuda_safe_call(cudaMemcpyToSymbol(cylYmin,   &f.ymin,   sizeof(float), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuYmin");

  cuda_safe_call(cudaMemcpyToSymbol(cylYmax,   &f.ymax,   sizeof(float), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuYmax");

  cuda_safe_call(cudaMemcpyToSymbol(cylDyi,    &f.dyi,    sizeof(float), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuDxi");

  cuda_safe_call(cudaMemcpyToSymbol(cylNumy,   &f.numy,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuNumy");

  cuda_safe_call(cudaMemcpyToSymbol(cylCmap,   &f.cmap,   sizeof(int),   size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cuCmap");
}


__global__
void testCoordCyl(dArray<float> mass, dArray<float> phi,
		  dArray<float> Xfac, dArray<float> Yfac,
		  dArray<int> IndX, dArray<int> IndY,
		  PII lohi)
{
  printf("**\n** Coordinate test\n");
  
  int N = lohi.second - lohi.first;

  printf("%8s %13s %13s %13s %13s %5s %5s\n", "#", "mass", "phi", "a", "c", "ix", "iy");

  for (int i=0; i<3; i++) {
    printf("%8d %13.7e %13.7e %13.7e %13.7e %5d %5d\n",
	   i, mass._v[i], phi._v[i], Xfac._v[i], Yfac._v[i],
	   IndX._v[i], IndY._v[i]);
  }

  for (int i=N-4; i<N; i++) {
    printf("%8d %13.7e %13.7e %13.7e %13.7e %5d %5d\n",
	   i, mass._v[i], phi._v[i], Xfac._v[i], Yfac._v[i],
	   IndX._v[i], IndY._v[i]);
  }

  printf("**\n");
}

__global__
void testTextureCyl(dArray<cudaTextureObject_t> tex, int nmax)
{
  printf("** DEVICE 2d texture compare\n");
  for (int k=0; k<10; k++) {
    for (int i : {0, 1, 126, 127}) 
      for (int j : {0, 1, 126, 127}) 
	printf("%5d %5d %5d %13.7e\n", k, i, j, tex2D<float>(tex._v[j], i, j));
  }
}

__global__ void testParticles(dArray<cudaParticle> in)
{
  for (int k=0; k<4; k++)
    printf("%5d %13.5e %13.5e %13.5e %13.5e\n",
	   k, in._v[k].mass, in._v[k].pos[0], in._v[k].pos[1], in._v[k].pos[2]);

  for (int k=in._s-4; k<in._s; k++)
    printf("%5d %13.5e %13.5e %13.5e %13.5e\n",
	   k, in._v[k].mass, in._v[k].pos[0], in._v[k].pos[1], in._v[k].pos[2]);
}

__global__ void coordKernelCyl
(dArray<cudaParticle> in, dArray<float> mass, dArray<float> phi,
 dArray<float> Xfac, dArray<float> Yfac,
 dArray<int> IndX, dArray<int> IndY,
 unsigned int stride, PII lohi, float rmax)
{
  // Thread ID
  //
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

  // Centering logic, TBD
  /*
    vector<double> ctr;
    if (mix) mix->getCenter(ctr);
  */
  float ctr[3] {0.0f, 0.0f, 0.0f};

  for (int n=0; n<stride; n++) {
    int i = tid*stride + n;	// Particle counter
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {

      if (npart>=in._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

      cudaParticle p = in._v[npart];
    
      float xx = p.pos[0] - ctr[0];
      float yy = p.pos[1] - ctr[1];
      float zz = p.pos[2] - ctr[2];
      
      float r2 = (xx*xx + yy*yy + zz*zz);
      float r  = sqrt(r2) + FSMALL;
      
      if (i>=mass._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

      mass._v[i] = -1.0;
      
      if (r<rmax) {
	
	mass._v[i] = p.mass;
	
	phi._v[i] = atan2(yy, xx);

	if (i>=phi._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	// Interpolation indices
	//
	float X  = (cu_r_to_xi_cyl(r) - cylXmin)/cylDxi;
	float Y  = (cu_z_to_y_cyl(zz) - cylYmin)/cylDyi;

	int indX = floor(X);
	int indY = floor(Y);
	
	if (indX<0) indX = 0;
	if (indX>cylNumx-2) indX = cylNumx - 2;
	
	if (indY<0) indY = 0;
	if (indY>cylNumy-2) indY = cylNumy - 2;
	
	Xfac._v[i] = float(indX+1) - X;
	IndX._v[i] = indX;

	Yfac._v[i] = float(indY+1) - Y;
	IndY._v[i] = indY;

	if (Xfac._v[i]<-0.5 or Xfac._v[i]>1.5) printf("X off grid: x=%f\n", X);
	if (Yfac._v[i]<-0.5 or Yfac._v[i]>1.5) printf("Y off grid: y=%f\n", Y);

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
 dArray<float> Mass, dArray<float> Phi,
 dArray<float> Xfac, dArray<float> Yfac,
 dArray<int> indX, dArray<int> indY,
 int stride, int m, unsigned int nmax, PII lohi)
{
  // Thread ID
  //
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  // Total number of particles to be evaluated
  //
  const unsigned int N = lohi.second - lohi.first;

  const float norm = -4.0*M_PI;	// Biorthogonality factor

  for (int istr=0; istr<stride; istr++) {

    int i = tid*stride + istr;	// Particle counter

    if (i<N) {			// Allow for grid padding

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
	if (i>=Yfac._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	float c00 = delx0*dely0;
	float c10 = delx1*dely0;
	float c01 = delx0*dely1;
	float c11 = delx1*dely1;

	int indx = indX._v[i];
	int indy = indY._v[i];

	if (i>=indX._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
	if (i>=indY._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	for (int n=0; n<nmax; n++) {

	  // Texture maps are packed sequentially as
	  // ---------------------------------------
	  //     potC, rforceC, zforceC
	  //
	  // for m=0, and
	  //
	  //     potC, rforceC, zforceC, potS, rforceS, zforceS
	  //
	  // for m>0
	  //
	  // The inner loop increments radial order and the outer loop
	  // increments harmonic order.
	  //
	  if (m==0) {

	    int k = 3*n;	// Indices: 0, 3, 6, ... (3*nmax-3) which pulls out potC

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

	    // Pulls out potC tables at indices
	    // --------------------------------
	    // For m=1, indicies:  3*nmax,  3*nmax+6,  3+nmax+12, ... ,  9*nmax - 6
	    // For m=2, indicies:  9*nmax,  9*nmax+6,  9*nmax+12, ... , 15*nmax - 6
	    // For m=3, indicies: 15*nmax, 15*nmax+6, 15*nmax+12, ... , 21*nmax - 6
	    // ...
	    int k = 3*(2*m - 1)*nmax + 6*n; 

	    const float d00  = tex2D<float>(tex._v[k  ], indx,   indy  );
	    const float d10  = tex2D<float>(tex._v[k  ], indx+1, indy  );
	    const float d01  = tex2D<float>(tex._v[k  ], indx,   indy+1);
	    const float d11  = tex2D<float>(tex._v[k  ], indx+1, indy+1);

	    if (k>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	    // potS tables are offset from potC tables by +3
	    //
	    const float e00  = tex2D<float>(tex._v[k+3], indx,   indy  );
	    const float e10  = tex2D<float>(tex._v[k+3], indx+1, indy  );
	    const float e01  = tex2D<float>(tex._v[k+3], indx,   indy+1);
	    const float e11  = tex2D<float>(tex._v[k+3], indx+1, indy+1);

	    if (k+1>=tex._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);

	    coef._v[(2*n+0)*N + i] = (c00*d00 + c10*d10 + c01*d01 + c11*d11) * cosp * norm * mass;
	    coef._v[(2*n+1)*N + i] = (c00*e00 + c10*e10 + c01*e01 + c11*e11) * sinp * norm * mass;

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
	       int stride, unsigned int mmax, unsigned int nmax, PII lohi,
	       float rmax, float cylmass)
{
  // Thread ID
  //
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

  // Maximum radius squared
  //
  const float rmax2 = rmax*rmax;

  // Centering logic, TBD
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
      
      float xx  = p.pos[0] - ctr[0];
      float yy  = p.pos[1] - ctr[1];
      float zz  = p.pos[2] - ctr[2];
      
      float phi = atan2(yy, xx);
      float R2  = xx*xx + yy*yy;
      float  R  = sqrt(R2) + FSMALL;
      
      const float ratmin = 0.75;
      const float maxerf = 3.0;
      const float midpt  = ratmin + 0.5*(1.0 - ratmin);
      const float rsmth  = 0.5*(1.0 - ratmin)/maxerf;

      float ratio = sqrt( (R2 + zz*zz)/rmax2 );
      float mfactor = 1.0, frac = 1.0, cfrac = 0.0;

      if (ratio >= 1.0) {
	cfrac      = 1.0 - mfactor;
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

	float X  = (cu_r_to_xi_cyl(R) - cylXmin)/cylDxi;
	float Y  = (cu_z_to_y_cyl(zz) - cylYmin)/cylDyi;

	int indX = floor(X);
	int indY = floor(Y);
	
	float delx0 = float(indX+1) - X;
	float dely0 = float(indY+1) - Y;

	if (delx0<0.0 or delx0>1.0) printf("X off grid: x=%f\n", delx0);
	if (dely0<0.0 or dely0>1.0) printf("Y off grid: y=%f\n", dely0);

	float delx1 = 1.0 - delx0;
	float dely1 = 1.0 - dely0;
      
	float c00 = delx0*dely0;
	float c10 = delx1*dely0;
	float c01 = delx0*dely1;
	float c11 = delx1*dely1;

	float cos1 = cos(phi);
	float sin1 = sin(phi);

	float ccos = 1.0;
	float ssin = 0.0;

	for (int mm=0; mm<=mmax; mm++) {

	  for (int n=0; n<nmax; n++) {
      
	    float fac0 = coef._v[Imn(mm, 'c', n, nmax)];
	    float fac1 = fac0 * ccos;
	    float fac2 = fac0 * ssin;
      
	    // Texture table index
	    //
	    int k = 3*n;
	    if (mm) k = 3*(2*mm - 1)*nmax + 6*n;

	    pp += fac1 *
	      (
	       tex2D<float>(tex._v[k  ], indX,   indY  ) * c00 +
	       tex2D<float>(tex._v[k  ], indX+1, indY  ) * c10 +
	       tex2D<float>(tex._v[k  ], indX,   indY+1) * c01 +
	       tex2D<float>(tex._v[k  ], indX+1, indY+1) * c11 
	       );
	    
	    fr += fac1 *
	      (
	       tex2D<float>(tex._v[k+1], indX,   indY  ) * c00 +
	       tex2D<float>(tex._v[k+1], indX+1, indY  ) * c10 +
	       tex2D<float>(tex._v[k+1], indX,   indY+1) * c01 +
	       tex2D<float>(tex._v[k+1], indX+1, indY+1) * c11 
	       );
      
	    fz += fac1 *
	      (
	       tex2D<float>(tex._v[k+2], indX,   indY  ) * c00 +
	       tex2D<float>(tex._v[k+2], indX+1, indY  ) * c10 +
	       tex2D<float>(tex._v[k+2], indX,   indY+1) * c01 +
	       tex2D<float>(tex._v[k+2], indX+1, indY+1) * c11 
	       );
	    
	    fp += fac2 * mm *
	      (
	       tex2D<float>(tex._v[k  ], indX,   indY  ) * c00 +
	       tex2D<float>(tex._v[k  ], indX+1, indY  ) * c10 +
	       tex2D<float>(tex._v[k  ], indX,   indY+1) * c01 +
	       tex2D<float>(tex._v[k  ], indX+1, indY+1) * c11 
	       );
      
      
	    if (mm) {
	
	      float fac0 =  coef._v[Imn(mm, 's', n, nmax)];
	      float fac1 =  fac0 * ssin;
	      float fac2 = -fac0 * ccos;

	      pp += fac1 *
		(
		 tex2D<float>(tex._v[k+3], indX,   indY  ) * c00 +
		 tex2D<float>(tex._v[k+3], indX+1, indY  ) * c10 +
		 tex2D<float>(tex._v[k+3], indX,   indY+1) * c01 +
		 tex2D<float>(tex._v[k+3], indX+1, indY+1) * c11 
		 );
	      
	      fr += fac1 *
		(
		 tex2D<float>(tex._v[k+4], indX,   indY  ) * c00 +
		 tex2D<float>(tex._v[k+4], indX+1, indY  ) * c10 +
		 tex2D<float>(tex._v[k+4], indX,   indY+1) * c01 +
		 tex2D<float>(tex._v[k+4], indX+1, indY+1) * c11 
		 );
	      
	      fz += fac1 *
		(
		 tex2D<float>(tex._v[k+5], indX,   indY  ) * c00 +
		 tex2D<float>(tex._v[k+5], indX+1, indY  ) * c10 +
		 tex2D<float>(tex._v[k+5], indX,   indY+1) * c01 +
		 tex2D<float>(tex._v[k+5], indX+1, indY+1) * c11 
		 );
	      
	      fp += fac2 * mm *
		(
		 tex2D<float>(tex._v[k+3], indX,   indY  ) * c00 +
		 tex2D<float>(tex._v[k+3], indX+1, indY  ) * c10 +
		 tex2D<float>(tex._v[k+3], indX,   indY+1) * c01 +
		 tex2D<float>(tex._v[k+3], indX+1, indY+1) * c11 
		 );
	      
	    }
	  }
	  
	  // Trig recursion to squeeze avoid internal FP fct call
	  //
	  float cosM = ccos;
	  float sinM = ssin;

	  ccos = cosM * cos1 - sinM * sin1;
	  ssin = sinM * cos1 + cosM * sin1;
	}

	in._v[npart].acc[0] += ( fr*xx/R - fp*yy/R2 ) * frac;
	in._v[npart].acc[1] += ( fr*yy/R + fp*xx/R2 ) * frac;
	in._v[npart].acc[2] += fz * frac;
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


void Cylinder::determine_coefficients_cuda()
{
  /*
  std::cout << " ** BEFORE initialize" << std::endl;
  testParticles<<<1, 1>>>(toKernel(cC->cuda_particles));
  std::cout << " **" << std::endl;
  */

  static bool cuda_firstime = true;
  if (cuda_firstime and myid==0) {
    initialize_cuda();
    initialize_mapping_constants();
    cuda_firstime = false;
  }

  /*
  std::cout << " ** AFTER initialize" << std::endl;
  testParticles<<<1, 1>>>(toKernel(cC->cuda_particles));
  std::cout << " **" << std::endl;
  */

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

  // unsigned int Nthread = gridSize*BLOCK_SIZE;

  std::cout << "**" << std::endl
	    << "** N      = " << N           << std::endl
	    << "** I low  = " << lohi.first  << std::endl
	    << "** I high = " << lohi.second << std::endl
	    << "** Stride = " << stride      << std::endl
	    << "** Block  = " << BLOCK_SIZE  << std::endl
	    << "** Grid   = " << gridSize    << std::endl
	    << "**" << std::endl;

  // Create space for coefficient reduction
  //
  thrust::device_vector<float> dN_coef(2*ncylorder*N);
  thrust::device_vector<float> dc_coef(2*ncylorder*gridSize);
  thrust::device_vector<float> df_coef(2*ncylorder);

  // Texture objects
  //
  thrust::device_vector<cudaTextureObject_t> t_d = tex;

  // Space for coordinate arrays
  //
  thrust::device_vector<float> m_d(N), X_d(N), Y_d(N), p_d(N);
  thrust::device_vector<int>   iX_d(N), iY_d(N);

  // Shared memory size for the reduction
  //
  int sMemSize = BLOCK_SIZE * sizeof(float);

  // For debugging (set to false to disable)
  //
  static bool firstime = false;

  if (firstime) {
    testConstantsCyl<<<1, 1>>>();
    cudaDeviceSynchronize();

    testTextureCyl<<<1, 1>>>(toKernel(t_d), ncylorder);
    cudaDeviceSynchronize();

    firstime = false;
  }

  std::vector<float> coefs((2*mmax+1)*ncylorder);

  thrust::counting_iterator<int> index_begin(0);
  thrust::counting_iterator<int> index_end(gridSize*2*ncylorder);

  // Maximum radius on grid
  //
  float rmax = rcylmax * acyl;

  // Do the work
  //
  // testParticles<<<1, 1>>>(toKernel(cC->cuda_particles));

				// Compute the coordinate
				// transformation
				// 
  coordKernelCyl<<<gridSize, BLOCK_SIZE>>>
    (toKernel(cC->cuda_particles), toKernel(m_d), toKernel(p_d),
     toKernel(X_d), toKernel(Y_d), toKernel(iX_d), toKernel(iY_d),
     stride, lohi, rmax);

  /*
  testCoordCyl<<<1, 1>>>(toKernel(m_d),  toKernel(p_d),
			 toKernel(X_d),  toKernel(Y_d),
			 toKernel(iX_d), toKernel(iY_d),
			 lohi);
  */
				// Compute the coefficient
				// contribution for each order
  int osize = ncylorder*2;	// 
  for (int m=0; m<=mmax; m++) {
    coefKernelCyl<<<gridSize, BLOCK_SIZE>>>
      (toKernel(dN_coef), toKernel(t_d), toKernel(m_d), toKernel(p_d),
       toKernel(X_d), toKernel(Y_d), toKernel(iX_d), toKernel(iY_d),
       stride, m, ncylorder, lohi);
    
				// Begin the reduction per grid block
				//
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
    for (size_t j=0; j<ncylorder; j++) {
      coefs[Imn(m, 'c', j, ncylorder)] = ret[2*j];
      if (m>0) coefs[Imn(m, 's', j, ncylorder)] = ret[2*j+1];
    }
  }

  // DEBUG
  //
  if (false) {
    std::cout << "M=0 coefficients" << std::endl;
    for (size_t n=0; n<ncylorder; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << Imn(0, 'c', n, ncylorder)
		<< std::setw(16) << coefs[Imn(0, 'c', n, ncylorder)]
		<< std::setw(16) << ortho->get_coef(0, n, 'c')
		<< std::endl;
    }

    std::cout << "M=1c coefficients" << std::endl;
    for (size_t n=0; n<ncylorder; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << Imn(1, 'c', n, ncylorder)
		<< std::setw(16) << coefs[Imn(1, 'c', n, ncylorder)]
		<< std::setw(16) << ortho->get_coef(1, n, 'c')
		<< std::endl;
    }

    std::cout << "M=1s coefficients" << std::endl;
    for (size_t n=0; n<ncylorder; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << Imn(1, 's', n, ncylorder)
		<< std::setw(16) << coefs[Imn(1, 's', n, ncylorder)]
		<< std::setw(16) << ortho->get_coef(1, n, 's')
		<< std::endl;
    }

    std::cout << "M=2c coefficients" << std::endl;
    for (size_t n=0; n<ncylorder; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << Imn(2, 'c', n, ncylorder)
		<< std::setw(16) << coefs[Imn(2, 'c', n, ncylorder)]
		<< std::setw(16) << ortho->get_coef(2, n, 'c')
		<< std::endl;
    }
    
    std::cout << "M=2s coefficients" << std::endl;
    for (size_t n=0; n<ncylorder; n++) {
      std::cout << std::setw(4)  << n
		<< std::setw(4)  << Imn(2, 's', n, ncylorder)
		<< std::setw(16) << coefs[Imn(2, 's', n, ncylorder)]
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
	for (int n=0; n<ncylorder; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 'c';
	  elem.d = ortho->get_coef(m, n, 'c');
	  elem.f = coefs[Imn(m, 'c', n, ncylorder)];
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-4) test /= fabs(elem.d);
	  
	  compare[test] = elem;
	    
	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 'c'
	      << std::setw( 5) << Imn(m, 'c', n, ncylorder)
	      << std::setw(14) << elem.d
	      << std::setw(14) << elem.f
	      << std::endl;
	}

      } else {
	for (int n=0; n<ncylorder; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 'c';
	  elem.d = ortho->get_coef(m, n, 'c');
	  elem.f = coefs[Imn(m, 'c', n, ncylorder)];

	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 'c'
	      << std::setw( 5) << Imn(m, 'c', n, ncylorder)
	      << std::setw(14) << elem.d
	      << std::setw(14) << elem.f
	      << std::endl;
	  
	  double test = fabs(elem.d - elem.f);
	  if (fabs(elem.d)>1.0e-4) test /= fabs(elem.d);

	  compare[test] = elem;
	}

	for (int n=0; n<ncylorder; n++) {
	  elem.m = m;
	  elem.n = n;
	  elem.cs = 's';
	  elem.d = ortho->get_coef(m, n, 'c');
	  elem.f = coefs[Imn(m, 's', n, ncylorder)];

	  out << std::setw( 5) << m
	      << std::setw( 5) << n
	      << std::setw( 5) << 's'
	      << std::setw( 5) << Imn(m, 's', n-1, ncylorder)
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

  // Maximum radius on grid
  //
  float rmax = rcylmax * acyl;

  // Do the work
  //
  forceKernelCyl<<<gridSize, BLOCK_SIZE, sMemSize>>>
    (toKernel(cC->cuda_particles), toKernel(dev_coefs), toKernel(t_d),
     stride, mmax, ncylorder, lohi, rmax, cylmass);
}

void Cylinder::HtoD_coefs()
{
  host_coefs.resize((2*mmax+1)*ncylorder);

  // m loop
  //
  for (int m=0; m<=mmax; m++) {
    
    // n loop
    //
    for (int n=0; n<ncylorder; n++) {
      host_coefs[Imn(m, 'c', n, ncylorder)] = ortho->get_coef(m, n, 'c');
      if (m>0) host_coefs[Imn(m, 's', n, ncylorder)] = ortho->get_coef(m, n, 's');
    }
  }

  dev_coefs = host_coefs;
}


void Cylinder::DtoH_coefs()
{
  host_coefs = dev_coefs;

  // m loop
  //
  for (int m=0; m<=mmax; m++) {
    
    // n loop
    //
    for (int n=0; n<ncylorder; n++) {
      ortho->get_coef(m, n, 'c') = host_coefs[Imn(m, 'c', n, ncylorder)];
      if (m>0)
	ortho->get_coef(m, n, 's') = host_coefs[Imn(m, 's', n, ncylorder)];
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

    
