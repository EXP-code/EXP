#include <Orient.H>
#include <cudaUtil.cuH>
#include <cudaParticle.cuH>

#include "expand.h"
#include "global.H"

__device__ __constant__
cuFP_t cuCtr[3];

__device__ __constant__
unsigned int cuFlags;

//! Simplified EL3 structure (no time) for use in CUDA kernel code
struct cudaEL3
{
  //! Mass
  cuFP_t M;
  //! Binding energy
  cuFP_t E;			
  //! Angular momentum
  cuFP_t L[3];
  //! Position
  cuFP_t R[3];
} __attribute__((aligned));


//! Copy device structure to host structure with time variable
EL3 cudaToEL3(const cudaEL3& p, double time)
{
  EL3 ret;

  ret.T = time;
  ret.M = p.M;
  ret.E = p.E;
  for (int k=0; k<3; k++) {
    ret.L[k+1] = p.L[k];
    ret.R[k+1] = p.R[k];
  }

  return ret;
}

//! For thrust sorts on sequence
struct LessCudaEL3
{
  __host__ __device__
  bool operator()(const cudaEL3& p1, const cudaEL3& p2) const
  {
    return (p1.E < p2.E);
  }
};

using PII = std::pair<unsigned int, unsigned int>;

__global__ void EL3Kernel
(dArray<cudaParticle> P, dArray<cudaEL3> el3, unsigned int stride, PII lohi)
{
  // The enum from the cpu code
  //
  enum ControlFlags {DIAG=1, KE=2, EXTERNAL=4};

  // Work array
  //
  cuFP_t psa[3];

  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i = tid*stride + n;
    int n = i + lohi.first;

    if (n < lohi.second) {

      cudaParticle & p = P._v[n];
      cudaEL3 & t = el3._v[i];

      cuFP_t v2 = 0.0;
      for (int k=0; k<3; k++) {
	psa[k] = p.pos[k] - cuCtr[k];
	v2 += p.vel[k]*p.vel[k];
      }

      cuFP_t energy = p.pot;
    
      if (cuFlags & KE) energy += 0.5*v2;

      if (cuFlags & EXTERNAL) energy += p.potext;

      cuFP_t mass = p.mass;

      t.E = energy;
      t.M = mass;

      t.L[0] = mass*(psa[1]*p.vel[2] - psa[2]*p.vel[1]);
      t.L[1] = mass*(psa[2]*p.vel[0] - psa[0]*p.vel[2]);
      t.L[2] = mass*(psa[0]*p.vel[1] - psa[1]*p.vel[0]);
      
      t.R[0] = mass*p.pos[0];
      t.R[1] = mass*p.pos[1];
      t.R[2] = mass*p.pos[2];
    }
  }
}

void Orient::accumulate_gpu(double time, Component *c)
{
  // Copy constants to constant memory
  //
  cuFP_t cen[3];
  for (int k=0; k<3; k++) cen[k] = center[k+1];
  
  cuda_safe_call(cudaMemcpyToSymbol(cuCtr, cen, 3*sizeof(cuFP_t), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying center");

  cuda_safe_call(cudaMemcpyToSymbol(cuFlags, &cflags, sizeof(unsigned int), size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying cflags");

  // Prepare bunch loop
  //
  unsigned nbodies = c->Number();
  unsigned tkeep = many/numprocs;

  const unsigned oBunchSize = 200000;
  unsigned int Npacks = nbodies/oBunchSize + 1;

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, c->cudaDevice);

  thrust::device_vector<cudaEL3> devEL3;
  thrust::host_vector<cudaEL3> hostEL3(tkeep);

  // Loop over bunches
  //
  for (int n=0; n<Npacks; n++) {

    PII cur;

    // Current bunch
    //
    cur. first = oBunchSize*n;
    cur.second = oBunchSize*(n+1);
    cur.second = std::min<unsigned int>(cur.second, nbodies);
    
    if (cur.second <= cur.first) break;
    
    // Compute grid
    //
    unsigned int N         = cur.second - cur.first;
    unsigned int stride    = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
    unsigned int gridSize  = N/BLOCK_SIZE/stride;
      
    if (N > gridSize*BLOCK_SIZE*stride) gridSize++;
      
    // Resize storage as needed
    //
    devEL3.resize(N);
      
    // Compute the coordinate transformation
    // 
    EL3Kernel<<<gridSize, BLOCK_SIZE>>>
      (toKernel(c->cuStream->cuda_particles), toKernel(devEL3),
       stride, cur);
    
    // Sort on the device
    //
    thrust::sort(devEL3.begin(), devEL3.end(), LessCudaEL3());

    // Get the first tkeep values
    //
    if (N<tkeep) {		// Bunch smaller than tkeep?
      hostEL3.resize(N);
      thrust::copy(devEL3.begin(), devEL3.end(), hostEL3.begin());
    } else			// First tkeep values
      thrust::copy(devEL3.begin(), devEL3.begin() + tkeep, hostEL3.begin());

    // Copy from cuda to host structure and load the std::set
    //
    for (auto v : hostEL3) angm.insert(cudaToEL3(v, time));
    if (angm.size() > tkeep) {	// Keep only tkeep smallest
      auto it = angm.begin();
      std::advance(it, tkeep);
      angm.erase(it, angm.end());
    }
  }
}
