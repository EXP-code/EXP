// -*- C++ -*-

#include <cudaUtil.cuH>
#include <cudaReduce.cuH>

#include <Component.H>
#include "UserSat.H"

// Global symbols
//
__device__ __constant__
cuFP_t userSatMass, userSatCore2, userSatCen[3], userSatPos[3];

__device__ __constant__
bool userSatShadow;

__global__ void
userSatForceKernel(dArray<cudaParticle> P, dArray<int> I,
		   int stride, PII lohi)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i     = tid*stride + n;	// Index in the stride
    int npart = i + lohi.first;	// Particle index

    if (npart < lohi.second) {
      
#ifdef BOUNDS_CHECK
      if (npart>=P._s) printf("out of bounds: %s:%d\n", __FILE__, __LINE__);
#endif
      cudaParticle & p = P._v[I._v[npart]];
      
      cuFP_t rr = userSatCore2;
      for (int k=0; k<3; k++) {
	cuFP_t f = p.pos[k] - userSatCen[k] - userSatPos[k];
	rr += f*f;
      }

      rr = pow(rr, -0.5);
    
      cuFP_t ffac = -userSatMass*rr*rr*rr;

      // Add acceration
      for (int k=0; k<3; k++)
	p.acc[k] += ffac*(p.pos[k] - userSatCen[k] - userSatPos[k]);

      p.potext += -userSatMass*rr;

      // Add the shadow satellite
      if (userSatShadow) {
	rr = userSatCore2;
	for (int k=0; k<3; k++) {
	  cuFP_t f = p.pos[k] - userSatCen[k] + userSatPos[k];
	  rr += f*f;
	}

	rr = pow(rr, -0.5);
	
	ffac = -userSatMass*rr*rr*rr;
	
	// Add acceration
	for (int k=0; k<3; k++)
	  p.acc[k] += ffac*(p.pos[k] - userSatCen[k] + userSatPos[k]);

	p.potext += -userSatMass*rr;
      }

    } // Particle index block

  } // END: stride loop

}


void UserSat::determine_acceration_and_potential_cuda()
{
  // Sanity check
  //
  int nbodies = cC->Number();
  if (nbodies != static_cast<int>(cC->Particles().size())) {
    std::cerr << "UserSat: ooops! number=" << nbodies
	      << " but particle size=" << cC->Particles().size() << endl;
    nbodies = static_cast<int>(cC->Particles().size());
  }
  
  if (nbodies==0) {		// Return if there are no particles
    if (verbose and zbflag) {
      cout << "Process " << myid << ": in UserSat, nbodies=0" 
	   << " for Component <" << cC->name << "> at T=" << tnow
	   << endl;
      zbflag = false;
    }
    return;
  }

  zbflag = true;

  double rs[3];

  if (traj_type==circ) {
    double phi = phase + omega*tnow;
    rs[0] = r0*cos(phi);
    rs[1] = r0*sin(phi);
    rs[2] = 0.0;
  }
  else
    traj->get_satellite_orbit(tnow - toffset, &rs[0]);

  double satmass = mass * 
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    
  if (shadow) satmass *= 0.5;

  if (orbit && myid==0 && mlevel==0 && tnow>tlast) {
    std::ofstream out (orbfile.c_str(), ios::app);
    if (out) {
      out << setw(15) << tnow;
      for (int k=0; k<3; k++) out << setw(15) << rs[k];
      out << endl;
      tlast = tnow;
    } else {
      std::cout << "Error opening trajectory file: " << orbfile << endl;
    }
  }

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cC->cudaDevice);
  cuda_check_last_error_mpi("cudaGetDeviceProperties", __FILE__, __LINE__, myid);

  // Stream structure iterators
  //
  auto cr = cC->cuStream;

  // Assign expansion center
  //
  std::vector<cuFP_t> ctr, sps;

  for (auto v : component->getCenter(Component::Inertial))
    ctr.push_back(v);

  for (int k=0; k<3; k++) {
    sps[k] = rs[k];
    if (pinning) ctr[k] += c0->com[k];
  }

  cuFP_t cuSatCore2 = core * core, cuSatMass = satmass;

  cuda_safe_call(cudaMemcpyToSymbol(userSatCore2, &cuSatCore2, sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userSatCore");
  
  cuda_safe_call(cudaMemcpyToSymbol(userSatMass, &cuSatMass, sizeof(cuFP_t),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userSatMass");
  
  cuda_safe_call(cudaMemcpyToSymbol(userSatCen, &ctr[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userSatCen");

  cuda_safe_call(cudaMemcpyToSymbol(userSatPos, &sps[0], sizeof(cuFP_t)*3,
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userSatPos");

  cuda_safe_call(cudaMemcpyToSymbol(userSatShadow, &shadow,  sizeof(bool),
				    size_t(0), cudaMemcpyHostToDevice),
		 __FILE__, __LINE__, "Error copying userSatShadow");

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

    // Shared memory size for the reduction
    //
    int sMemSize = BLOCK_SIZE * sizeof(cuFP_t);
    
    // Do the work
    //
    userSatForceKernel<<<gridSize, BLOCK_SIZE, sMemSize, cr->stream>>>
      (toKernel(cr->cuda_particles), toKernel(cr->indx1), stride, lohi);
  }
}
