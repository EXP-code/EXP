#include <Component.H>

// Define for checking cuda_particles: values and counts
//
// #define BIG_DEBUG

std::pair<unsigned int, unsigned int>
Component::CudaSortByLevel(int minlev, int maxlev)
{
  std::pair<unsigned, unsigned> ret;

 try
   {
     thrust::sort(cuda_particles.begin(), cuda_particles.end(), LessCudaLev());

     cudaParticle temp;

     thrust::device_vector<cudaParticle>::iterator
       pbeg = cuda_particles.begin(),
       pend = cuda_particles.end();

     // Get positions of level boundaries
     //
     temp.level = minlev;
     ret.first  = thrust::lower_bound(pbeg, pend, temp, LessCudaLev()) - pbeg;
     
     temp.level = maxlev;
     ret.second = thrust::upper_bound(pbeg, pend, temp, LessCudaLev()) - pbeg;
   }
 catch(std::bad_alloc &e)
  {
    std::cerr << "Ran out of memory while sorting" << std::endl;
    exit(-1);
  }
 catch(thrust::system_error &e)
   {
     std::cerr << "Some other error happened during sort, lower_bound, or upper_bound:" << e.what() << std::endl;
     exit(-1);
   }
 
  return ret;
}

void Component::CudaSortBySequence()
{
  thrust::device_vector<cudaParticle>::iterator
    pbeg = cuda_particles.begin(),
    pend = cuda_particles.end();

  thrust::sort(pbeg, pend, LessCudaSeq());
}

#ifdef BIG_DEBUG
__global__ void testParticles(cudaParticle* in, int N)
{
  for (int k=0; k<4; k++)
    printf("%5d %13.5e %13.5e %13.5e %13.5e\n",
	   k, in[k].mass, in[k].pos[0], in[k].pos[1], in[k].pos[2]);

  for (int k=N-4; k<N; k++)
    printf("%5d %13.5e %13.5e %13.5e %13.5e\n",
	   k, in[k].mass, in[k].pos[0], in[k].pos[1], in[k].pos[2]);
}
#endif

void Component::ParticlesToCuda()
{
#ifdef BIG_DEBUG
  static unsigned count = 0;
  std::cout << std::string(72, '-') << std::endl
	    << "---- Component name: " << name
	    << " [#" << count++ << "]" << std::endl
	    << "---- Comp particle size: " << particles.size() << std::endl
	    << "---- Host particle size: " << host_particles.size()
	    << " [" << std::hex << thrust::raw_pointer_cast(host_particles.data()) << "] before" << std::endl << std::dec
	    << "---- Cuda particle size: " << cuda_particles.size()
	    << " [" << hex << thrust::raw_pointer_cast(cuda_particles.data()) << "] before" << std::endl << std::dec;
#endif
  
  host_particles.resize(particles.size());
  thrust::host_vector<cudaParticle>::iterator it = host_particles.begin();
  for (auto v : particles) {
    ParticleHtoD(v.second, *(it++));
  }
  
  cuda_particles = host_particles;

#ifdef BIG_DEBUG
  /*
  testParticles<<<1, 1>>>(thrust::raw_pointer_cast(cuda_particles.data()),
			  cuda_particles.size());
  */

  std::cout << std::string(72, '-') << std::endl
	    << "---- Host particle size: " << host_particles.size()
	    << " [" << std::hex << thrust::raw_pointer_cast(host_particles.data()) << "] after" << std::endl << std::dec
	    << "---- First mass: " << host_particles[0].mass << std::endl
	    << "---- Last mass: " << host_particles[host_particles.size()-1].mass << std::endl
	    << "---- Cuda particle size: " << cuda_particles.size()
	    << " [" << hex << thrust::raw_pointer_cast(cuda_particles.data()) << "] after" << std::endl << std::dec
	    << std::string(72, '-') << std::endl;
#endif
}


void Component::CudaToParticles()
{
  host_particles = cuda_particles;
  for (auto v : host_particles) ParticleDtoH(v, particles[v.indx]);
}

struct cudaZeroAcc : public thrust::unary_function<cudaParticle, cudaParticle>
{
  __host__ __device__
  cudaParticle operator()(cudaParticle& p)
  {
    for (size_t k=0; k<3; k++) p.acc[k] = 0.0;
    p.pot = p.potext = 0.0;
    return p;
  }
};

void Component::ZeroPotAccel(int minlev)
{
  if (multistep) {
    std::pair<unsigned int, unsigned int> ret = CudaSortByLevel(minlev, multistep);

    thrust::transform(cuda_particles.begin()+ret.first, cuda_particles.end(),
		      cuda_particles.begin()+ret.first, cudaZeroAcc());
  } else {
    thrust::transform(cuda_particles.begin(), cuda_particles.end(),
		      cuda_particles.begin(), cudaZeroAcc());
  }
}
