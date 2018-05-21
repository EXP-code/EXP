#include <Component.H>

// Define >0 for checking cuda_particles: values and counts
// Set to zero for production
//
#define BIG_DEBUG 1

std::pair<unsigned int, unsigned int>
Component::CudaSortByLevel(int minlev, int maxlev)
{
  std::pair<unsigned, unsigned> ret;

  try {
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
  catch(std::bad_alloc &e) {
    std::cerr << "Ran out of memory while sorting" << std::endl;
    exit(-1);
  }
  catch(thrust::system_error &e) {
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

void Component::ParticlesToCuda()
{
  std::cout << std::scientific;

#if BIG_DEBUG > 0
  static unsigned count = 0;
  std::cout << std::string(72, '-')        << std::endl
	    << "---- Component name: "     << name
	    << " [#" << count++ << "]"     << std::endl
	    << "---- Particle count: "     << particles.size() << std::endl
	    << std::string(72, '-')        << std::endl
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
  
#if BIG_DEBUG > 1
  static unsigned cnt = 0;
  std::ostringstream sout;
  sout << "test." << cnt++;

  std::ofstream tmp(sout.str());
  std::copy(host_particles.begin(), host_particles.end(),
	    std::ostream_iterator<cudaParticle>(tmp, "\n") );

  std::cout << "[host] BEFORE copy" << std::endl;
  std::copy(host_particles.begin(), host_particles.begin()+5,
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );

  std::copy(host_particles.end()-5, host_particles.end(),
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );
#endif

  cuda_particles = host_particles;


#if BIG_DEBUG > 1
  std::cout << "[cuda] AFTER copy" << std::endl;
  std::copy(cuda_particles.begin(), cuda_particles.begin()+5,
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );

  std::copy(cuda_particles.end()-5, cuda_particles.end(),
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );

  std::cout << "SORT particles by seq" << std::endl;
  CudaSortBySequence();

  std::cout << "[cuda] AFTER sort" << std::endl;

  std::copy(cuda_particles.begin(), cuda_particles.begin()+5,
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );

  std::copy(cuda_particles.end()-5, cuda_particles.end(),
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );
#endif

#if BIG_DEBUG > 0
  std::cout << std::string(72, '-') << std::endl
	    << "---- Host particle size: " << host_particles.size()
	    << " [" << std::hex << thrust::raw_pointer_cast(host_particles.data()) << "] after" << std::endl << std::dec
#if BIG_DEBUG > 1
	    << "---- First mass: " << host_particles[0].mass << std::endl
	    << "---- Last mass:  " << host_particles[host_particles.size()-1].mass << std::endl
#endif
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
