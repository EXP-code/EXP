#include <boost/make_shared.hpp>
#include <Component.H>

// Define >0 for checking cuda_particles: values and counts
// Set to zero for production
//
#define BIG_DEBUG 0

unsigned Component::cudaStreamData::totalInstances=0;


void Component::cuda_initialize()
{
  // Initialize 3 streams
  //
  cuRingData.resize(3);
  cuRing = boost::make_shared<cuRingType>(cuRingData);
}

std::pair<unsigned int, unsigned int>
Component::CudaSortByLevel(Component::cuRingType cr, int minlev, int maxlev)
{
  std::pair<unsigned, unsigned> ret;

  try {
    thrust::sort(thrust::cuda::par.on(cr->stream),
		 cr->cuda_particles.begin(), cr->cuda_particles.end(),
		 LessCudaLev());

    cudaParticle temp;

    thrust::device_vector<cudaParticle>::iterator
      pbeg = cr->cuda_particles.begin(),
      pend = cr->cuda_particles.end();
    
    // Get positions of level boundaries
    //
    temp.level = minlev;
    ret.first  = thrust::lower_bound(thrust::cuda::par.on(cr->stream),
				     pbeg, pend, temp, LessCudaLev()) - pbeg;
    
    temp.level = maxlev;
    ret.second = thrust::upper_bound(thrust::cuda::par.on(cr->stream),
				     pbeg, pend, temp, LessCudaLev()) - pbeg;
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

void Component::CudaSortBySequence(Component::cuRingType cr)
{
  thrust::device_vector<cudaParticle>::iterator
    pbeg = cr->cuda_particles.begin(),
    pend = cr->cuda_particles.end();

  thrust::sort(thrust::cuda::par.on(cr->stream), pbeg, pend, LessCudaSeq());
}

void Component::ParticlesToCuda(Component::cuRingType cr,
				PartMap::iterator first, PartMap::iterator last)
{
  std::cout << std::scientific;

  auto npart = std::distance(first, last);
  
#if BIG_DEBUG > 0
  static unsigned count = 0;
  std::cout << std::string(72, '-')        << std::endl
	    << "---- Component name: "     << name
	    << " [#" << count++ << "]"     << std::endl
	    << "---- Particle count: "     << npart  << std::endl
	    << std::string(72, '-')        << std::endl
	    << "---- Host particle size: " << cr->host_particles.size()
	    << " [" << std::hex << thrust::raw_pointer_cast(cr->host_particles.data()) << "] before" << std::endl << std::dec
	    << "---- Cuda particle size: " << cuda_particles.size()
	    << " [" << hex << thrust::raw_pointer_cast(cr->cuda_particles.data()) << "] before" << std::endl << std::dec
	    << "---- Size of real: " << sizeof(cuFP_t) << std::endl
	    << "---- Size of cudaParticle: " << sizeof(cudaParticle)
	    << ", native=" << 12*sizeof(cuFP_t) + 2*sizeof(unsigned) << std::endl;
#endif

  if (cr->host_particles.capacity()<npart) cr->host_particles.reserve(npart);
  cr->host_particles.resize(npart);

  if (cr->cuda_particles.capacity()<npart) cr->cuda_particles.reserve(npart);
  cr->cuda_particles.resize(npart);
  
  thrust::host_vector<cudaParticle>::iterator hit = cr->host_particles.begin();
  for (auto pit=first; pit!=last; pit++) {
    ParticleHtoD(pit->second, *(hit++));
  }
  
#if BIG_DEBUG > 1
  static unsigned cnt = 0;
  std::ostringstream sout;
  sout << "test." << cnt++;

  std::ofstream tmp(sout.str());
  std::copy(cr->host_particles.begin(), cr->host_particles.end(),
	    std::ostream_iterator<cudaParticle>(tmp, "\n") );

  std::cout << "[host] BEFORE copy" << std::endl;
  std::copy(cr->host_particles.begin(), cr->host_particles.begin()+5,
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );

  std::copy(cr->host_particles.end()-5, cr->host_particles.end(),
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );
#endif

  cudaMemcpyAsync(thrust::raw_pointer_cast(cr->cuda_particles.data()),
		  thrust::raw_pointer_cast(cr->host_particles.data()),
		  cr->cuda_particles.size()*sizeof(cudaParticle),
		  cudaMemcpyHostToDevice, cr->stream);

#if BIG_DEBUG > 1
  std::cout << "[cuda] AFTER copy" << std::endl;
  std::copy(cr->cuda_particles.begin(), cr->cuda_particles.begin()+5,
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );

  std::copy(cr->cuda_particles.end()-5, cr->cuda_particles.end(),
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );

  std::cout << "SORT particles by seq" << std::endl;
  CudaSortBySequence(cp);

  std::cout << "[cuda] AFTER sort" << std::endl;

  std::copy(cr->cuda_particles.begin(), cr->cuda_particles.begin()+5,
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );

  std::copy(cr->cuda_particles.end()-5, cr->cuda_particles.end(),
	    std::ostream_iterator<cudaParticle>(std::cout, "\n") );
#endif

#if BIG_DEBUG > 0
  std::cout << std::string(72, '-') << std::endl
	    << "---- Host particle size: " << cr->host_particles.size()
	    << " [" << std::hex << thrust::raw_pointer_cast(cr->host_particles.data()) << "] after" << std::endl << std::dec
#if BIG_DEBUG > 1
	    << "---- First mass: " << cr->host_particles[0].mass << std::endl
	    << "---- Last mass:  " << cr->host_particles[cr->host_particles.size()-1].mass << std::endl
#endif
	    << "---- Cuda particle size: " << cr->cuda_particles.size()
	    << " [" << hex << thrust::raw_pointer_cast(cr->cuda_particles.data()) << "] after" << std::endl << std::dec
	    << std::string(72, '-') << std::endl;
#endif
}


void Component::CudaToParticles(Component::cuRingType cr)
{
  auto csize = cr->cuda_particles.size();
  if (csize > cr->host_particles.size()) cr->host_particles.reserve(csize);
  cr->host_particles.resize(csize);  

  cudaMemcpyAsync(thrust::raw_pointer_cast(cr->host_particles.data()),
		  thrust::raw_pointer_cast(cr->cuda_particles.data()),
		  cr->host_particles.size()*sizeof(cudaParticle),
		  cudaMemcpyDeviceToHost, cr->stream);


  for (auto v : cr->host_particles) ParticleDtoH(v, particles[v.indx]);
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
  // Loop over bunches
  //
  size_t psize  = Particles().size();
  
  PartMap::iterator begin = Particles().begin();
  PartMap::iterator first = begin;
  PartMap::iterator last  = begin;
  PartMap::iterator end   = Particles().end();

  if (psize <= bunchSize) last = end;
  else std::advance(last, bunchSize);

  Component::cuRingType cr = *cuRing.get();

  while (std::distance(first, last)) {
    
    // Copy bunch to device
    //
    ParticlesToCuda(cr, first, last);

    if (multistep) {
      std::pair<unsigned int, unsigned int>
	ret = CudaSortByLevel(cr, minlev, multistep);
      
      thrust::transform(thrust::cuda::par.on(cr->stream),
			cr->cuda_particles.begin()+ret.first, cr->cuda_particles.end(),
			cr->cuda_particles.begin()+ret.first, cudaZeroAcc());
    } else {
      thrust::transform(thrust::cuda::par.on(cr->stream),
			cr->cuda_particles.begin(), cr->cuda_particles.end(),
			cr->cuda_particles.begin(), cudaZeroAcc());
    }

    // Copy device to host
    //
    CudaToParticles(cr);

    // Advance iterators
    //
    first = last;
    size_t nadv = std::distance(first, end);
    if (nadv <= bunchSize) last = end;
    else std::advance(last, bunchSize);

    // Advance stream iterators
    //
    cr++;
  }


}
