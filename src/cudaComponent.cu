#include <Component.H>
#include <Bounds.cuH>
#include "expand.h"

#include <boost/make_shared.hpp>

unsigned Component::cudaStreamData::totalInstances=0;

Component::cudaStreamData::cudaStreamData()
{
  // cuda_safe_call(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking),
  cuda_safe_call(cudaStreamCreate(&stream),
		 __FILE__, __LINE__,
		 "Component::cudaStreamData: error creating stream");
  instance = totalInstances++;
}

Component::cudaStreamData::~cudaStreamData()
{
  cuda_safe_call(cudaStreamDestroy(stream), __FILE__, __LINE__,
		 "Component::cudaStreamData: error destroying stream");
  totalInstances--;
}

void Component::cuda_initialize()
{
  // Initialize streams
  //
  cuRingData.resize(cuStreams);
  cuRing = boost::make_shared<cuRingType>(cuRingData);
}

struct LevelFunctor
{
  __host__ __device__
  unsigned operator()(const cudaParticle &p) const
  {
    return p.level;
  }
};


std::pair<unsigned int, unsigned int>
Component::CudaSortByLevel(Component::cuRingType cr, int minlev, int maxlev)
{
  std::pair<unsigned, unsigned> ret;

  try {
    auto exec = thrust::cuda::par.on(cr->stream);

    thrust::sort(exec, cr->cuda_particles.begin(), cr->cuda_particles.end(), LessCudaLev());

    cudaParticle temp;

    thrust::device_vector<cudaParticle>::iterator
      pbeg = cr->cuda_particles.begin(),
      pend = cr->cuda_particles.end();
    
    // Get positions of level boundaries
    //
    temp.level = minlev;
    ret.first  = thrust::lower_bound(thrust::cuda::par, pbeg, pend, temp, LessCudaLev()) - pbeg;

    temp.level = maxlev;
    ret.second = thrust::upper_bound(thrust::cuda::par, pbeg, pend, temp, LessCudaLev()) - pbeg;
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

void Component::ParticlesToCuda()
{
  auto npart = Particles().size();
  
  if (host_particles.capacity()<npart) host_particles.reserve(npart);
  host_particles.resize(npart);

  hostPartItr hit = host_particles.begin();
  for (auto pit : Particles()) {
    ParticleHtoD(pit.second, *(hit++));
  }
}

void Component::HostToDev(Component::cuRingType cr)
{
  auto npart = thrust::distance(cr->first, cr->last);
  
  if (npart) {		  // Don't bother trying to copy zero particles

    if (cr->cuda_particles.capacity()<npart) cr->cuda_particles.reserve(npart);
    cr->cuda_particles.resize(npart);
  
    cudaMemcpyAsync(thrust::raw_pointer_cast(&cr->cuda_particles[0]),
		    thrust::raw_pointer_cast(&(*cr->first)),
		    npart*sizeof(cudaParticle),
		    cudaMemcpyHostToDevice, cr->stream);
  }
}

void Component::DevToHost(Component::cuRingType cr)
{
  auto npart = thrust::distance(cr->first, cr->last);
  
  if (npart) {		  // Don't bother trying to copy zero particles

    cudaMemcpyAsync(thrust::raw_pointer_cast(&(*cr->first)),
		    thrust::raw_pointer_cast(&cr->cuda_particles[0]),
		    npart*sizeof(cudaParticle),
		    cudaMemcpyDeviceToHost, cr->stream);
  }
}


void Component::CudaToParticles()
{
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
  return;			// Don't need this now, since this
				// duplicates zeroing performed on
				// host
  // Copy particles to host vector
  //
  ParticlesToCuda();

  // Loop over bunches
  //
  size_t psize  = host_particles.size();

  Component::hostPartItr begin = host_particles.begin();
  Component::hostPartItr first = begin;
  Component::hostPartItr last  = begin;
  Component::hostPartItr end   = host_particles.end();

  if (psize <= bunchSize) last = end;
  else std::advance(last, bunchSize);

  Component::cuRingType cr = *cuRing.get();

  while (thrust::distance(first, last)) {
    
    cr->first = first;
    cr->last  = last;

    // Copy bunch to device
    //
    HostToDev(cr);

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
    DevToHost(cr);

    // Advance iterators
    //
    first = last;
    size_t nadv = thrust::distance(first, end);
    if (nadv <= bunchSize) last = end;
    else thrust::advance(last, bunchSize);

    // Advance stream iterators
    //
    cr++;
  }

  CudaToParticles();
}
