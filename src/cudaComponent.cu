// -*- C++ -*-

#include <Component.H>
#include "expand.h"
#include "cudaParticle.cuH"

#include <boost/make_shared.hpp>

unsigned Component::cudaStreamData::totalInstances=0;

Component::cudaStreamData::cudaStreamData()
{
  // Not sure why this breaks thrust, but it does . . .
  /*
  cuda_safe_call(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking),
		 __FILE__, __LINE__,
		 "Component::cudaStreamData: error creating stream");
  */

  // Need blocking until thrust bug in binary search is fixed
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
    
    thrust::device_vector<cudaParticle>::iterator
      pbeg = cr->cuda_particles.begin(),
      pend = cr->cuda_particles.end();
    
    thrust::sort(exec, pbeg, pend, LessCudaLev());
    
    // Convert from cudaParticle to a flat vector to prevent copying
    // the whole particle structure in getBound.  Perhaps this
    // temporary should be part of the data storage structure?
    //
    thrust::device_vector<unsigned> lev(cr->cuda_particles.size());

    thrust::transform(exec, pbeg, pend, lev.begin(), cuPartToLevel());

    // Perform in the sort on the int vector of levels on the GPU
    //
    thrust::device_vector<unsigned>::iterator it;
    if (binary_search_workaround) {
      cudaStreamSynchronize(cr->stream);
      it  = thrust::lower_bound(lev.begin(), lev.end(), minlev);
    } else {
      // This should work but doesn't.  See:
      // https://github.com/NVIDIA/thrust/pull/1104
      //
      it  = thrust::lower_bound(exec, lev.begin(), lev.end(), minlev);
    }
      
    ret.first = thrust::distance(lev.begin(), it);

				// Wait for completion before memcpy
    cudaStreamSynchronize(cr->stream);
				// If maxlev==multistep: upper bound
				// is at end, so skip explicit computation
    if (maxlev < multistep) {

      if (binary_search_workaround) {
	it  = thrust::upper_bound(lev.begin(), lev.end(), maxlev);
      } else {
	it  = thrust::upper_bound(exec, lev.begin(), lev.end(), maxlev);
      }

      ret.second = thrust::distance(lev.begin(), it);

      cudaStreamSynchronize(cr->stream);
    } else {
      ret.second = thrust::distance(pbeg, pend);
    }
    
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

void Component::ParticlesToCuda(PartMap::iterator beg, PartMap::iterator fin)
{
  auto npart = std::distance(beg, fin);
  
  if (host_particles.capacity()<npart) host_particles.reserve(npart);
  host_particles.resize(npart);

  hostPartItr hit = host_particles.begin();
  for (auto pit=beg; pit!=fin; pit++) {
    ParticleHtoD(pit->second, *(hit++));
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


void Component::CudaToParticles(hostPartItr beg, hostPartItr end)
{
  for (hostPartItr v=beg; v!=end; v++) {
    cudaParticle & p = *v;
    ParticleDtoH(p, particles[p.indx]);
  }
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
  if (false) {
				// Don't need this now, since this
				// duplicates zeroing performed on
				// host
    // Loop over bunches
    //
    size_t psize  = particles.size();

    PartMap::iterator begin = particles.begin();
    PartMap::iterator first = begin;
    PartMap::iterator last  = begin;
    PartMap::iterator end   = particles.end();

    if (psize <= bunchSize) last = end;
    else std::advance(last, bunchSize);

    Component::cuRingType cr = *cuRing.get();

    while (thrust::distance(first, last)) {
    
      // Copy particles to host vector
      //
      ParticlesToCuda(first, last);
      
      cr->first = host_particles.begin();
      cr->last  = host_particles.end();
      
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

      // Do copy from host to component
      //
      CudaToParticles(cr->first, cr->last);
      
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
  }

}
