#include <Component.H>

std::pair<unsigned, unsigned> Component::CudaSortByLevel(int minlev, int maxlev)
{
  std::pair<unsigned, unsigned> ret;

  std::cout << "Size of cuda particle vector: " << cuda_particles.size()
	    << std::endl
	    << "Begin sort by level . . . " << std::flush;

  thrust::sort(cuda_particles.begin(), cuda_particles.end(), LessCudaLev());
  std::cout << "done" << std::endl;

  cudaParticle temp;

  thrust::device_vector<cudaParticle>::iterator
    pbeg = cuda_particles.begin(),
    pend = cuda_particles.end();

  // Get positions of level boundaries
  //
  temp.level = minlev;
  ret.first  =
    thrust::lower_bound(pbeg, pend, temp, LessCudaLev()) - pbeg;

  temp.level = maxlev;
  ret.second = thrust::upper_bound(pbeg, pend, temp, LessCudaLev()) - pbeg;

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
  host_particles.resize(particles.size());
  unsigned cnt = 0;
  for (auto v : particles) {
    ParticleHtoD(v.second, host_particles[cnt++]);
  }
  cuda_particles = host_particles;
}


void Component::CudaToParticles()
{
  host_particles = cuda_particles;
  for (auto v : host_particles) {
    ParticleDtoH(v, particles[v.indx]);
  }
}
