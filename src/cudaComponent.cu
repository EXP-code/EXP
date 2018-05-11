#include <Component.H>

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


void Component::ParticlesToCuda()
{
  host_particles.resize(particles.size());
  thrust::host_vector<cudaParticle>::iterator it = host_particles.begin();
  for (auto v : particles) ParticleHtoD(v.second, *(it++));
  cuda_particles = host_particles;
}


void Component::CudaToParticles()
{
  host_particles = cuda_particles;
  for (auto v : host_particles) ParticleDtoH(v, particles[v.indx]);
}
