#include <Component.H>
#include "expand.h"
#include "cudaParticle.cuH"

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>

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
  cuStream = boost::make_shared<cudaStreamData>();
}

struct LevelFunctor
{
  int _t;

  LevelFunctor(int t=0) : _t(t) {}

  __host__ __device__
  int operator()(const cudaParticle &p) const
  {
    return p.lev[_t];
  }
};


Component::I2vec Component::CudaSortLevelChanges(Component::cuSharedStream cr)
{
  // The plan: for the current active level search above and below for
  // particles for correction to coefficient matrix
  //
  // 1. Sort all particles by current level
  // 2. Get indices to range for each level L
  // 3. Within each level L, compute the ranges for changes,
  //    delta L = [-L, multistep-L]
  // 4. For each (L, delta L), compute the coefficient changes and
  //    apply to the appropriate coefficient matrices

  I2vec ret(multistep+1);
  for (auto & v : ret) v.resize(multistep+1);

  try {
    auto exec = thrust::cuda::par.on(cuStream->stream);
    
    thrust::device_vector<cudaParticle>::iterator
      pbeg = cuStream->cuda_particles.begin(),
      pend = cuStream->cuda_particles.end();
    
    if (thrust_binary_search_workaround) {
      cudaStreamSynchronize(cuStream->stream);
      thrust::sort(pbeg, pend, LessCudaLev2());
    } else {
      thrust::sort(exec, pbeg, pend, LessCudaLev2());
    }
    
    pbeg = cuStream->cuda_particles.begin();
    pend = cuStream->cuda_particles.end();

    cudaParticle trg;

    for (int target=0; target<=multistep; target++) {

      trg.lev[0] = target;

      for (int del=0; del<=multistep; del++) {

	if (del==target) {
	  ret[target][del] = {0, 0};
	  continue;
	}
	
	trg.lev[1] = del;

	thrust::device_vector<cudaParticle>::iterator lo, hi;

	if (thrust_binary_search_workaround) {
	  cudaStreamSynchronize(cuStream->stream);
	  lo  = thrust::lower_bound(pbeg, pend, trg, LessCudaLev2());
	} else {
	  lo = thrust::lower_bound(exec, pbeg, pend, trg, LessCudaLev2());
	}
	
	cudaStreamSynchronize(cuStream->stream);

	if (thrust_binary_search_workaround) {
	  hi = thrust::upper_bound(pbeg, pend, trg, LessCudaLev2());
	} else {
	  hi = thrust::upper_bound(exec, pbeg, pend, trg, LessCudaLev2());
	}

	cudaStreamSynchronize(cuStream->stream);

	ret[target][del] = {thrust::distance(pbeg, lo), 
			    thrust::distance(pbeg, hi)};
      }
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
 
  if (false) {
    std::cout << std::string(15*(multistep+1), '-') << std::endl;
    std::cout << "--- " << name << " [" << myid << "]" << std::endl;
    std::cout << std::string(15*(multistep+1), '-') << std::endl;
    for (int m1=0; m1<=multistep; m1++) {
      for (int m2=0; m2<=multistep; m2++) {
	std::ostringstream sout;
	sout << "(" << ret[m1][m2].first << ", " << ret[m1][m2].second << ")";
	std::cout << std::setw(15) << sout.str();
      }
      std::cout << std::endl;
    }
    std::cout << std::string(15*(multistep+1), '-') << std::endl;
  }

  return ret;
}


void Component::CudaSortByLevel()
{
  try {
    auto exec = thrust::cuda::par.on(cuStream->stream);
    
    thrust::device_vector<cudaParticle>::iterator
      pbeg = cuStream->cuda_particles.begin(),
      pend = cuStream->cuda_particles.end();
    
    if (thrust_binary_search_workaround) {
      cudaStreamSynchronize(cuStream->stream);
      thrust::sort(pbeg, pend, LessCudaLev());
    } else {
      thrust::sort(exec, pbeg, pend, LessCudaLev());
    }
    
    // Convert from cudaParticle to a flat vector to prevent copying
    // the whole particle structure in getBound.  Perhaps this
    // temporary should be part of the data storage structure?
    //
    cuStream->levList.resize(cuStream->cuda_particles.size());

    if (thrust_binary_search_workaround) {
      cudaStreamSynchronize(cuStream->stream);
      thrust::transform(pbeg, pend, cuStream->levList.begin(), cuPartToLevel());
    } else {
      thrust::transform(exec, pbeg, pend, cuStream->levList.begin(), cuPartToLevel());
    }
  }
  catch(thrust::system_error &e) {
    std::cerr << "Some other error happened during sort, lower_bound, or upper_bound:" << e.what() << std::endl;
    exit(-1);
  }

}

std::pair<unsigned int, unsigned int>
Component::CudaGetLevelRange(int minlev, int maxlev)
{
  std::pair<unsigned, unsigned> ret;

  try {
    auto exec = thrust::cuda::par.on(cuStream->stream);

    // Get unsigned from input
    //
    unsigned int minl = static_cast<unsigned>(minlev);
    unsigned int maxl = static_cast<unsigned>(maxlev);

    thrust::device_vector<int>::iterator lbeg = cuStream->levList.begin();
    thrust::device_vector<int>::iterator lend = cuStream->levList.end();
    thrust::device_vector<int>::iterator lo, hi;

    if (thrust_binary_search_workaround) {
      cudaStreamSynchronize(cuStream->stream);
      lo  = thrust::lower_bound(lbeg, lend, minl);
    } else {
      lo = thrust::lower_bound(exec, lbeg, lend, minl);
    }
	
    if (thrust_binary_search_workaround) {
      cudaStreamSynchronize(cuStream->stream);
      hi = thrust::upper_bound(lbeg, lend, maxl);
    } else {
      hi = thrust::upper_bound(exec, lbeg, lend, maxl);
    }

    ret.first  = thrust::distance(lbeg, lo);
    ret.second = thrust::distance(lbeg, hi);
  }
  catch(thrust::system_error &e) {
    std::cerr << "Some other error happened during sort, lower_bound, or upper_bound:" << e.what() << std::endl;
    exit(-1);
  }
 
  return ret;
}

void Component::CudaSortBySequence(Component::cuSharedStream cr)
{
  thrust::device_vector<cudaParticle>::iterator
    pbeg = cr->cuda_particles.begin(),
    pend = cr->cuda_particles.end();

  thrust::sort(thrust::cuda::par.on(cr->stream), pbeg, pend, LessCudaSeq());
}

void Component::ParticlesToCuda(PartMap::iterator beg, PartMap::iterator fin)
{
  if (step_timing and use_cuda) comp->timer_cuda.start();

  auto npart = std::distance(beg, fin);
  
  if (host_particles.capacity()<npart) host_particles.reserve(npart);
  host_particles.resize(npart);

  hostPartItr hit = host_particles.begin();
  for (auto pit=beg; pit!=fin; pit++) {
    ParticleHtoD(pit->second, *(hit++));
  }

  if (step_timing and use_cuda) comp->timer_cuda.stop();
}

void Component::HostToDev(Component::cuSharedStream cr)
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

void Component::DevToHost(Component::cuSharedStream cr)
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
  if (step_timing and use_cuda) comp->timer_cuda.start();

  // DEBUG PRINTING
  const int imax = 0;
  int icnt = 0;

  for (hostPartItr v=beg; v!=end; v++) {
    cudaParticle & p = *v;
    if (icnt < imax) {
      std::cout << "[" << icnt++ << ", " << myid << "] " << p << std::endl;
    }
    ParticleDtoH(p, particles[p.indx]);
  }

  if (step_timing and use_cuda) comp->timer_cuda.stop();
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
  size_t psize  = particles.size();
  
  if (multistep) {
    std::pair<unsigned int, unsigned int>
      ret = CudaGetLevelRange(minlev, multistep);
    
    thrust::transform(// thrust::cuda::par.on(cuStream->stream),
		      thrust::cuda::par,
		      cuStream->cuda_particles.begin()+ret.first, cuStream->cuda_particles.end(),
		      cuStream->cuda_particles.begin()+ret.first, cudaZeroAcc());
  } else {
    thrust::transform(// thrust::cuda::par.on(cuStream->stream),
		      thrust::cuda::par,
		      cuStream->cuda_particles.begin(), cuStream->cuda_particles.end(),
		      cuStream->cuda_particles.begin(), cudaZeroAcc());
  }
  
}

struct getMass : public thrust::unary_function<cudaParticle, cuFP_t>

{
 __host__ __device__
 cuFP_t operator()(const cudaParticle& p) const
  {
    return p.mass;
  }
};

struct getPos : public thrust::unary_function<cudaParticle, cuFP_t>
{
  int _t;
  getPos(int t) : _t(t) {}

 __host__ __device__
 cuFP_t operator()(const cudaParticle& p) const
  {
    return p.mass * p.pos[_t];
  }
};

struct getVel : public thrust::unary_function<cudaParticle, cuFP_t>
{
  int _t;
  getVel(int t) : _t(t) {}

 __host__ __device__
 cuFP_t operator()(const cudaParticle& p) const
  {
    return p.mass * p.vel[_t];
  }
};

struct getAcc : public thrust::unary_function<cudaParticle, cuFP_t>
{
  int _t;
  getAcc(int t) : _t(t) {}

 __host__ __device__
 cuFP_t operator()(const cudaParticle& p) const
  {
    return p.mass * p.acc[_t];
  }
};

void Component::fix_positions_cuda(unsigned mlevel)
{
				// Zero center
  for (int i=0; i<3; i++) center[i] = 0.0;

  				// Zero variables
  mtot = 0.0;
  for (int k=0; k<dim; k++) com[k] = cov[k] = coa[k] = 0.0;

				// Zero multistep counters at and
				// above this level
  try {
    auto exec = thrust::cuda::par.on(cuStream->stream);
    
    for (unsigned mm=mlevel; mm<=multistep; mm++) {

      cudaParticle trg;
      trg.lev[0] = mm;
      
      thrust::device_vector<cudaParticle>::iterator
	pbeg = cuStream->cuda_particles.begin(),
	pend = cuStream->cuda_particles.end();
    
      thrust::device_vector<cudaParticle>::iterator lo, hi;

      cudaStreamSynchronize(cuStream->stream);

      if (thrust_binary_search_workaround) {
	cudaStreamSynchronize(cuStream->stream);
	lo  = thrust::lower_bound(pbeg, pend, trg, LessCudaLev());
      } else {
	lo = thrust::lower_bound(exec, pbeg, pend, trg, LessCudaLev());
      }
      
      cudaStreamSynchronize(cuStream->stream);

      if (thrust_binary_search_workaround) {
	hi = thrust::upper_bound(pbeg, pend, trg, LessCudaLev());
      } else {
	hi = thrust::upper_bound(exec, pbeg, pend, trg, LessCudaLev());
      }
      
      com_mas[mm] = thrust::transform_reduce(lo, hi, getMass(), 0.0, thrust::plus<cuFP_t>());
      for (unsigned k=0; k<3; k++)  {
	com_lev[3*mm+k] = thrust::transform_reduce(lo, hi, getPos(k), 0.0, thrust::plus<cuFP_t>());
	cov_lev[3*mm+k] = thrust::transform_reduce(lo, hi, getVel(k), 0.0, thrust::plus<cuFP_t>());
	coa_lev[3*mm+k] = thrust::transform_reduce(lo, hi, getAcc(k), 0.0, thrust::plus<cuFP_t>());
      }
    }
  }
  catch(thrust::system_error &e) {
    std::cerr << "Some other error happened during sort, lower_bound, or upper_bound:" << e.what() << std::endl;
    exit(-1);
  }
 
  std::vector<double> com1(3, 0.0), cov1(3, 0.0), coa1(3, 0.0);
  double              mtot1 = 0.0;

  for (unsigned mm=0; mm<=multistep; mm++) {
    for (int k=0; k<3; k++) {
      com1[k] += com_lev[3*mm + k];
      cov1[k] += cov_lev[3*mm + k];
      coa1[k] += coa_lev[3*mm + k];
    }
    mtot1 += com_mas[mm];
  }

  MPI_Allreduce(&mtot1, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&com1[0], com, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&cov1[0], cov, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&coa1[0], coa, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
  if (VERBOSE>5) {
				// Check for NaN
    bool com_nan = false, cov_nan = false, coa_nan = false;
    for (int k=0; k<3; k++)
      if (std::isnan(com[k])) com_nan = true;
    for (int k=0; k<3; k++)
      if (std::isnan(cov[k])) cov_nan = true;
    for (int k=0; k<3; k++)
      if (std::isnan(coa[k])) coa_nan = true;
    if (com_nan && myid==0)
      cerr << "Component [" << name << "] com has a NaN" << endl;
    if (cov_nan && myid==0)
      cerr << "Component [" << name << "] cov has a NaN" << endl;
    if (coa_nan && myid==0)
      cerr << "Component [" << name << "] coa has a NaN" << endl;
  }
				// Compute component center of mass and
				// center of velocity, and center of accel

  if (mtot > 0.0) {
    for (int k=0; k<dim; k++) com[k]  /= mtot;
    for (int k=0; k<dim; k++) cov[k]  /= mtot;
    for (int k=0; k<dim; k++) coa[k]  /= mtot;
  }

  if (com_system and not consp) {
    for (int k=0; k<dim; k++) com0[k] = com[k];
    for (int k=0; k<dim; k++) cov0[k] = cov[k];
  }

  if (com_system) {	   // Use local center of accel for com update
    for (int k=0; k<dim; k++) acc0[k]  = coa[k];
  } else {			// No mass, no acceleration?
    for (int k=0; k<dim; k++) acc0[k]  = 0.0;
  }

  if ((EJ & Orient::CENTER) && !EJdryrun) {
    Vector ctr = orient->currentCenter();
    bool ok    = true;
    for (int i=0; i<3; i++) {
      if (std::isnan(ctr[i+1])) ok = false;
    } 
    if (ok) {
      for (int i=0; i<3; i++) center[i] += ctr[i+1];
    } else if (myid==0) {
      cout << "Orient: center failure, T=" << tnow 
	   << ", adjustment skipped" << endl;
    }
  }
}

struct testCountLevel :  public thrust::unary_function<cudaParticle, int>
{
  int _l;

  testCountLevel(int l) : _l(l) {}

  __host__ __device__
  int operator()(const cudaParticle& p) const
  {
    if (p.lev[0] == _l) return 1;
    return 0;
  }
};


void Component::print_level_lists_cuda(double T)
{
				// Retrieve counts per level
  std::vector<int> cntr(multistep+1);
  for (int m=0; m<=multistep; m++) {
    cntr[m] = thrust::transform_reduce(cuStream->cuda_particles.begin(),
				       cuStream->cuda_particles.end(),
				       testCountLevel(m),
				       0, thrust::plus<int>());
  }

  if (myid==0) {
				// Sum reduce to root
    MPI_Reduce(MPI_IN_PLACE, &cntr[0], multistep+1, MPI_INT, MPI_SUM,
	       0, MPI_COMM_WORLD);

    int tot=0;
    for (int n=0; n<=multistep; n++) tot += cntr[n];

    if (tot) {

      std::ostringstream ofil;
      ofil << runtag << ".levels";
      std::ofstream out(ofil.str().c_str(), ios::app);

      int sum=0;
      out << setw(60) << setfill('-') << '-' << endl;
      std::ostringstream sout;
      sout << "--- Component <" << name 
	   << ", " << id  << ">, T=" << T;
      out << std::setw(60) << std::left << sout.str().c_str() << std::endl;
      out << std::setw(60) << '-' << std::endl << std::setfill(' ');
      out << std::setw(3)  << "L" 
	  << std::setw(10) << "Number" 
	  << std::setw(10) << "dN/dL" 
	  << std::setw(10) << "N(<=L)"
	  << std::endl
	  << std::setw(60) << std::setfill('-') << '-'
	  << std::endl << std::setfill(' ');
      for (int n=0; n<=multistep; n++) {
	sum += cntr[n];
	out << std::setw(3)  << n 
	    << std::setw(10) << cntr[n] << std::setprecision(3) << std::fixed
	    << std::setw(10) << static_cast<double>(cntr[n])/tot
	    << std::setw(10) << static_cast<double>(sum)    /tot;
	out << std::endl;
      }
      out << std::endl << std::setw(3) << "T" << std::setw(10) << tot
	  << std::endl << std::endl << std::right;
    } else {
      std::cout << "print_level_lists [" << name 
		<< ", T=" << tnow << "]: tot=" << tot << std::endl;
    }

  } else {
				// Sum reduce counts to root
    MPI_Reduce(&cntr[0], 0, multistep+1, MPI_INT, MPI_SUM,
	       0, MPI_COMM_WORLD);
  }

}

// -*- C++ -*-
