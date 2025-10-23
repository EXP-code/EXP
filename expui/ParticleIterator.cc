#include "ParticleIterator.H"
#include "localmpi.H"

namespace Utility
{
  void particleIterator(PR::PRptr reader, const Callback& func)
  {
    int use_mpi;
    MPI_Initialized(&use_mpi);

    // Fall back sanity for MPI. Appears to allow MPI calls without
    // 'mpirun' but may be an accident.
    //
    if (use_mpi) {
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    } 

    std::vector<double> pp, vv;
    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {
      pp.assign(p->pos, p->pos+3);
      vv.assign(p->vel, p->vel+3);
      func(p->mass, pp, vv, p->indx);
    }
  }

}
// END namespace Utility
