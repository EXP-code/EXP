#include <iostream>
#include <iomanip>

#include <algorithm>
#include <numeric>
#include <random>

#include <ParticleReader.H>
#include <localmpi.H>
#include <KDtree.H>

/** Make a permutation index.  Share will MPI nodes if MPI is active.

    This will initialze the Mersenne Twister from the random device.
    This initialization is not cryptographically sound, I know, but
    that doesn't matter here.
*/
struct permutation
{
  permutation(unsigned n) : perm(n), g(std::random_device{}())
  {
    if (myid==0) std::iota(perm.begin(), perm.end(), unsigned(0));
  }

  void shuffle() {
    if (myid==0) std::shuffle(perm.begin(), perm.end(), g);
    if (numprocs>1)
      MPI_Bcast(perm.data(), perm.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  }

  //! The permutation operator
  size_t operator[](size_t n) const { return perm[n]; }

private:
  std::vector<size_t> perm;
  std::mt19937 g;
};


namespace Utility
{
  std::vector<double> getDensityCenter(PR::PRptr reader, int Ndens,
				       int stride)
  {
    int flag;
    MPI_Initialized(&flag);

    bool use_mpi;
    if (flag) use_mpi = true;
    else      use_mpi = false;
    
    // Fall back sanity (works for me but this needs to be fixed
    // generally)
    //
    if (use_mpi) {
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    } else {
      int argc = 0; char **argv = 0;
      MPI_Init(&argc, &argv);
    }

    std::vector<double> ctr(3, 0.0);

    typedef point <double, 3> point3;
    typedef kdtree<double, 3> tree3;

    std::vector<point3> points;

    double KDmass = 0.0, dentot = 0.0;
    
    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {
      KDmass += p->mass;
      points.push_back(point3({p->pos[0], p->pos[1], p->pos[2]}, p->mass));
    }
	
    std::vector<double> dd;
    int sz0 = points.size();
    
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	MPI_Bcast(&sz0, 1, MPI_INT, n, MPI_COMM_WORLD);
	dd.resize(sz0*4);
	for (int i=0; i<sz0; i++) {
	  dd[i*4+0] = points[i].get(0);
	  dd[i*4+1] = points[i].get(1);
	  dd[i*4+2] = points[i].get(2);
	  dd[i*4+3] = points[i].mass();
	}
	MPI_Bcast(dd.data(), sz0*4, MPI_DOUBLE, n, MPI_COMM_WORLD);
      } else {
	int sz;
	MPI_Bcast(&sz, 1, MPI_INT, n, MPI_COMM_WORLD);
	dd.resize(sz*4);
	MPI_Bcast(dd.data(), sz*4, MPI_DOUBLE, n, MPI_COMM_WORLD);
	for (int i=0; i<sz; i++) {
	  points.push_back(point3({dd[i*4+0], dd[i*4+1], dd[i*4+2]}, dd[i*4+3]));
	  KDmass += dd[i*4+3];
	}
      }
    }
	
    tree3 tree(points.begin(), points.end());
  
    int nbods = points.size();

    std::shared_ptr<permutation> sigma;

    // Generate the permutation.  We need this because the phase space
    // may have order imposed by various domain decomposition schemes
    // (e.g.).
    //
    if (stride>1) {
      nbods /= stride;
      sigma = std::make_shared<permutation>(nbods);
      sigma->shuffle();
    }

    // Run through the bodies or subsampled bodies
    //
    for (int j=0; j<nbods; j++) {
      if (j % numprocs == myid) {
	int i = j;
	if (sigma) i = (*sigma)[j]; // The permutation if required

	auto ret = tree.nearestN(points[i], Ndens);

	double volume = 4.0*M_PI/3.0*std::pow(std::get<2>(ret), 3.0);
	double density = 0.0;
	if (volume>0.0 and KDmass>0.0)
	  density = std::get<1>(ret)/volume/KDmass;

	for (int k=0; k<3; k++) ctr[k] += density * points[i].get(k);
	dentot += density;
      }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, ctr.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &dentot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (dentot>0.0) {
      for (int k=0; k<3; k++) ctr[k] /= dentot;
    }

    return ctr;
  }

  //! Brute force center of mass computation
  std::vector<double> getCenterOfMass(PR::PRptr reader)
  {
    std::vector<double> ctr(3, 0.0);
    double mastot = 0.0;
  
    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {
      for (int k=0; k<3; k++) ctr[k] += p->mass * p->pos[k];
      mastot += p->mass;
    }
    
    MPI_Allreduce(MPI_IN_PLACE, ctr.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &mastot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if (mastot>0.0) {
      for (int k=0; k<3; k++) ctr[k] /= mastot;
    }

    return ctr;
  }

}
// END namespace Utility

