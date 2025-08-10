#include <ParticleReader.H>
#include <permutation.H>
#include <localmpi.H>
#include <KDtree.H>

namespace Utility
{
  std::vector<double> getDensityCenter(PR::PRptr reader, int stride,
				       int Nsort, int Ndens)
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

    std::vector<double> ctr(3, 0.0);

    using point3 = KDtree::point<double, 3>;
    using tree3  = KDtree::kdtree<double, 3>;

    std::vector<point3> points;

    double KDmass = 0.0, dentot = 0.0;
    
    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {
      KDmass += p->mass;
      points.push_back(point3({p->pos[0], p->pos[1], p->pos[2]}, p->mass));
    }
	
    if (use_mpi) {

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
      // END mpi procs loop

    }
    // END MPI block
	
    tree3 tree(points.begin(), points.end());
  
    int nbods = points.size();

    std::shared_ptr<permutation> sigma;

    // Generate the permutation if the stride is greater than 1.  We
    // need this because the phase space from some codes may have
    // order imposed by their domain decomposition schemes.
    //
    if (stride>1) {
      sigma = std::make_shared<permutation>(nbods);
      sigma->shuffle();
      nbods /= stride;		// Compute the new subsample size
    }

    // Density stack for Nsort>0
    //
    std::multimap<double, int> stack;

    // Run through the bodies or subsampled bodies
    //
    for (int j=0; j<nbods; j++) {
      if (j % numprocs == myid) {
	int i = j;		    // Default to no permutation
	if (sigma) i = (*sigma)[j]; // The permutation if required

	auto ret = tree.nearestN(points[i], Ndens);

	double volume = 4.0*M_PI/3.0*std::pow(std::get<2>(ret), 3.0);
	double density = 0.0;
	if (volume>0.0 and KDmass>0.0) {
	  density = std::get<1>(ret)/volume/KDmass;
	  if (Nsort>0) {
	    stack.insert({density, i});
	    if (stack.size()>Nsort) stack.erase(stack.begin());
	  } else {
	    for (int k=0; k<3; k++) ctr[k] += density * points[i].get(k);
	    dentot += density;
	  }
	}
      }
    }

    if (Nsort>0) {
      auto combo = stack;

      if (use_mpi) {

	for (int n=0; n<numprocs; n++) {
	  if (myid==n) {
	    int sz = stack.size();
	    std::vector<double> dens;
	    std::vector<int>    indx;
	    for (auto v : stack) {
	      dens.push_back(v.first);
	      indx.push_back(v.second);
	    }
	    MPI_Bcast(&sz, 1, MPI_INT, n, MPI_COMM_WORLD);
	    MPI_Bcast(dens.data(), sz, MPI_DOUBLE, n, MPI_COMM_WORLD);
	    MPI_Bcast(indx.data(), sz, MPI_INT,    n, MPI_COMM_WORLD);
	  } else {
	    int sz;
	    MPI_Bcast(&sz, 1, MPI_INT, n, MPI_COMM_WORLD);
	    std::vector<double> dens(sz);
	    std::vector<int>    indx(sz);
	    MPI_Bcast(dens.data(), sz, MPI_DOUBLE, n, MPI_COMM_WORLD);
	    MPI_Bcast(indx.data(), sz, MPI_INT,    n, MPI_COMM_WORLD);
	    for (int i=0; i<sz; i++) combo.insert({dens[i], indx[i]});
	  }
	}
	// END mpi procs loop

      }
      // END MPI block

      int count = 0;
      for (auto it=combo.rbegin(); it!=combo.rend(); it++) {
	if (count++ > Nsort) break;
	for (int k=0; k<3; k++)
	  ctr[k] += it->first * points[it->second].get(k);
	dentot += it->first;
      }
      
    } else {
      if (use_mpi) {
	MPI_Allreduce(MPI_IN_PLACE, ctr.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &dentot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
    }

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
    
    int use_mpi;
    MPI_Initialized(&use_mpi);

    // Fall back sanity for MPI. Appears to allow MPI calls without
    // 'mpirun' but may be an accident.
    //
    if (use_mpi) {
      MPI_Allreduce(MPI_IN_PLACE, ctr.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mastot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    if (mastot>0.0) {
      for (int k=0; k<3; k++) ctr[k] /= mastot;
    }

    return ctr;
  }

}
// END namespace Utility

