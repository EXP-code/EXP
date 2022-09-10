				// C++/STL headers
#include <filesystem>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

#include <ParticleReader.H>
#include <localmpi.H>
#include <KDtree.H>

namespace Utility
{
  std::vector<double> getDensityCenter(PR::PRptr reader, int Ndens)
  {
    bool use_mpi;
    int flag;
    MPI_Initialized(&flag);
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

    for (int i=0; i<nbods; i++) {
      if (i % numprocs == myid) {
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

