#include "KDdensity.H"

namespace Utility
{
  KDdensity::KDdensity(PR::PRptr reader, int Ndens) : Ndens(Ndens)
  {
    // Check for MPI and initialize if needed
    //
    int flag;
    MPI_Initialized(&flag);
    bool mpi_enabled = (flag != 0);

    // Get number of particles
    //
    int nbod = reader->CurrentNumber();

    // Every node needs to make the tree (a parallel share could be
    // implemented in KDtree.H).  Build the ID map.
    //
    KDmass = 0.0;
    for (auto part=reader->firstParticle(); part!=0; part=reader->nextParticle()) {
      KDmass += part->mass;
      points.push_back(point3({part->pos[0], part->pos[1], part->pos[2]}, part->mass));
      indx[part->indx] = points.size()-1;
    }
    
    // Build the k-d tree
    //
    kdtree_ = std::make_shared<tree3>(points.begin(), points.end());


    // Share the density computation among the nodes
    //
    KDdens.resize(nbod, 0.0);
    unsigned long badVol = 0;
    for (int k=0; k<points.size(); k++) {
      if (k % numprocs == myid) {
	auto ret = kdtree_->nearestN(points[k], Ndens);
	double volume = 4.0*M_PI/3.0*std::pow(std::get<2>(ret), 3.0);
	if (volume>0.0 and KDmass>0.0)
	  KDdens[k] = std::get<1>(ret)/volume/KDmass;
	else badVol++;
      }
    }

    if (mpi_enabled)
      MPI_Allreduce(MPI_IN_PLACE, KDdens.data(), nbod,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (myid==0) {
      std::cout << "---- KDdensity: finished density estimate with " << badVol
		<< " undetermined densities, mass=" << KDmass << std::endl;
      std::cout << "---- KDdensity: a few densities are " << KDdens.front()
		<< ", " << KDdens[nbod/2] << ", " << KDdens.back() << std::endl;
      std::cout << "---- KDdensity: bad density count is " << badVol
		<< " of " << nbod << " particles." << std::endl;
    }
  }

  KDdensity::KDdensity(const Eigen::VectorXd& m,
		       const KDdensity::RowMatrixXd& pos,
		       int Ndens) : Ndens(Ndens)
  {
    // Check for MPI and initialize if needed
    //
    int flag;
    MPI_Initialized(&flag);
    bool mpi_enabled = (flag != 0);

    // Get number of particles
    //
    int nbod = m.size();

    // Sanity check
    //
    if (pos.rows() != nbod or pos.cols() != 3)
      throw std::runtime_error("KDdensity: position array has wrong dimensions");

    // Every node needs to make the tree (a parallel share could be
    // implemented in KDtree.H).  Build the ID map.
    //
    KDmass = 0.0;
    for (int i=0; i<nbod; i++) {
      KDmass += m(i);
      points.push_back(point3({pos(i, 0), pos(i, 1), pos(i, 2)}, m(i)));
      indx[i] = i;
    }
    
    // Build the k-d tree
    //
    kdtree_ = std::make_shared<tree3>(points.begin(), points.end());


    // Share the density computation among the nodes
    //
    KDdens.resize(nbod, 0.0);
    unsigned long badVol = 0;
    for (int k=0; k<points.size(); k++) {
      if (k % numprocs == myid) {
	auto ret = kdtree_->nearestN(points[k], Ndens);
	double volume = 4.0*M_PI/3.0*std::pow(std::get<2>(ret), 3.0);
	if (volume>0.0 and KDmass>0.0)
	  KDdens[k] = std::get<1>(ret)/volume/KDmass;
	else badVol++;
      }
    }

    if (mpi_enabled)
      MPI_Allreduce(MPI_IN_PLACE, KDdens.data(), nbod,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (myid==0) {
      std::cout << "---- KDdensity: finished density estimate with " << badVol
		<< " undetermined densities, mass=" << KDmass << std::endl;
      std::cout << "---- KDdensity: a few densities are " << KDdens.front()
		<< ", " << KDdens[nbod/2] << ", " << KDdens.back() << std::endl;
      std::cout << "---- KDdensity: bad density count is " << badVol
		<< " of " << nbod << " particles." << std::endl;
    }
  }

  double KDdensity::getDensity(unsigned long index)
  {
    return KDdens[indx[index]];
  }
}
