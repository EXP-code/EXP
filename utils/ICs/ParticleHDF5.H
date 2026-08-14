#ifndef EXP_PARTICLE_HDF5_H
#define EXP_PARTICLE_HDF5_H

// Shared writer for the simple EXP phase-space HDF5 schema.  The
// schema unchanged from the Component update and is the same one tested 
// by the hdf5bods and slabics implementations.

#include <algorithm>
#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <H5Zpublic.h>
#include <mpi.h>

namespace EXP
{

  namespace ParticleHDF5
  {

    //! Collects the mass, position, and velocity of a particle container into a vector of 7-element arrays.
    template <typename Particles>
    inline std::vector<std::array<double, 7>> records(const Particles& particles)
    {
      std::vector<std::array<double, 7>> result;
      result.reserve(particles.size());
      for (const auto& p : particles)
	result.push_back({p.mass, p.pos[0], p.pos[1], p.pos[2],
	    p.vel[0], p.vel[1], p.vel[2]});
      return result;
    }

    //! Create a HighFive::DataSetCreateProps object with chunking and optional filter.
    inline HighFive::DataSetCreateProps properties(hsize_t size, unsigned filter)
    {
      if (size == 0) throw std::invalid_argument("cannot write an empty particle file");

      HighFive::DataSetCreateProps props;
      // chunk size is min(size, 256kB) but at least 1kB
      props.add(HighFive::Chunking({std::min<hsize_t>(size, std::min<hsize_t>(262144, std::max<hsize_t>(1024, size)))}));
      if (filter != 3 && filter != 4 && filter != 32001) props.add(HighFive::Shuffle());
      
      // The filter IDs are from HDF5 1.12.0 and later, see https://support.hdfgroup.org/HDF5/doc/Advanced/Filters.html
      // Filter types:
      // 1: GZIP (DEFLATE)
      // 2: SZIP
      // 3: SHUFFLE
      // 4: FLETCHER32
      // 307: LZ4
      // 32004: ZSTD
      // 32001: BLOSC
      //
      // The filter options are specific to each filter type, see the HDF5 documentation for details.
      //
      hid_t plist = props.getId();
      switch (filter) {
      case 1:
      case 307: { const unsigned level = 5; H5Pset_filter(plist, filter, H5Z_FLAG_OPTIONAL, 1, &level); break; }
      case 2: { const unsigned options[] = {141, 16}; H5Pset_filter(plist, filter, H5Z_FLAG_OPTIONAL, 2, options); break; }
      case 3: return props;                 // shuffle only
      case 4: H5Pset_filter(plist, filter, H5Z_FLAG_OPTIONAL, 0, nullptr); return props;
      case 32004: { const unsigned options[] = {0}; H5Pset_filter(plist, filter, H5Z_FLAG_OPTIONAL, 1, options); break; }
      case 32001: { const unsigned options[] = {0, 0, 0, 0, 5, 1, 1}; H5Pset_filter(plist, filter, H5Z_FLAG_OPTIONAL, 7, options); break; }
      default:
	throw std::invalid_argument("unsupported HDF5 filter ID");
      }
      return props;
    }

    //! Write a particle phase-space file with optional auxiliary integer and floating-point data.
    template <typename T>
    inline void write(const std::string& filename,
		      const std::vector<T>& m, const std::vector<T>& x,
		      const std::vector<T>& y, const std::vector<T>& z,
		      const std::vector<T>& u, const std::vector<T>& v,
		      const std::vector<T>& w,
		      const std::vector<std::vector<int>>& aux_ints,
		      const std::vector<std::vector<T>>& aux_floats,
		      const std::optional<std::vector<unsigned long>>& index,
		      unsigned filter)
    {
      const auto count = m.size();
      if (count == 0) return;
      if (x.size() != count || y.size() != count || z.size() != count ||
	  u.size() != count || v.size() != count || w.size() != count)
	throw std::invalid_argument("inconsistent phase-space vector sizes");
      
      HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

      const int n = static_cast<int>(count);
      file.createAttribute<int>("num_particles", HighFive::DataSpace::From(n)).write(n);

      const int num_aux_ints = static_cast<int>(aux_ints.size());
      const int num_aux_floats = static_cast<int>(aux_floats.size());
      file.createAttribute<int>("num_aux_ints", HighFive::DataSpace::From(num_aux_ints)).write(num_aux_ints);
      file.createAttribute<int>("num_aux_floats", HighFive::DataSpace::From(num_aux_floats)).write(num_aux_floats);

      auto group = file.createGroup("particles");
      auto props = properties(count, filter);
      group.createDataSet("m", m, props); group.createDataSet("x", x, props); group.createDataSet("y", y, props);
      group.createDataSet("z", z, props); group.createDataSet("u", u, props); group.createDataSet("v", v, props);
      group.createDataSet("w", w, props);
      if (index) group.createDataSet("index", *index, props);
      for (size_t j = 0; j < aux_ints.size(); ++j)
	group.createDataSet("aux_int_" + std::to_string(j), aux_ints[j], props);
      for (size_t j = 0; j < aux_floats.size(); ++j)
	group.createDataSet("aux_float_" + std::to_string(j), aux_floats[j], props);
    }

    //! Write a particle phase-space file from a vector of 7-element arrays with optional auxiliary integer and floating-point data.
    template <typename T>
    inline void write(const std::string& filename,
		      const std::vector<std::array<double, 7>>& particles,
		      int num_aux_ints, int num_aux_floats, unsigned filter)
    {
      std::vector<T> m, x, y, z, u, v, w;
      m.reserve(particles.size()); x.reserve(particles.size()); y.reserve(particles.size()); z.reserve(particles.size());
      u.reserve(particles.size()); v.reserve(particles.size()); w.reserve(particles.size());
      for (const auto& p : particles) {
	m.push_back(p[0]); x.push_back(p[1]); y.push_back(p[2]); z.push_back(p[3]);
	u.push_back(p[4]); v.push_back(p[5]); w.push_back(p[6]);
      }
      std::vector<std::vector<int>> auxi(num_aux_ints, std::vector<int>(particles.size(), 0));
      std::vector<std::vector<T>> auxf(num_aux_floats, std::vector<T>(particles.size(), 0));
      write(filename, m, x, y, z, u, v, w, auxi, auxf, std::nullopt, filter);
    }

    //! Collect rank-local phase space records and have rank zero write one file.
    inline void gather_and_write(const std::string& filename,
				 const std::vector<std::array<double, 7>>& local,
				 int num_aux_ints, int num_aux_floats,
				 unsigned filter, bool double_precision,
				 MPI_Comm comm = MPI_COMM_WORLD)
    {
      // Check if MPI is initialized
      //
      int is_init = 0;
      MPI_Initialized(&is_init);

      // Gather the local counts to rank 0 only if MPI is initialized, otherwise just use the local data
      //
      if (is_init) {
	int myid=0, ranks=1, local_count = static_cast<int>(local.size());
	// Get the rank and size of the communicator
	MPI_Comm_rank(comm, &myid );
	MPI_Comm_size(comm, &ranks);

	// Gather the local counts to rank 0
	std::vector<int> counts(myid == 0 ? ranks : 0);

	MPI_Gather(&local_count, 1, MPI_INT,
		   myid == 0 ? counts.data() : nullptr,
		   1, MPI_INT, 0, comm);

	// Gather the local data to rank 0
	std::vector<double> send(7 * local.size());
	for (size_t i = 0; i < local.size(); ++i)
	  std::copy(local[i].begin(), local[i].end(), send.begin() + 7*i);

	std::vector<int> counts7, offsets;
	std::vector<double> received;

	// Only Rank 0 needs to know the counts and offsets for the gather
	//
	if (myid == 0) {
	  counts7.resize(ranks);
	  offsets.resize(ranks);
	  int total = 0;
	  for (int i = 0; i < ranks; ++i) {
	    offsets[i] = total;
	    counts7[i] = 7 * counts[i];
	    total += counts7[i];
	  }
	  received.resize(total);
	}

	// Do the gather of the 7-element arrays to Rank 0
	//
	MPI_Gatherv(send.data(), static_cast<int>(send.size()), MPI_DOUBLE,
		    myid == 0 ? received.data() : nullptr,
		    myid == 0 ? counts7. data() : nullptr,
		    myid == 0 ? offsets. data() : nullptr,
		    MPI_DOUBLE, 0, comm);

	if (myid != 0) return;

	// Rank 0 reconstructs the all vector from the received data
	//
	std::vector<std::array<double, 7>> all(received.size()/7);
	for (size_t i = 0; i < all.size(); ++i)
	  std::copy(received.begin() + 7*i, received.begin() + 7*(i+1), all[i].begin());

	// Write the gathered data to the HDF5 file
	//
	if (double_precision) write<double>(filename, all, num_aux_ints, num_aux_floats, filter);
	else                  write<float >(filename, all, num_aux_ints, num_aux_floats, filter);
	
      } else {
	// Write the local/only data to the HDF5 file
	//
	if (double_precision) write<double>(filename, local, num_aux_ints, num_aux_floats, filter);
	else                  write<float >(filename, local, num_aux_ints, num_aux_floats, filter);
      }
    }

  } // namespace ParticleHDF5

} // namespace EXP

#endif
