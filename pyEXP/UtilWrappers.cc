#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <Centering.H>
#include <ParticleIterator.H>

void UtilityClasses(py::module &m) {

  m.doc() = "Utility class bindings\n\n"
    "This module provides routines for BFE tasks that do not naturally fit\n"
    "into the main categories.  The current list of utilities is:\n"
    "  1. Report the current EXP version, GIT branch and commit\n"
    "  2. Compute the center of the particle distribution from its center\n"
    "     of mass.  Very fast but easily biased.\n"
    "  3. Compute the mean density weighted center of the particle distribu-\n"
    "     tion from KD density estimator at each particle position. Very is\n"
    "     very slow.  One can change the stride to decrease the sample size\n"
    "     to speed this up.\n\n"
    "     Note on centers: the EXP n-body code does this automatically.  The\n"
    "     density weighted center is an alternative for snapshots without\n"
    "     center estimates.  Only use COM if know your simulation remains\n"
    "     close to bisymmetric.\n\n"
    "  4. Apply a user-defined Python function to all particles in a phase-\n"
    "     space reader.  This may be used to do calculations using all or a\n"
    "     user-determined subset of snapshot particles.  The functor has no\n"
    "     return type; it is up to the user to put accumulated values in\n"
    "     the scope.  The functor needs the arguments of mass, position\n"
    "     vector, velocity vector, and index.  For example, the following\n"
    "     Python code computes the center of mass:\n"
    "     #---------------------------------------------------------------\n"
    "     # Variables in the scope of myFunctor\n"
    "     #\n"
    "     totalMass = 0.0\n"
    "     centerOfMass = [0.0, 0.0, 0.0]\n"
    "     #\n"
    "     # Definition of the functor\n"
    "     #\n"
    "     def myFunctor(m, pos, vel, index):\n"
    "        global totalMass, centerOfMass\n"
    "        totalMass += m\n"
    "        for i in range(3): centerOfMass[i] += m*pos[i]\n\n"
    "     #\n"
    "     # Now iterate through the particles provided by a reader instance\n"
    "     #\n"
    "     pyEXP.util.particleIterator(reader, myFunctor)\n"
    "     #\n"
    "     # Print the COM\n"
    "     #\n"
    "     for i in range(3): centerOfMass[i] /= totalMass\n"
    "     #\n"
    "     #---------------------------------------------------------------\n\n";

  using namespace Utility;

  m.def("getDensityCenter", &getDensityCenter,
	R"(
        Compute the center of the particle component

        This implementation uses the density weighted position using KD 
        N-nearest neighbor estimator.  

        Parameters
        ----------
        reader : ParticleReader
            the particle-reader class instance
        stride : int, default=1
             stride >1 will generate a subsample of every nth particle 
             over a random permutation
        Ndens : int, default=32
             number of particles per sample ball (32 is a good 
             compromise between accuracy and runtime; 16 is okay if 
             you are trying to shave off runtime. 
        Nsort : int, default=0
             Nsort >0 keeps the particles of the Nsort densest samples

        Returns
        -------
        list(float)
            Computed center
        )",
	py::arg("reader"), py::arg("stride")=1,
	py::arg("Nsort")=0, py::arg("Ndens")=32);

  m.def("getCenterOfMass", &getCenterOfMass,
	R"(
        Compute the center of mass for the particle component

        Parameters
        ----------
        reader : ParticleReader
            the particle-reader class instance

        Returns
        -------
        list(float)
            Computed center
        )", py::arg("reader"));

  m.def("particleIterator", &particleIterator,
	R"(
        Apply a user-defined functor to every particle in phase space

        Parameters
        ----------
        reader : ParticleReader
            the particle-reader class instance
        functor :
            the callback function to compute a phase-space quantity

        Returns
        -------
        None

        Notes
        -----
        The callback function must have the signature:

        void(float, list, list, int)

        where the first argument is mass, the second is position,
        the third is velocity, and the fourth is index.  Not all
        values need to be used in the function, of course.
        )",
	py::arg("reader"), py::arg("functor"));

  m.def("getVersionInfo",
	[]() {
	  const int W = 80;		// Full linewidth
	  std::ostringstream sout;	// Get a std::string from the string
	  // literal
	  sout << "%%%%% This is " << PACKAGE_STRING << " ";
				// Print the info block
	  std::cout << std::endl
		    << std::setw(W) << std::setfill('%') << '%' << std::endl
		    << std::left << setw(W) << sout.str() << std::endl
		    << std::setw(W) << std::setfill('%') << '%' << std::endl
		    << std::setfill(' ')
		    << std::setw(20) << "%%%%% Repository URL" << " | "
		    << std::setw(W-24) << PACKAGE_URL << '%' << std::endl
		    << std::setw(20) << "%%%%% Current branch" << " | "
		    << std::setw(W-24) << GIT_BRANCH << '%' << std::endl
		    << std::setw(20) << "%%%%% Current commit" << " | "
		    << std::setw(W-24) << GIT_COMMIT << '%' << std::endl
		    << std::setw(20) << "%%%%% Compile time"   << " | "
		    << std::setw(W-24) << COMPILE_TIME << '%' << std::endl
		    << std::setfill('%')
		    << std::setw(W) << '%' << std::setfill(' ') << std::endl
		    << std::endl;
	},
	R"(
        Report on the version and git commit.  

        This is the same version information reported by the EXP N-body code.)");
  m.def("Version",
    []() {
      std::string version_str = VERSION;
      std::istringstream iss(version_str);
      std::string token;
      int parts[3] = {0, 0, 0};
      int i = 0;

      while (std::getline(iss, token, '.') && i < 3) {
          parts[i++] = std::stoi(token);
      }

      return std::make_tuple(parts[0], parts[1], parts[2]);
  },
    R"(
          Return the version.  
  
          This is the same version information reported by the EXP N-body code.)");

  m.def("setMPI",
    [](bool verbose) {
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
      MPI_Get_processor_name(processor_name, &proc_namelen);

      workers = numprocs - 1;

      if (workers) {
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	std::vector<int> nworkers (workers);
	
	for (int n=1; n<numprocs; n++) nworkers[n-1] = n;
	MPI_Group_incl(world_group, workers, &nworkers[0], &worker_group);
	MPI_Comm_create(MPI_COMM_WORLD, worker_group, &MPI_COMM_WORKER);
      }

      if (verbose) {
	if (myid==0)
	  std::cout << std::string(80, '%') << std::endl
		    << std::setfill('%') << std::setw(80) << std::left
		    << "%%%%% Node, process, and communicator assignment " 
		    << std::endl << std::string(80, '%') << std::endl
		    << std::setfill(' ')
		    << std::right << std::setw(8)  << "Node #"
		    << " | " << std::setw(20) << "Hostname"
		    << " | " << std::setw(8)  << "PID"
		    << " | " << std::setw(10) << "Status"
		    << std::endl << std::string(80, '%') << std::endl << std::left;
	
	for (int n=0; n<numprocs; n++) {

	  if (n == myid) {

	    std::cout << std::right << std::setw(8)  << myid 
		      << " | " << std::setw(20) << processor_name 
		      << " | " << std::setw(8)  << getpid();
	  
	    if (n) {
	      int m; MPI_Group_rank ( worker_group, &m );
	      std::cout << " | " << std::setw(10) << "WORKER " << m << std::endl;
	    } else {
	      std::cout << " | " << std::setw(10) << "ROOT" << std::endl;
	    }
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	}
      }
    },
    R"(
      Initialize EXP's C++ MPI communicators and variables from pyEXP
  
      Sets all of the internal MPI-specific library values necessary
      to use the full complement of the EXP C++ classes from Python
      scripts.  For example, this is needed for ParticleReader, which
      does not provide its own MPI initialization. The BasisFactory,
      Basis, ParticleReader, and FieldGenerator classes can do this
      automatically.

      Parameters
      ----------
      verbose : bool
         Print the node and process assignments.  Default: false.

      Returns
      -------
      None

   )", py::arg("verbose") = false);
}
