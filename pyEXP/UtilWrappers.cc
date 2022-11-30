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
	"Compute the center of the particle component using the density\n"
	"weighted position using KD NN estimator.  A stride >1 will\n"
	"generate a subsample of every nth particle over a random\n"
	"permutation. Ndens is the number of particles per sample ball\n"
	"(32 is a good compromise between accuracy and runtime; 16 is\n"
	"okay if you are trying to shave off runtime. Nsort>0 keeps\n"
	"the particles of the Nsort densest samples.",
	py::arg("reader"), py::arg("stride")=1,
	py::arg("Nsort")=0, py::arg("Ndens")=32);

  m.def("getCenterOfMass", &getCenterOfMass,
	"Compute the center of mass for the particle component",
	py::arg("reader"));

  m.def("particleIterator", &particleIterator,
	"Apply a user-defined functor to every particle in phase space",
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
	"Report on the version and git commit.  This is the same version\n"
	"information reported by the EXP N-body code.");
}
