#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <Centering.H>

void UtilityClasses(py::module &m) {

  using namespace Utility;

  m.doc() = "Utility class bindings\n\n"
    "This module provides routines for BFE tasks that do not naturally fit\n"
    "into the main categories.  The current list of utilities is:\n"
    "  1. Compute the center of the particle distribution from its center\n"
    "     of mass.  Very fast but easily biased.\n"
    "  2. Compute the mean density weighted center of the particle distribu-\n"
    "     tion from KD density estimator at each particle position. Very is\n"
    "     very slow.  One can change the stride to decrease the sample size\n"
    "     to speed this up.\n\n"
    "     Note on centers: the EXP n-body code does this automatically.  The\n"
    "     density weighted center is an alternative for snapshots without\n"
    "     center estimates.  Only use COM if know your simulation remains\n"
    "     close to bisymmetric.\n\n";

  using namespace Utility;

  m.def("getDensityCenter", &getDensityCenter,
	"Compute the center of the particle component using the density "
	"weighted position using KD NN estimator.  Ndens is the number "
	"of particles per sample ball (32 is a good choice; 16 is okay if "
	"you are trying to shave off runtime. A stride >1 will generate a "
	"subsample of every nth particle over a random permutation",
	py::arg("reader"), py::arg("Ndens")=32, py::arg("stride")=1);

  m.def("getCenterOfMass", &getCenterOfMass,
	"Compute the center of mass for the particle component",
	py::arg("reader"));
}
