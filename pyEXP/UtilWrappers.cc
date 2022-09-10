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
    "     very slow.  The EXP n-body code does this automatically.  The\n"
    "     density weighted center is an alternative for snapshots without\n"
    "     center estimates.\n\n";

  using namespace Utility;

  m.def("getDensityCenter", &getDensityCenter,
	"Compute the center of the particle component using the density "
	"weighted position from KD density estimator.  Ndens is the number ",
	"of particles per sample region.",
	py::arg("reader"), py::arg("Ndens")=32);

  m.def("getCenterOfMass", &getCenterOfMass,
	"Compute the center of mass for the particle component",
	py::arg("reader"));
}
