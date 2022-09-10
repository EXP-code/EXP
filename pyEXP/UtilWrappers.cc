#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <Centering.H>

void UtilityClasses(py::module &m) {

  using namespace Utility;

  m.doc() = "Utility class bindings\n\n"
    "This module provides classes that are useful in support of expansion\n"
    "tasks that do not naturally fit into the main categories.  The list\n"
    "uiltiies are:\n"
    "  1. Compute the center of the particle distribution from its center\n"
    "     of mass.  Very fast.\n"
    "  2. Compute the density weighted center of the particle distribution\n"
    "     from KD density estimator at each particle position. Very slow.\n"
    "     While EXP has native methods for doing this, the KD estimator is\n"
    "     an alternative for snapshots without native or prior expansion\n"
    "     centers. The COM is much faster but can be very biased.\n\n";

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
