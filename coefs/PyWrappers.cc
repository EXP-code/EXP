#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

extern void CoefFactoryClasses   (py::module &m);
extern void BasisFactoryClasses  (py::module &m);
extern void FieldGeneratorClasses(py::module &m);
extern void ParticleReaderClasses(py::module &m);

PYBIND11_MODULE(pyEXP, m)
{
  m.doc() = "A collection of EXP tools for processing and analysing simulation data using BFE techniques";

  auto mcoefs = m.def_submodule("coefs");
  auto mbasis = m.def_submodule("basis");
  auto mfield = m.def_submodule("field");
  auto mpread = m.def_submodule("pread");

  mcoefs.doc() = "Coefficient containers and manipulators for all EXP geometries (spherical, cylindrical, cubic, etc.)";

  mbasis.doc() = "Create and apply specific biorthogonal bases to generate coefficients from particle data and evaluate potential, density, and force fields";

  mfield.doc() = "Create two- and three-dimension rectangular grids of fields for visualization";

  mpread.doc() = "Read particle snapshots of various types.  Currently EXP, Gadget, and Tipsy types are supported.";

  CoefFactoryClasses(mcoefs);
  BasisFactoryClasses(mbasis);
  FieldGeneratorClasses(mfield);
  ParticleReaderClasses(mpread);
  
}

