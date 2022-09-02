#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

extern void MSSAtoolkitClasses   (py::module &m);
extern void CoefFactoryClasses   (py::module &m);
extern void BasisFactoryClasses  (py::module &m);
extern void FieldGeneratorClasses(py::module &m);
extern void ParticleReaderClasses(py::module &m);

PYBIND11_MODULE(pyEXP, m)
{
  m.doc() = "A collection of EXP tools for processing and analysing simulation\ndata using BFE techniques.  Run 'help' on the submodules below for\nmore detailed information...";

  auto mod_coefs = m.def_submodule("coefs",
				   "Classes for reading, passing and "
				   "manipulating coefficient sets");

  auto mod_basis = m.def_submodule("basis",
				   "Create and apply specific biorthogonal "
				   "bases to generate coefficients from "
				   "particle data and evaluate potential, "
				   "density, and force fields");

  auto mod_field = m.def_submodule("field",
				   "Create two- and three-dimension rectangular "
				   "grids of fields for visualization");

  auto mod_read = m.def_submodule("read", "Read particle snapshots of various "
				  "types.  Currently EXP, Gadget, and Tipsy "
				  "types are supported.");

  auto mod_mssa = m.def_submodule("mssa",  "Tools to apply Multivariate Singular "
				  "Spectrum Analysis (MSSA) to the coefficients "
				  "computed using the 'basis' classes");
  
  CoefFactoryClasses(mod_coefs);
  BasisFactoryClasses(mod_basis);
  FieldGeneratorClasses(mod_field);
  ParticleReaderClasses(mod_read);
  MSSAtoolkitClasses(mod_mssa);
}

