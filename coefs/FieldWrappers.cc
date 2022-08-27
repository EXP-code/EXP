#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <FieldGenerator.H>

namespace py = pybind11;

void FieldGeneratorClasses(py::module &m) {

  m.doc() = "FieldGenerator class bindings";

  using namespace Field;

  py::class_<Field::FieldGenerator, std::shared_ptr<Field::FieldGenerator>>
    f(m, "FieldGenerator");

  f.def(py::init<const std::vector<double>, const std::vector<double>,
	const std::vector<double>, const std::vector<int>>());

  f.def("slices", &Field::FieldGenerator::slices,
	"Return a dictionary of grids (2d numpy arrays) indexed by "
	"time and field type");
  
  f.def("file_slices", &Field::FieldGenerator::file_slices,
	"Write 2d field grids to files using the supplied string prefix");

  f.def("volumes", &Field::FieldGenerator::volumes,
	"Returns a dictionary of volume grids (3d numpy arrays) indexed by "
	"time and field type");

  f.def("file_volumes", &Field::FieldGenerator::file_volumes,
	"Write 3d field grids to files using the supplied string prefix");
}

