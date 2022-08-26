#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <FieldGenerator.H>

namespace py = pybind11;

void FieldGeneratorClasses(py::module &m) {

  m.doc() = "FieldGenerator class bindings";

  using namespace Field;

  py::class_<Field::FieldGenerator>(m, "FieldGenerator")
    .def(py::init<const std::vector<double>, const std::vector<double>,
	 const std::vector<double>, const std::vector<int>>())
    .def("slices", &Field::FieldGenerator::slices)
    .def("file_slices", &Field::FieldGenerator::file_slices)
    .def("volumes", &Field::FieldGenerator::volumes)
    .def("file_volumes", &Field::FieldGenerator::file_volumes);
}

