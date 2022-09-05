#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <FieldGenerator.H>

namespace py = pybind11;

//! Helper function that maps the Eigen::Tensor<T, 3> into an numpy.ndarray
template <typename T>
py::array_t<T> make_ndarray(Eigen::Tensor<T, 3>& mat)
{
  // Get the tensor dimenions
  auto dims = mat.dimensions();

  // Check rank
  if (dims.size() != 3) {
    std::ostringstream sout;
    sout << "make_ndarray: tensor rank must be 3, found " << dims.size();
    throw std::runtime_error(sout.str());
  }
  
  // Make the memory mapping
  return py::array_t<T>
    (
     // shape
     {dims[0], dims[1], dims[2]},
     // C-style contiguous strides for double
     {sizeof(T), dims[0]*sizeof(T), dims[0]*dims[1]*sizeof(T)},
     // the data pointer
     mat.data()
     );
}

void FieldGeneratorClasses(py::module &m) {

  m.doc() = "FieldGenerator class bindings.  This class computes surfaces\n"
    "and volumes for visualizing the physical quantities implied\nby your"
    "coefficients.\n\n"
    "The generator is constructed by passing a vector of desired times\n"
    "that must be in the coefficient object and list of lower bounds,\n"
    "a list of upper bounds, and a list of knots per dimension.  These\n"
    "lists all have rank 3 for (x, y, z).  For a two-dimensional surface,\n"
    "one of the knot array values must be zero.  The member functions\n"
    "slices and volumes, called with the basis and coefficient objects,\n"
    "return a numpy.ndarray containing the field evaluations.  See\n",
    "help(pyEXP.basis) and help(pyEXP.coefs) for info on the basis\n"
    "and coefficient objects.\n\n";

  using namespace Field;

  py::class_<Field::FieldGenerator, std::shared_ptr<Field::FieldGenerator>>
    f(m, "FieldGenerator");

  f.def(py::init<const std::vector<double>, const std::vector<double>,
	const std::vector<double>, const std::vector<int>>(),
	"Create fields for given times and lower and upper bounds and "
	"grid sizes", py::arg("times"), py::arg("lower"), py::arg("upper"),
	py::arg("gridsize"));

  f.def("slices", &Field::FieldGenerator::slices,
	"Return a dictionary of grids (2d numpy arrays) indexed by "
	"time and field type", py::arg("basis"), py::arg("coefs"));
  
  f.def("file_slices", &Field::FieldGenerator::file_slices,
	"Write 2d field grids to files using the supplied string prefix",
	py::arg("basis"), py::arg("coefs"), py::arg("filename"),
	py::arg("dir")="");

  f.def("volumes", [](FieldGenerator& A,
		      Basis::BasisPtr basis, Coefs::CoefsPtr coefs)
  {
    std::map<double, std::map<std::string, py::array_t<float>>> ret;
    auto vols = A.volumes(basis, coefs);
    for (auto & v : vols) {
      for (auto & u : v.second) {
	ret[v.first][u.first] = make_ndarray<float>(u.second);
      }
    }

    return ret;
  },
    "Returns a dictionary of volume grids (3d numpy arrays) indexed by "
    "time and field type");

  f.def("file_volumes", &Field::FieldGenerator::file_volumes,
	"Write 3d field grids to files using the supplied string prefix");
}
