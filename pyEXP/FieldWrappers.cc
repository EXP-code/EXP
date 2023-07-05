#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <FieldGenerator.H>

namespace py = pybind11;

#include "TensorToArray.H"

void FieldGeneratorClasses(py::module &m) {

  m.doc() = "FieldGenerator class bindings\n\n"
    "FieldGenerator\n"
    "==============\n"
    "This class computes surfaces and volumes for visualizing the physical\n"
    "quantities implied by your basis and associated coefficients.\n\n"
    "The generator is constructed by passing a vector of desired times\n"
    "that must be in the coefficient object and list of lower bounds,\n"
    "a list of upper bounds, and a list of knots per dimension.  These\n"
    "lists all have rank 3 for (x, y, z).  For a two-dimensional surface,\n"
    "one of the knot array values must be zero.  The member functions\n"
    "lines, slices and volumes, called with the basis and coefficient\n"
    "objects, return a numpy.ndarray containing the field evaluations.\n"
    "Each of these functions returns a dictionary of times to a dictionary\n"
    "of field names to numpy.ndarrays at each time.  There are also members\n"
    "which will write these generated fields to files. The linear probe\n"
    "members, 'lines' and 'file_lines', evaluate 'num' field points along\n"
    "a user-specified segment between the 3d points 'beg' and 'end'.  See\n"
    "help(pyEXP.basis) and help(pyEXP.coefs) for info on the basis and\n"
    "coefficient objects.\n"
    "Data packing\n"
    "------------\n"
    "All slices and volumes are returned as numpy.ndarrays in row-major\n"
    "order.  That is, the first index is x, the second is y, and the\n"
    "third is z.  The ranks are specified by the 'gridsize' array with\n"
    "(nx, ny, nz) as input to the FieldGenerator constructor.\n\n";

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
  
  f.def("lines", &Field::FieldGenerator::lines,
	"Return a dictionary of arrays (1d numpy arrays) indexed by "
	"time and field type", py::arg("basis"), py::arg("coefs"),
	py::arg("beg"), py::arg("end"), py::arg("num"));
  
  f.def("histo2d", &Field::FieldGenerator::histogram2d,
	"Return a density histogram (2d numpy arrays)",
	py::arg("reader"),
	py::arg("center") = std::vector<double>(3, 0.0));

  f.def("histo1d", &Field::FieldGenerator::histogram1d,
	"Return a radial density histogram (array) for a chosen projection "
	"indicated by a string: \"xy\", \"xz\", \"yz\", \"r\"",
	py::arg("reader"), py::arg("rmax"), py::arg("nbins"),
	py::arg("projection"),
	py::arg("center") = std::vector<double>(3, 0.0));

  f.def("file_lines", &Field::FieldGenerator::file_lines,
	"Write field arrays to files using the supplied string prefix. "
	"The files are ascii tables with column headers.", py::arg("basis"),
	py::arg("coefs"), py::arg("beg"), py::arg("end"),
	py::arg("num")=1000, py::arg("filename"), py::arg("dir")=".");

  f.def("file_slices", &Field::FieldGenerator::file_slices,
	"Write 2d field grids to files using the supplied string prefix. "
	"The files will be ascii or VTK rectangular grid files (if you "
	"compiled with VTK)",
	py::arg("basis"), py::arg("coefs"), py::arg("filename"),
	py::arg("dir")=".");


  f.def("volumes", [](FieldGenerator& A,
		      BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs)
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
	"Write 3d field grids to files using the supplied string prefix. "
	"The files will be ascii or VTK rectangular grid files (if you "
	"compiled with VTK)",
	py::arg("basis"), py::arg("coefs"), py::arg("filename"),
	py::arg("dir")=".");
}
