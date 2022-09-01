#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <expMSSA.H>

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

void MSSAtoolkitClasses(py::module &m) {

  m.doc() = "Multivariate Singular Spectrum Analysis (MSSA) class bindings\n"
    "The expMSSA class analyzes and separates the temporal patterns\n"
    "in your coefficients and auxilliary data.  It can return coeffi\n"
    "cients sets for one or more individual patterns that may be\n"
    "used with the FieldGenerator for visualiation.";

  using namespace MSSA;

  py::class_<MSSA::expMSSA, std::shared_ptr<MSSA::expMSSA>> f(m, "expMSSA");

  f.def(py::init<const mssaConfig&, int, int, const std::string>(),
	"Constructor\nconfig\tis the input database of components"
	"\nwindow\tis the length of the cross-correlation interval"
	"\nnumpc\tis the default number of eigenvalues to compute"
	"\nflag\tis a YAML stanza of parameter values",
	py::arg("config"),
	py::arg("window"),
	py::arg("numpc"),
	py::arg("flags") = "");

  f.def("eigenvalues", &expMSSA::eigenvalues,
	"Return the vector of eigenvalues from the MSSA analysis");
  
  f.def("getPC", &expMSSA::getPC,
	"Return the principle component (left-singular) vectors from the MSSA "
	"analysis which describe the key temporal variation");

  f.def("pcDFT", &expMSSA::pcDFT,
	"Return the DFT of the principal component vectors for quantifying "
	"temporal power distribution", py::arg("freq"), py::arg("period"));

  f.def("channelDFT", &expMSSA::channelDFT,
	"Return the DFT of the selected data channels for comparison with "
	"the PC power",	py::arg("freq"), py::arg("period"));

  f.def("reconstruct", &expMSSA::reconstruct,
	"Reconstruct the data channels with the provided list of eigenvalue "
	"indices", py::arg("evlist"));

  f.def("getReconstructed", &expMSSA::getReconstructed,
	"Return the reconstucted time series in the orginal coefficient form "
	"that may be used in basis classes", py::arg("zero") = false);

  f.def("wCorr", &expMSSA::wCorr,
	"Get the w-correlation matrix for the selected component and channel key");

  f.def("wCorrAll", &expMSSA::wCorrAll,
	"Get the w-correlation matrix for all channels in the reconstruction");

}
