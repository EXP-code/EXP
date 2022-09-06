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

  m.doc() = "Multivariate Singular Spectrum Analysis (MSSA) class bindings.\n"
    "The expMSSA class analyzes and separates the temporal patterns\n"
    "in your coefficients and auxiliary data.  It can return coeffi\n"
    "cients sets for one or more individual patterns that may be\n"
    "used with the FieldGenerator for visualization\n\n"
    "expMSSA has a number of internal configuration parameters, listed\n"
    "below. Most people will not need them, but they do exist should\n"
    "you really want to dig in.\n\n"
    "Configuration parameters\n"
    "------------------------\n"
    "expMSSA can be configured with a YAML database of parameters.\n"
    "They are all optional but a few of these will be useful. They\n"
    "are all YAML key-value pairs.  The first group below are simple\n"
    "boolean toggles.  That is there presence turns on the feature\n"
    "independent of the value.  The Python scripts samples show how\n"
    "to construct these simple YAML configurations on the fly. A\n"
    "Simple example is also given below.  The boolean parameters\n"
    "are listed below in order of my guess of their usefulness to\n"
    "most people:\n\n"
    "  chatty: true          Prints eigenvalues, PC, etc. to stdout.\n"
    "                        This one is useful if you want to see a nice\n"
    "                        formatted output of eigenvalues, PCs, k-means\n"
    "                        output, etc, rather than formating your own.\n"
    "  writeFiles: true      Writes most of the 'chatty' output to files\n"
    "  writeCov: true        Write the covariance matrix to a file\n"
    "  Jacobi: true          Use the Jacobi SVD rather than the Random\n"
    "                        approximation algorithm from Halko, Martinsson,\n"
    "                        and Tropp (RedSVD). This is quite accurate but\n" 
    "                        _very_ slow\n"
    "  BDCSVD: true          Use the Binary Divide and Conquer SVD rather\n"
    "                        rather than RedSVD; this is faster and more\n"
    "                        accurate than the default RedSVD but slower\n"
    "  Traj: true            Perform the SVD of the trajectory matrix\n"
    "                        rather than the more computationally SVD\n"
    "                        of the trajectory matrix.  Do not use this\n"
    "                        with the default RedSVD algorithm.  Use either\n"
    "                        either Jacobi or BDCSVD for accuracy.\n"
    "  allchan: true         Perform k-mean clustering analysis using all\n"
    "                        channels simultaneously\n"
    "  distance: true        Compute w-correlation matrix PNG images using\n"
    "                        w-distance rather than correlation\n"
    "  flip: true            Exchanges the x-y axes in PNG plots\n\n"
    "The following parameters take values, defaults are given in ()\n\n"
    "  evtol: double(0.01)   Truncate by the given cumulative p-value\n"
    "  output: str(exp_mssa) Prefix name for output files\n\n"
    "The 'output' value is only used if 'writeFiles' is specified, too.\n"
    "A simple YAML configuration for expMSSA might look like this:\n"
    "---\n"
    "chatty: true\n"
    "output: my_test_run\n"
    "evtol: 0.01\n"
    "...\n\n"
    "Grouping is an important part of MSSA analysis.  MSSA often spreads\n"
    "the same signal between more than one PC.  The w-correlation matrix,\n"
    "the k-means analyses and similarities of the PC frequency spectra are\n"
    "are ways of diagnosing groups of related components. You may group\n"
    "using reconstruct([list]), where [list] is a Python list of component\n"
    "indices in eigenvalue/PC index order. The default, an empty argument,\n"
    "is no grouping; that is, reconstruct for all eigenvalues.  MSSA will\n"
    "often need to be used iteratively.  First, to learn about the primary\n"
    "signals in your data and catagorize the PCs into groups.  Second, a\n"
    "new analysis that provides separation of your main distinct signals\n"
    "into individual reconstructions\n\n";
  
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
  
  f.def("getU", &expMSSA::getU,
	"Return the right-singular) vectors from the MSSA analysis "
	"which describe the contribution of each channel to each PC");

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
	"indices (a group).", py::arg("evlist")=std::vector<int>());

  f.def("getReconstructed", &expMSSA::getReconstructed,
	"Return the reconstucted time series in the orginal coefficient form "
	"that may be used in basis classes", py::arg("zero") = false);

  f.def("wCorr", &expMSSA::wCorr,
	"Get the w-correlation matrix for the selected component and channel "
	"key");

  f.def("wCorrAll", &expMSSA::wCorrAll,
	"Get the w-correlation matrix for all channels in the reconstruction");

  f.def("wcorrPNG", &expMSSA::wcorrPNG,
	"Create wcorrlation matricies and output PNG image representations");

  f.def("kmeans", &expMSSA::kmeans,
	"Perform a k-means analysis on the reconstructed trajectory matrices "
	"to provide grouping insight.  The file name will be derived from the "
	"'output' parameter", py::arg("clusters")=4);

  f.def("contrib", &expMSSA::contributions,
	"Computes the relative contribution of each PC to the coefficient "
	"series and the breakdown of the coefficient series to each PC. "
	"The results are rendered as PNG images");
}
