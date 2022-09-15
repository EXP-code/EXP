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

  m.doc() =
    "Multivariate Singular Spectrum Analysis (MSSA) class bindings\n\n"
    "expMSSA\n"
    "=======\n"
    "The expMSSA class analyzes and separates the temporal patterns\n"
    "in your coefficients and auxiliary data.  Instances can return\n"
    "coefficient sets for one or more individual patterns that may be\n"
    "used with the FieldGenerator for visualization\n\n"
    "A note on memory usage\n"
    "----------------------\n"
    "A key feature of MSSA inference is the field generation from\n"
    "the coefficients reconstructed from dominant eigenvalues. After\n"
    "identifying interesting PCs and associated eigenvalues, this\n"
    "is done in 'expMSSA' in two steps:\n"
    "  1. A call to 'reconstruction(list)' with a list of eigenvalue\n"
    "     indices reconstruct the data series and saves those series\n"
    "     to working vectors\n"
    "  2 'getReconstructed()' returns a dictionary of name strings that\n"
    "     point to coefficient objects (Coefs). These may be used for\n"
    "     diagnostics (e.g. visualizing fields).  These new sets are\n"
    "     are shallow copies of the original Coefs data updated by from\n"
    "     the working vectors produced by the last call to the function\n"
    "     'reconstruction(list)'.\n"
    "This strategy prevents your stack from growing large by efficiently\n"
    "reusing previously allocated storage. If you want to save a copy of\n"
    "the original coefficients, use the 'deepcopy()' member for Coefs\n"
    "before the reconstruction.  The original 'factory' call and each\n"
    "call to 'deepcopy()' will make a new set that will be freed when\n"
    "the reference to the coefficient instance disappears.  This allows\n"
    "some explicit control over your memory use.\n\n"
    "expMSSA has a number of internal configuration parameters, listed\n"
    "below. Most people will not need them, but they do exist should\n"
    "you really want to dig in.\n\n"
    "Configuration parameters\n"
    "------------------------\n"
    "expMSSA can be configured with a YAML database of parameters.\n"
    "I wrote this is a convenient way to have control over some\n"
    "tuning parameters or experimental features without having to\n"
    "change the API.  They are all optional but a few of these will\n"
    "be useful.  The most useful might be 'writeFiles' which turns on\n"
    "writing chatty MSSA results to files and 'output' to specify the\n"
    "output-file prefix. They are all YAML key-value pairs.  The first\n"
    "group below are simple boolean toggles.  That is there presence\n"
    "turns on the feature independent of the value.  The Python scripts\n"
    "samples show how to construct these simple YAML configurations on\n"
    "the fly. A simple example is also given below.  The boolean para-\n"
    "meters are listed below in order of my guess of their usefulness to\n"
    "most people:\n\n"
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
    "  RedSym: true          Use the randomized symmetric eigenvalue solver\n"
    "                        RedSym rather rather than RedSVD\n"
    "  allchan: true         Perform k-mean clustering analysis using all\n"
    "                        channels simultaneously\n"
    "  distance: true        Compute w-correlation matrix PNG images using\n"
    "                        w-distance rather than correlation\n"
    "  flip: true            Exchanges the x-y axes in PNG plots\n\n"
    "The following parameters take values, defaults are given in ()\n\n"
    "  evtol: double(0.01)   Truncate by the given cumulative p-value in\n"
    "                        chatty mode\n"
    "  output: str(exp_mssa) Prefix name for output files\n\n"
    "The 'output' value is only used if 'writeFiles' is specified, too.\n"
    "A simple YAML configuration for expMSSA might look like this:\n"
    "---\n"
    "chatty: true\n"
    "output: my_test_run\n"
    "evtol: 0.01\n"
    "...\n\n"
    "If you find that you are using these YAML parameters often, let me\n"
    "know and I can promote them to the main API\n\n"
    "Grouping\n"
    "--------\n"
    "MSSA often spreads the same signal between more than one PC. We try to\n"
    "group eigenvalue-PC pairs that seem to have similar temporatl behavior\n"
    "in the final stage before reconstruction.  The w-correlation matrix,\n"
    "the k-means analyses and similarities of the PC frequency spectra are\n"
    "are ways of diagnosing groups of related components. You may group\n"
    "using 'reconstruct(list)', where 'list' is a Python list of component\n"
    "indices in eigenvalue-PC index order. The default, an empty argument,\n"
    "is no grouping; that is, reconstruct for all eigenvalues.  MSSA will\n"
    "often need to be used iteratively.  First, to learn about the primary\n"
    "signals in your data and catagorize the PCs into groups.  Second, a\n"
    "new analysis that provides separation of your main distinct signals\n"
    "into individual reconstructions\n\n"
    "Save/Restore\n"
    "------------\n"
    "MSSA analyses can be computationally intensive so we include a way\n"
    "to serialize the state to an HDF5 file.  The saveState(prefix) member\n"
    "takes a prefix file parameter and creates the file <prefix>_mssa.h5.\n"
    "You may restore the state by recreating your expMSSA instance with\n"
    "the identical input data and keys and then using restoreState(prefix)\n"
    "to read the MSSA analysis from the HDF5 file.  The restore step will\n"
    "check that your data has the same dimension, same parameters, and key\n"
    "list so it should be pretty hard to fool, but not impossible.\n\n";
  
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
  
  f.def("cumulative", &expMSSA::cumulative,
	"Return a cumulatively summed vector of eigenvalues from the\n"
	"MSSA analysis");
  
  f.def("getU", &expMSSA::getU,
	"Return the right-singular) vectors from the MSSA analysis\n"
	"which describe the contribution of each channel to each PC");

  f.def("getPC", &expMSSA::getPC,
	"Return the principle component (left-singular) vectors from the MSSA\n"
	"analysis which describe the key temporal variation");

  f.def("pcDFT", &expMSSA::pcDFT,
	"Return the DFT of the principal component vectors for quantifying\n"
	"temporal power distribution", py::arg("freq"), py::arg("period"));

  f.def("channelDFT", &expMSSA::channelDFT,
	"Return the DFT of the selected data channels for comparison with\n"
	"the PC power",	py::arg("freq"), py::arg("period"));

  f.def("reconstruct", &expMSSA::reconstruct,
	"Reconstruct the data channels with the provided list of eigenvalue "
	"indices (a group).", py::arg("evlist")=std::vector<int>());

  f.def("getReconstructed", &expMSSA::getReconstructed,
	"Return the reconstucted time series in the orginal coefficient form\n"
	"that may be used in basis classes.  Note: the reconstructed data\n"
	"will overwrite the memory of the original coefficient data.");

  f.def("wCorr", &expMSSA::wCorr,
	"Get the w-correlation matrix for the selected component and channel\n"
	"key.  Returns the combined cosine+sine correlation for complex types",
	py::arg("name"), py::arg("key"));

  f.def("wCorrKey", &expMSSA::wCorrKey,
	"Get the w-correlation matrix for the selected component and channel\n"
	"key extended by the cosine/sine index if the channel is complex and\n"
	"the component index", py::arg("key"));

  f.def("wCorrAll", &expMSSA::wCorrAll,
	"Get the w-correlation matrix for all channels in the reconstruction");

  f.def("wcorrPNG", &expMSSA::wcorrPNG,
	"Create wcorrlation matricies and output PNG image representations");

  f.def("kmeans", &expMSSA::kmeans,
	"Perform a k-means analysis on the reconstructed trajectory matrices\n"
	"to provide grouping insight.  This will write to the standard output\n"
	"by default.  Set toFile=True to write ot a file.  The file name will\n"
	"be derived from the 'output' parameter",
	py::arg("clusters")=4,
	py::arg("toTerm")=true,
	py::arg("toFile")=false);

  f.def("contrib", &expMSSA::contributions,
	"Computes the relative contribution of each PC to the coefficient "
	"series and the breakdown of the coefficient series to each PC. "
	"The results are rendered as PNG images");

  f.def("saveState", &expMSSA::saveState,
	"Save current MSSA state to an HDF5 file with the given prefix",
	py::arg("prefix"));

  f.def("restoreState", &expMSSA::restoreState,
	"Restore current MSSA state from an HDF5 file with the given prefix.\n"
	"To use this, the expMSSA instance must be constructed with the same\n"
	"data and parameters as the save stated.  The restoreState routine will\n"
	"check for the same data dimension and trend state but can not sure\n"
	"complete consistency.\n",
	py::arg("prefix"));

  f.def("getTotVar", &expMSSA::getTotVar,
	"Variance value used for normalizing coefficient series");

  f.def("getTotPow", &expMSSA::getTotPow,
	"Power value used for normalizing coefficient series");

  f.def("zeroReconstructed", &expMSSA::zeroReconstructed,
	"Zero reconstruction for all keys. Designed for removing\n"
	"the influence of coefficients deemed to be null and\n"
	"rewriting for playback");

}
