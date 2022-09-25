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
    "Coefficient update and memory usage\n"
    "-----------------------------------\n"
    "A key feature of MSSA inference is the field generation from\n"
    "the coefficients reconstructed from dominant eigenvalues. After\n"
    "identifying interesting PCs and associated eigenvalues, this\n"
    "is done in 'expMSSA' in two steps:\n"
    "  1. A call to 'reconstruct(list)' with a list of eigenvalue\n"
    "     indices reconstruct the data series and saves those series\n"
    "     to working (possibly detrended) vectors\n"
    "  2 'getReconstructed()' returns a dictionary of name strings that\n"
    "     point to coefficient objects (Coefs). These may be used for\n"
    "     diagnostics (e.g. visualizing fields).\n"
    "The split between the the 'reconstruct()' and the 'getReconstructed()\n"
    "members allows one to examine the detrended series with the selected\n"
    "eigenvalues before updating the coefficient data.  The 'reconstruct()\n"
    "call creates the internal working channel data, leaving the initial\n"
    "coefficient set untouched. The 'getReconstructed()' call replaces the\n"
    "coefficient data in the Coefs object passed to expMSSA, as well as\n"
    "returning a dictionary of name strings and coefficient containers.\n"
    "Subsequent calls to 'getReconstructed' will overwrite the previous\n"
    "updates.  In practice, the coefficients object that you made in Python\n"
    "will always contain the latest updates from the last 'getReconstructed'\n"
    "call.  This allows you to update channels incrementally and the updated\n"
    "copy is available to Python as the originally passed reference.  This\n"
    "strategy prevents your stack from growing very large by efficiently\n"
    "reusing previously allocated storage. If you want to save a copy of\n"
    "the original coefficients, use the 'deepcopy()' member for Coefs\n"
    "before the reconstruction.  The original 'factory' call and each\n"
    "call to 'deepcopy()' will make a new set that will be freed when\n"
    "the reference to the coefficient instance disappears.  This allows\n"
    "some explicit control over your memory use.\n\n"
    "Configuration parameters\n"
    "------------------------\n"
    "expMSSA has a number of internal configuration parameters, listed\n"
    "below. Most people will not need them, but they do exist should\n"
    "you really want to dig in.\n\n"
    "To change the default expMSSA, pass a YAML database of parameters\n"
    "to the expMSSA constructor.  This is a convenient way to have control\n"
    "over some fine-tuning parameters or experimental features without a\n"
    "a change to the API.  They are all optional but a few of these will\n"
    "be useful.  They are all YAML key-value pairs.  The first group\n"
    "below are simple boolean toggles; their presence turns on a feature\n"
    "independent of the value.  The Python script samples show how to\n"
    "construct these simple YAML configurations on the fly. A simple example\n"
    "is also given below.  The boolean parameters are listed below by my\n"
    "guess of their usefulness to most people:\n\n"
    "  writeCov: true        Write the covariance matrix to a file for\n"
    "                        diagnostics since this is not directly\n"
    "                        available from the interface\n"
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
    "                        RedSym rather rather than RedSVD. The main use\n"
    "                        for these various SVD algorithm toggles is for\n"
    "                        checking the accuracy of the default randomized\n"
    "                        matrix methods.\n"
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
    "BDCSVD: true\n"
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
    "indices in eigenvalue-PC index order. An empty argument implies no\n"
    "that the will be no PCs in the reconstruction, only the mean values.\n"
    "The two expMSSA steps will often need to be used iteratively.  First,\n"
    "to learn about the primary signals in your data and catagorize the PCs\n"
    "into groups.  Second, a new analysis that provides separation of your\n"
    "main distinct signals into individual reconstructions.\n\n"
    "Save/Restore\n"
    "------------\n"
    "MSSA analyses can be computationally intensive so we include a way\n"
    "to serialize the state to an HDF5 file.  The saveState(prefix) member\n"
    "takes a prefix file parameter and creates the file <prefix>_mssa.h5.\n"
    "You may restore the state by recreating your expMSSA instance with\n"
    "the identical input data and keys and then using restoreState(prefix)\n"
    "to read the MSSA analysis from the HDF5 file.  The restore step will\n"
    "check that your data has the same dimension, same parameters, and key\n"
    "list so it should be pretty hard to fool, but not impossible.\n\n"
    "Computational notes\n"
    "-------------------\n"
    "Some of the linear algebra operations and the MSSA reconstruction\n"
    "can parallelize itself using OpenMP threads.  If your compilter has\n"
    "OpenMP, please enable it (i.e. -fopenmp when compiling in GNU or\n"
    "and set the appropriate environment (e.g. 'export OMP_NUM_THREADS=X'\n"
    "for Bash).  This will be especially useful for the 'reconstruct()'\n"
    "step.  The initial MSSA analysis can also be memory intensive,\n"
    "depending on the size of the time series and the number of channels.\n"
    "One possible strategy is to run the analyses on a larger memory\n"
    "compute node and save the results using the 'saveState()' member to\n"
    "an HDF5 file and reread those files on a local machine using the\n"
    "'restoreState()' member.\n\n";
  
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
	"indices (a group).", py::arg("evlist"));

  f.def("getReconstructed", &expMSSA::getReconstructed,
	"Return the reconstucted time series in the orginal coefficient form\n"
	"that may be used in basis classes.  Note: the reconstructed data\n"
	"will overwrite the memory of the original coefficient data.");

  f.def("wCorr", &expMSSA::wCorr,
	"Get the w-correlation matrix for the selected component and channel\n"
	"key.  Returns the combined cosine+sine correlation for complex types\n"
	"for viewing (e.g.) with 'imshow'",
	py::arg("name"), py::arg("key"));

  f.def("wCorrKey", &expMSSA::wCorrKey,
	"Get the w-correlation matrix for the selected component and channel\n"
	"key extended by the cosine/sine index if the channel is complex and\n"
	"the component index. Try plotting using 'imshow'.\n", py::arg("key"));

  f.def("wCorrAll", &expMSSA::wCorrAll,
	"Get the w-correlation matrix for all channels in the reconstruction.\n"
	"These can be nicely plotted using 'imshow'.");

  f.def("wcorrPNG", &expMSSA::wcorrPNG,
	"Create wcorrlation matricies and output PNG image representations");

  f.def("kmeans", &expMSSA::kmeans,
	"Perform a k-means analysis on the reconstructed trajectory matrices\n"
	"to provide grouping insight.  This will write to the standard output\n"
	"by default.  Set toFile=True to write ot a file.  The file name will\n"
	"be derived from the 'output' parameter.",
	py::arg("clusters")=4,
	py::arg("toTerm")=true,
	py::arg("toFile")=false);

  f.def("contrib", &expMSSA::contributions,
	"Computes the relative contribution of each PC to the coefficient\n"
	"series and the breakdown of the coefficient series to each PC.\n"
	"The views normed on columns and rows are returned as a tuple\n"
	"of 2d arrays. These are intended to be plotted using 'imshow'.");

  f.def("saveState", &expMSSA::saveState,
	"Save current MSSA state to an HDF5 file with the given prefix",
	py::arg("prefix"));

  f.def("restoreState", &expMSSA::restoreState,
	"Restore current MSSA state from an HDF5 file with the given prefix.\n"
	"To use this, the expMSSA instance must be constructed with the same\n"
	"data and parameters as the save stated.  The restoreState routine will\n"
	"check for the same data dimension and trend state but can not sure\n"
	"complete consistency.",
	py::arg("prefix"));

  f.def("getTotVar", &expMSSA::getTotVar,
	"Variance value used for normalizing coefficient series");

  f.def("getTotPow", &expMSSA::getTotPow,
	"Power value used for normalizing coefficient series");

  f.def("getRC", &expMSSA::getRC,
	"Access to detrended reconstructed channel series by "
	"internal key");

  f.def("getRCkeys", &expMSSA::getRCkeys,
	"Provides a list of internal keys for accessing the "
	"detrended channel series using getRC()");

}
