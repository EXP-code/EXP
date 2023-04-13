#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Koopman.H>

namespace py = pybind11;
#include <TensorToArray.H>

void EDMDtoolkitClasses(py::module &m) {

  m.doc() =
    "Extended Dynamical Mode Decomposition (EDMD) class bindings\n\n"
    "Koopman\n"
    "=======\n"
    "This Koopman class implements EDMD, using the exact DMD algorith\n"
    "from Tu et al. (2014) using the EXP coefficients and your optional\n"
    "auxiliary data.  Just as in expMSSA, instances return coefficient\n"
    "sets for one or more individual EMD modes.  These may, in turn, be\n"
    "used with the FieldGenerator for visualization\n\n"
    "Coefficient update and memory usage\n"
    "-----------------------------------\n"
    "A key feature of Koopman operator theory is the identification of\n"
    "the underlying dynamical structure in the measurement of state space.\n"
    "In this case, the measurement of phase space is the collection of\n"
    "scalars that represent the density-potential fields.  There is the\n"
    "possibility and maybe even the likelihood that the mismatch of the\n"
    "fields to the true nature of the dynamics will generate spurious\n"
    "models.  This is true with mSSA as well.  In fact, mSSA is also an"
    "approximation to the Koopman operator, so performance will be similar.\n"
    "As as in in mSSA, this class allows you to gain insight from a set of\n"
    "coefficients reconstructed from dominant EDMD modes. After\n"
    "identifying interesting modes and associated eigenvalues, this\n"
    "is done in 'Koopman' in two steps:\n"
    "  1. A call to 'reconstruct(list)' with a list of eigenvalue\n"
    "     indices reconstruct the data series and saves those series\n"
    "     to working (possibly detrended) vectors\n"
    "  2 'getReconstructed()' returns a dictionary of name strings that\n"
    "     point to coefficient objects (Coefs). These may be used for\n"
    "     diagnostics (e.g. visualizing fields).\n"
    "This structure exactly parallels the implementation in 'expMSSA'\n"
    "The split between the the 'reconstruct()' and the 'getReconstructed()\n"
    "members allows one to examine the detrended series with the selected\n"
    "eigenvalues before updating the coefficient data.  The 'reconstruct()\n"
    "call creates the internal working channel data, leaving the initial\n"
    "coefficient set untouched. The 'getReconstructed()' call replaces the\n"
    "coefficient data in the Coefs object passed to Koopman, as well as\n"
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
    "Koopman has a number of internal configuration parameters, listed\n"
    "below. Most people will not need them, but they do exist should\n"
    "you really want to dig in.  These are a subset of the same parameter\n"
    "available in expMSSA.\n\n"
    "To change the default values, pass a YAML database of parameters\n"
    "to the Koopman constructor.  This is a convenient way to have control\n"
    "over some fine-tuning parameters or experimental features without a\n"
    "a change to the API.  They are all optional but a few of these will\n"
    "be useful.  They are all YAML key-value pairs.  The first group\n"
    "below are simple boolean toggles; their presence turns on a feature\n"
    "independent of the value.  The Python script samples show how to\n"
    "construct these simple YAML configurations on the fly. A simple example\n"
    "is also given below.  The boolean parameters are listed below by my\n"
    "guess of their usefulness to most people:\n\n"
    "  verbose: false        Whether there is report or not\n"
    "  Jacobi: true          Use the Jacobi SVD rather than the Random\n"
    "                        approximation algorithm from Halko, Martinsson,\n"
    "                        and Tropp (RedSVD). This is quite accurate but\n"
    "                        _very_ slow\n"
    "  BDCSVD: true          Use the Binary Divide and Conquer SVD rather\n"
    "                        rather than RedSVD; this is faster and more\n"
    "                        accurate than the default RedSVD but slower\n"
    "  RedSym: true          Use the randomized symmetric eigenvalue solver\n"
    "                        RedSym rather rather than RedSVD. The main use\n"
    "                        for these various SVD algorithm toggles is for\n"
    "                        checking the accuracy of the default randomized\n"
    "                        matrix methods.\n"
    "The following parameters take values, defaults are given in ()\n\n"
    "  output: str(exp_mssa) Prefix name for output files\n\n"
    "The 'output' value is only used if 'writeFiles' is specified, too.\n"
    "A simple YAML configuration for Koopman might look like this:\n"
    "---\n"
    "BDCSVD: true\n"
    "project: true\n"
    "...\n\n"
    "If you find that you are using these YAML parameters often, let me\n"
    "know and I can promote them to the main API\n\n"
    "Save/Restore\n"
    "------------\n"
    "Koopman analyses can be computationally intensive so we include a way\n"
    "to serialize the state to an HDF5 file.  The saveState(prefix) member\n"
    "takes a prefix file parameter and creates the file <prefix>_edmd.h5.\n"
    "You may restore the state by recreating your Koopman instance with\n"
    "the identical input data and keys and then using restoreState(prefix)\n"
    "to read the Koopman analysis from the HDF5 file.  The restore step will\n"
    "check that your data has the same dimension, same parameters, and key\n"
    "list so it should be pretty hard to fool, but not impossible.  The\n"
    "computation is much less expensive, generally, than MSSA so saving the\n"
    "analysis is not strictly necessary, even on a laptop.\n\n"
    "Computational notes\n"
    "-------------------\n"
    "Some of the linear algebra operations can parallelize itself using\n"
    "OpenMP threads.  If your compilter has OpenMP, please enable it (i.e.\n"
    "-fopenmp when compiling in GNU or and set the appropriate environment\n"
    "(e.g. 'export OMP_NUM_THREADS=X' for Bash). The initial Koopman analysis\n"
    "can also be memory intensive, depending on the size of the time series\n"
    "and the number of channels.  Although it requires much less memory than\n"
    "the equivalent mSSA analysis. If memory is a problem, try running the\n"
    "analyses on a larger memory compute node and save using the 'saveState()'\n"
    "member to an HDF5 file and reread those files on a local machine using\n"
    "the 'restoreState()' member.\n\n";

  using namespace MSSA;

  py::class_<MSSA::Koopman, std::shared_ptr<MSSA::Koopman>> f(m, "Koopman");

  f.def(py::init<const mssaConfig&, int, const std::string>(),
	"Constructor\nconfig\tis the input database of components"
	"\numev\tis the default number of eigenvalues to compute"
	"\nflag\tis a YAML stanza of parameter values",
	py::arg("config"),
	py::arg("numev"),
	py::arg("flags") = "");

  f.def("eigenvalues", &Koopman::eigenvalues,
	"Return the vector of eigenvalues from the EDMD analysis. Note that\n"
	"these eigenvalues are complex.  It may be helpful to separate the\n"
	"magnitude and phase using Numpy's 'absolute' and 'angle' functions.\n");

  f.def("channelDFT", &Koopman::channelDFT,
	"Returns the frequency and the DFT of the selected data channels.\n"
	"Unlike mSSA, there is no meaningful counterpart to DFT analysis of\n"
	"PCs.");

  f.def("reconstruct", &Koopman::reconstruct,
	"Reconstruct the data channels with the provided list of eigenvalue "
	"indices", py::arg("evlist"));

  f.def("getReconstructed", &Koopman::getReconstructed,
	"Return the reconstucted time series in the orginal coefficient form\n"
	"that may be used in basis classes.  Note: the reconstructed data\n"
	"will overwrite the memory of the original coefficient data.");

  f.def("contrib", &Koopman::contributions,
	"Computes the relative contribution of each mode to the coefficient\n"
	"series and the breakdown of the coefficient series to each mode.\n"
	"The views normed on columns and rows are returned as a tuple\n"
	"of 2d arrays. These are intended to be plotted using 'imshow'.");

  f.def("saveState", &Koopman::saveState,
	"Save current EDMD state to an HDF5 file with the given prefix",
	py::arg("prefix"));

  f.def("restoreState", &Koopman::restoreState,
	"Restore current EDMD state from an HDF5 file with the given prefix.\n"
	"To use this, the Koopman instance must be constructed with the same\n"
	"data and parameters as the save stated.  The restoreState routine will\n"
	"check for the same data dimension and trend state but can not sure\n"
	"complete consistency.",
	py::arg("prefix"));

  f.def("getModes", &Koopman::getModes,
	"Access to detrended reconstructed channel series.  Also see the member\n"
	"'contributions' to visualize the support from each EDMD mode to the\n"
	"coefficient series.");

  f.def("getAllKeys", &Koopman::getAllKeys,
	"As in mSSA, provides a list of all internal channel keys (for reference)");

}
