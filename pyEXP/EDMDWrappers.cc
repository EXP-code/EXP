#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Koopman.H>
#include <KoopmanRKHS.H>
#include <LiouvilleRKHS.H>

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
    "models.  This is true with mSSA as well.  In fact, mSSA is also an\n"
    "approximation to the Koopman operator, so the philosophy is similar.\n"
    "My experience to date suggests that mSSA has much better performance\n"
    "and separation of signals.  If you try this and have good success,\n"
    "please do contact me.  As as in in mSSA, this class allows you to gain\n"
    "insight from a set of coefficients reconstructed from dominant EDMD\n"
    "modes. After identifying interesting modes and associated eigenvalues,\n"
    "this is done in 'Koopman' in two steps:\n"
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
    "  project: true         Use the classic DMD projected modes rather than\n"
    "                        the exact DMD modes, as described by Tu et al.\n"
    "                        2014.\n"
    "  power: true           Write partial power contributions into a file in\n"
    "                        a ascii table format if set to 'true'.  Default\n"
    "                        is 'false'\n"
    "The following parameters take values, defaults are given in ()\n\n"
    "  output: sting         Prefix name for output files.  The default is\n"
    "                        'exp_edmd'.\n\n"
    "The 'output' value is used by 'getContributions()' and 'channelDFT()'\n"
    "if the 'power' options is set.\n"
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

  py::class_<MSSA::Koopman, std::shared_ptr<MSSA::Koopman>>
    f(m, "Koopman");

  f.def(py::init<const mssaConfig&, int, const std::string>(),
	R"(
        Koopman operator approximatation class

        Parameters
        ----------
        config : mssaConfig
	     the input database of components
	numev : int
             the default number of eigenvalues to compute"
	flag : str
            YAML stanza of parameter values

        Returns
        -------
        Koopman instance

        Notes
        -----
        The configuration should be in the format:

        {'example': (coefs, keylst, [])}

        where keylst is a list of selected PCs.  With a SphericalBasis
        for example, the keylst should have the format:

        [[l1, m1, n1], [l2, m2, n2], ...]

        Each sublist represents a PC, where l, m, and n are the 
        spherical harmonics basis parameters.
        )",
	py::arg("config"),
	py::arg("numev"),
	py::arg("flags") = "");

  f.def("eigenvalues", &Koopman::eigenvalues,
	R"(
        Vector of eigenvalues from the EDMD analysis. 

        Returns
        -------
        numpy.ndarray
            the vector of eigenvalues

        Notes
        -----
        Note that these eigenvalues are complex.  It may be helpful 
        to separate the magnitude and phase using Numpy's 'absolute' 
        and 'angle' functions.
       )");

  f.def("channelDFT", &Koopman::channelDFT,
	R"(
        Returns the frequency and the DFT of the selected data channels.

        Returns
        -------
        tuple(numpy.ndarray, numpy.ndarray) : frequencies and spectral power
        
        Notes
        -----
	Unlike mSSA, there is no meaningful counterpart to DFT analysis of PCs.
        )");

  f.def("reconstruct", &Koopman::reconstruct,
	R"(
        Reconstruct the data channels with the provided list of eigenvalue indices

        Parameters
        ----------
        evlist : list(int)
            eigenvalue indices to include in the reconstruction
        
        Returns
        -------
        None
        )", py::arg("evlist"));

  f.def("getReconstructed", &Koopman::getReconstructed,
	R"(
        Reconstruct the channels

        Reconstucted time series in the orginal coefficient form
	that may be used in basis classes.  

        Notes
        -----
	the reconstructed data will overwrite the memory of the original coefficient data
        )");

  f.def("contrib", &Koopman::contributions,
	R"(
        Contributions per channel

        Computes the relative contribution of each mode to the coefficient
	series and the breakdown of the coefficient series to each mode.

        Notes
        -----
	The views normed on columns and rows are returned as a tuple
	"of 2d arrays. These are intended to be plotted using 'imshow'.
        )");

  f.def("saveState", &Koopman::saveState,
	R"(
        Save current EDMD state to an HDF5 file with the given prefix

        Parameters
        ----------
        prefix : str
            output filename prefix

        Returns
        -------
        None
        )", py::arg("prefix"));

  f.def("restoreState", &Koopman::restoreState,
	R"(
        Restore current EDMD state from an HDF5 file

        Parameters
        ----------
        prefix : str
            input filename prefix

        Notes
        -----
	The Koopman instance must be constructed with the same data and parameters
        as the saved state.  The restoreState routine will check for the same 
        data dimension and trend state but can not sure	complete consistency.
        )", py::arg("prefix"));

  f.def("getModes", &Koopman::getModes,
	R"(
        Access to detrended reconstructed channel series.

        Returns
        -------
        numpy.ndarray
            the EDMD modes

        See also
        --------
        Use in conjunction with 'contributions' to visualize the support 
        from each EDMD mode to the coefficient series.
        )");

  f.def("getAllKeys", &Koopman::getAllKeys,
	R"(
        Provides a list of all internal channel keys (for reference)

        Returns
        -------
        list(Key)
            list of keys in the format described in config file

        See also
        --------
        Koopman
        )");

  py::class_<MSSA::KoopmanRKHS, std::shared_ptr<MSSA::KoopmanRKHS>>
    g(m, "KoopmanRKHS");

  g.def(py::init<const mssaConfig&, double, int, const std::string>(),
	R"(
        Koopman RKHS operator approximatation class

        Parameters
        ----------
        config : mssaConfig
	    the input database of components
        tol : double
            tolerance for eigenvalue decomposition
        count : int
            top count of eigenvalues
	flag : str
            YAML stanza of parameter values

        Returns
        -------
        Koopman RKHS instance

        Notes
        -----
        The configuration should be in the format:

        {'example': (coefs, keylst, [])}

        where keylst is a list of selected PCs.  With a SphericalBasis
        for example, the keylst should have the format:

        [[l1, m1, n1], [l2, m2, n2], ...]

        Each sublist represents a PC, where l, m, and n are the 
        spherical harmonics basis parameters.
        )",
	py::arg("config"),
	py::arg("tol"),
	py::arg("count"),
	py::arg("flags") = "");

  g.def("eigenvalues", &KoopmanRKHS::eigenvalues,
	R"(
        Vector of eigenvalues from the EDMD RKHS analysis. 

        Returns
        -------
        numpy.ndarray
            the vector of eigenvalues

        Notes
        -----
        Note that these eigenvalues are complex.  It may be helpful 
        to separate the magnitude and phase using Numpy's 'absolute' 
        and 'angle' functions.
       )");

  g.def("saveState", &KoopmanRKHS::saveState,
	R"(
        Save current EDMD RKHS state to an HDF5 file with the given prefix

        Parameters
        ----------
        prefix : str
            output filename prefix

        Returns
        -------
        None
        )", py::arg("prefix"));

  g.def("restoreState", &KoopmanRKHS::restoreState,
	R"(
        Restore current EDMD RKHS state from an HDF5 file

        Parameters
        ----------
        prefix : str
            input filename prefix

        Notes
        -----
	The Koopman RKHS instance must be constructed with the same data and
        parameters as the saved state.  The restoreState routine will check
        for the same data dimension and trend state but can not sure
        complete consistency.
        )", py::arg("prefix"));

  g.def("getModes", &KoopmanRKHS::getModes,
	R"(
        Get the RKHS mode coefficients for all triples

        Returns
        -------
        numpy.ndarray
            the RKHS mode coefficients

        See also
        --------
        Use in conjunction with 'contributions' to visualize the support 
        from each EDMD RKHS mode to the coefficient series.
        )");

  g.def("modeEval", &KoopmanRKHS::modeEval,
	R"(
        Evaluate the contribution from the index triple

        Parameters
        ----------
        index : int
            the triple index in eigenvalue order
        value : ndarray
            the input point

        Returns
        -------
        numpy.ndarray
            the contribution to the trajectory from the indexed triple

        See also
        --------
        Use in conjunction with 'contributions' to visualize the support 
        from each EDMD RKHS mode to the coefficient series.
        )", py::arg("index"), py::arg("value"));

  g.def("evecEval", &KoopmanRKHS::evecEval,
	R"(
        Evaluate the Koopman eigenfunction with given index

        Parameters
        ----------
        index : int
            the triple index in eigenvalue order
        value : ndarray
            the input point

        Returns
        -------
        numpy.ndarray
            the Koopman eigenfunction

        See also
        --------
        Use in conjunction with 'contributions' to visualize the support 
        from each EDMD RKHS mode to the coefficient series.
        )", py::arg("index"), py::arg("value"));

  g.def("getAllKeys", &KoopmanRKHS::getAllKeys,
	R"(
        Provides a list of all internal channel keys (for reference)

        Returns
        -------
        list(Key)
            list of keys in the format described in config file

        See also
        --------
        KoopmanRKHS
        )");

  py::class_<MSSA::LiouvilleRKHS, std::shared_ptr<MSSA::LiouvilleRKHS>>
    h(m, "LiouvilleRKHS");

  h.def(py::init<const mssaConfig&, double, int, const std::string>(),
	R"(
        Liouville RKHS operator approximatation class

        Parameters
        ----------
        config : mssaConfig
	    the input database of components
        tol : double
            tolerance for eigenvalue decomposition
        count : int
            top count of eigenvalues
	flag : str
            YAML stanza of parameter values

        Returns
        -------
        Liouville RKHS instance

        Notes
        -----
        The configuration should be in the format:

        {'example': (coefs, keylst, [])}

        where keylst is a list of selected PCs.  With a SphericalBasis
        for example, the keylst should have the format:

        [[l1, m1, n1], [l2, m2, n2], ...]

        Each sublist represents a PC, where l, m, and n are the 
        spherical harmonics basis parameters.
        )",
	py::arg("config"),
	py::arg("tol"),
	py::arg("count"),
	py::arg("flags") = "");

  h.def("eigenvalues", &LiouvilleRKHS::eigenvalues,
	R"(
        Vector of eigenvalues from the EDMD RKHS analysis. 

        Returns
        -------
        numpy.ndarray
            the vector of eigenvalues

        Notes
        -----
        Note that these eigenvalues are complex.  It may be helpful 
        to separate the magnitude and phase using Numpy's 'absolute' 
        and 'angle' functions.
       )");

  h.def("saveState", &LiouvilleRKHS::saveState,
	R"(
        Save current EDMD RKHS state to an HDF5 file with the given prefix

        Parameters
        ----------
        prefix : str
            output filename prefix

        Returns
        -------
        None
        )", py::arg("prefix"));

  h.def("restoreState", &LiouvilleRKHS::restoreState,
	R"(
        Restore current EDMD RKHS state from an HDF5 file

        Parameters
        ----------
        prefix : str
            input filename prefix

        Notes
        -----
	The Liouville RKHS instance must be constructed with the same data and
        parameters as the saved state.  The restoreState routine will check
        for the same data dimension and trend state but can not sure
        complete consistency.
        )", py::arg("prefix"));

  h.def("getModes", &LiouvilleRKHS::getModes,
	R"(
        Get the RKHS mode coefficients for all triples

        Returns
        -------
        numpy.ndarray
            the RKHS mode coefficients

        See also
        --------
        Use in conjunction with 'contributions' to visualize the support 
        from each EDMD RKHS mode to the coefficient series.
        )");

  h.def("modeEval", &LiouvilleRKHS::modeEval,
	R"(
        Evaluate the contribution from the index triple

        Parameters
        ----------
        index : int
            the triple index in eigenvalue order
        value : ndarray
            the input point

        Returns
        -------
        numpy.ndarray
            the contribution to the trajectory from the indexed triple

        See also
        --------
        Use in conjunction with 'contributions' to visualize the support 
        from each EDMD RKHS mode to the coefficient series.
        )", py::arg("v"));

  h.def("evecEval", &LiouvilleRKHS::evecEval,
	R"(
        Evaluate the Liouville eigenfunction with given index

        Parameters
        ----------
        index : int
            the triple index in eigenvalue order
        value : ndarray
            the input point

        Returns
        -------
        numpy.ndarray
            the Liouville eigenfunction

        See also
        --------
        Use in conjunction with 'contributions' to visualize the support 
        from each EDMD RKHS mode to the coefficient series.
        )", py::arg("index"), py::arg("v"));

  h.def("basisEval", &LiouvilleRKHS::basisEval,
	R"(
        Return the Liouville basis function

        Parameters
        ----------
        None

        Returns
        -------
        numpy.ndarray
            the Liouville eigenbasis

        See also
        --------
        Use in conjunction with 'contributions' to visualize the support 
        from each EDMD RKHS mode to the coefficient series.
        )");

  h.def("getAllKeys", &LiouvilleRKHS::getAllKeys,
	R"(
        Provides a list of all internal channel keys (for reference)

        Returns
        -------
        list(Key)
            list of keys in the format described in config file

        See also
        --------
        LiouvilleRKHS
        )");

}
