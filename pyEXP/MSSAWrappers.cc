#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <expMSSA.H>

namespace py = pybind11;
#include <TensorToArray.H>

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
    "  verbose: false        Whether there is report or not\n"
    "  noMean: false         If true, do not subtract the mean when\n"
    "                        reading in channels. Valid only for totPow\n"
    "                        detrending method.\n"
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
    "  power: false          Compute and output Fourier power\n"
    "  totVar: false         Detrend according to the total variance\n"
    "                        in all channel\n"
    "  totPow: false         Detrend according to the total power in\n"
    "                        all channels\n\n"
    "The following parameters take values,\ndefaults are given in ()\n\n"
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

  f.def(py::init<const mssaConfig&, int, int, const std::string&>(),
    R"(
    Constructor.

    Args:
        config (mssaConfig): The input database of components. The configuration should be in the format:
                             {'example': (coefs, keylst, [])}
                             where keylst is a list of selected PCs in the format: 
                             [[l1, m1, n1], [l2, m2, n2], ...]
                             Each sublist represents a PC, where l, m, and n are the spherical harmonics basis parameters.
        window (int): The length of the cross-correlation interval. It is suggested to use half of the time series length.
        numpc (int): The default number of eigenvalues to compute.
        flags (str): A YAML stanza of parameter values. Default is an empty string.

    Returns:
        None
    )",
    py::arg("config"),
    py::arg("window"),
    py::arg("numpc"),
    py::arg("flags") = "");


  f.def("eigenvalues", &expMSSA::eigenvalues,
    R"(
    Return the vector of eigenvalues from the MSSA analysis.

    Returns:
        list: The vector of eigenvalues obtained from the MSSA analysis.
    )");

  f.def("cumulative", &expMSSA::cumulative,
    R"(
    Return a cumulatively summed vector of eigenvalues from the MSSA analysis.

    Returns:
        list: A cumulatively summed vector of eigenvalues from the MSSA analysis.
    )");

  
  f.def("getU", &expMSSA::getU,
    R"(
    Return the right-singular vectors from the MSSA analysis
    which describe the contribution of each channel to each PC.

    Returns:
        ndarray: The right-singular vectors describing the contribution
        of each channel to each principal component.
    )");


  f.def("getPC", &expMSSA::getPC,
    R"(
    Return the principal component (left-singular) vectors from the MSSA
    analysis which describe the key temporal variation.

    Returns:
        ndarray: The principal component vectors describing the key temporal variation.
    )");

  f.def("pcDFT", &expMSSA::pcDFT,
    R"(
    Return the DFT of the principal component vectors for quantifying temporal power distribution.

    The frequency values are given in terms of angular frequency with 2 * np.pi / T.

    Returns:
        ndarray: The DFT of the principal component vectors for temporal power distribution.
    )");

  f.def("channelDFT", &expMSSA::channelDFT,
    R"(
    Returns the frequency values and the DFT of the selected data channels.

    This information can be used to compare with the power of the principal components (PC).

    Returns:
        tuple: A tuple containing the angular frequency values and the DFT of the selected data channels.
    )");
    
  f.def("singleDFT", &expMSSA::singleDFT, py::arg("key"),
    R"(
    Returns the frequency, the DFT of the selected data channel
    with partial power for each PC.

    Args:
        key (str): The identifier of the selected data channel.

    Returns:
        tuple: A tuple containing the frequency values and the DFT of
        the selected data channel with partial power for each PC.
    )");
  

  f.def("reconstruct", &expMSSA::reconstruct, py::arg("evlist"),
    R"(
    Reconstruct the data channels with the provided list of eigenvalue indices (a group).

    Args:
        evlist (list): The list of eigenvalue indices specifying the group.

    Returns:
        ndarray: The reconstructed data channels.
    )");


  f.def("getReconstructed", &expMSSA::getReconstructed, py::arg("reconstructmean")=true,
    R"(
    Return the reconstructed time series in the original coefficient form
    that may be used in basis classes.   Setting 'reconstructmean=False' may 
    be helpful to compare the variance between data channels.

    Note: The reconstructed data will overwrite the memory of the original 
          coefficient data.

    Args:
        reconstructmean (bool): If True (default), the data channel includes 
                                the mean value that was subtracted during 
                                SSA detrending.

    Returns:
        ndarray: The reconstructed time series in the original coefficient form.
    )");

  f.def("background", &expMSSA::background,
    R"(
    Copy the background data streams back to the working coefficient
    database.

    This can be used after a zerodata() call to include the background in the reconstruction.
    )");

  f.def("wCorr", &expMSSA::wCorr, py::arg("name"), py::arg("key"),
    R"(
    Get the w-correlation matrix for the selected component and channel key.

    Returns the combined cosine+sine correlation for complex types for viewing (e.g., with 'imshow').

    The w-correlation values range from 0 to 1, where a higher value corresponds to a stronger correlation.

    Args:
        name (str): The name of the selected component.
        key (str): The identifier of the selected channel.

    Returns:
        ndarray: The w-correlation matrix for the selected component and channel.
    )");

  f.def("wCorrKey", &expMSSA::wCorrKey, py::arg("key"),
    R"(
    Get the w-correlation matrix for the selected component and channel key,
    extended by the cosine/sine index if the channel is complex and the component index.

    This matrix can be visualized using 'imshow' for plotting.

    The w-correlation values range from 0 to 1, where a higher value corresponds to a stronger correlation.

    Args:
        key (str): The identifier of the selected channel.

    Returns:
        ndarray: The w-correlation matrix for the selected component and channel key.
    )");

  f.def("wCorrAll", &expMSSA::wCorrAll,
    R"(
    Get the w-correlation matrix for all channels in the reconstruction.

    These matrices can be nicely plotted using 'imshow'.

    The w-correlation values range from 0 to 1, where a higher value corresponds to a stronger correlation.

    Returns:
        ndarray: The w-correlation matrices for all channels in the reconstruction.
    )");

  f.def("wcorrPNG", &expMSSA::wcorrPNG,
    R"(
    Create w-correlation matrices and output PNG image representations.

    The w-correlation values range from 0 to 1, where a higher value corresponds to a stronger correlation.
    )");

  f.def("kmeans", &expMSSA::kmeans, py::arg("clusters") = 4, py::arg("toTerm") = true, py::arg("toFile") = false,
    R"(
    Perform a k-means analysis on the reconstructed trajectory matrices to provide grouping insight.

    By default, this will write the output to the standard output. Set `toFile=True` to write the output to a file.
    The file name will be derived from the 'output' parameter.

    Note: The output format and how to interpret it is currently work in progress. [TODO: Update explanation]

    Args:
        clusters (int): The number of clusters for the k-means analysis.
        toTerm (bool): Flag indicating whether to output to the terminal (standard output).
        toFile (bool): Flag indicating whether to write the output to a file.

    Returns:
        None
    )");

  f.def("contrib", &expMSSA::contributions,
    R"(
    Computes the relative contribution of each principal component (PC) to the coefficient series
    and the breakdown of the coefficient series to each PC.

    The views are L2 normed on columns and rows and returned as a 2D tuple.

    Returns:
        tuple: A 2D tuple (F, G) representing the contributions:
            - F: Each PC's contribution to each channel. The columns are L2 normed.
            - G: Each channel's contribution to each PC. The rows are L2 normed.

    By default, channels for non-zero 'm' are split into cosine and sine components from the real+imaginary values.

    The L2 norm, or Euclidean norm, computes the length of a vector in a multi-dimensional space.
    For a vector v = [v1, v2, ..., vn], the L2 norm is calculated as sqrt(v1^2 + v2^2 + ... + vn^2).

    The L2 normed views provide a measure of the relative contribution of each PC to each channel
    and the relative contribution of each channel to each PC. These contributions can be plotted using 'imshow'.
    )");


  f.def("saveState", &expMSSA::saveState,
    R"(
    Save the current MSSA state to an HDF5 file with the given prefix.

    Args:
        prefix (str): The prefix used for the HDF5 file.

    Returns:
        None
    )", py::arg("prefix"));

  f.def("restoreState", &expMSSA::restoreState,
    R"(
    Restore the current MSSA state from an HDF5 file with the given prefix.

    Note: To use this method, the expMSSA instance must be constructed with the same
    data and parameters as the saved state. The restoreState routine will check for
    the same data dimension and trend state but cannot ensure complete consistency.

    Args:
        prefix (str): The prefix used for the HDF5 file.

    Returns:
        None
    )", py::arg("prefix"));


  f.def("getTotVar", &expMSSA::getTotVar,
	R"(
	Returns the variance value used for normalizing the coefficient series.

	Returns:
	    The variance value.
	)");

  f.def("getTotPow", &expMSSA::getTotPow,
	R"(
	Returns the power value used for normalizing the coefficient series.

	Returns:
	    The power value.
	)");


  f.def("getRC", &expMSSA::getRC,
    R"(
    Access the detrended reconstructed channel series by internal key.

    Args:
        key (str): The internal key for the desired channel.

    Returns:
        ndarray: The detrended reconstructed channel series.
    )", py::arg("key"));

  f.def("getRCkeys", &expMSSA::getRCkeys,
	R"(
	Provides a list of internal keys for accessing the detrended channel series using getRC().

	Returns:
	    A list of internal keys representing the detrended channel series.
	)");

  f.def("getAllKeys", &expMSSA::getAllKeys,
	R"(
	Provides a list of all internal channel keys for reference.

	Returns:
	    A list of all internal channel keys.
	)");


}
