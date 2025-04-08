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
    "  verbose: false        Report on internal progress for debugging only\n"
    "  noMean: false         If true, do not subtract the mean when\n"
    "                        reading in channels. Used only for totPow\n"
    "                        detrending method.\n"
    "  writeCov: true        Write the covariance matrix to a file for\n"
    "                        diagnostics since this is not directly\n"
    "                        available from the interface. Used only if\n"
    "                        Traj: false.\n"
    "  Jacobi: true          Use the Jacobi SVD rather than the Random\n"
    "                        approximation algorithm from Halko, Martinsson,\n"
    "                        and Tropp (RedSVD). This is quite accurate but\n"
    "                        _very_ slow\n"
    "  BDCSVD: true          Use the Bidiagonal Divide and Conquer SVD\n"
    "                        rather than RedSVD; this is faster and more\n"
    "                        accurate than the default RedSVD but slower.\n"
    "  Traj: true            Perform the SVD of the trajectory matrix\n"
    "                        rather than the more computationally intensive\n"
    "                        but more stable SVD of the covariance matrix.\n"
    "                        Set to false to get standard covariance SVD.\n"
    "                        In practice, 'Traj: true' is sufficiently\n"
    "                        accurate for all PC orders that are typically\n"
    "                        found to have dynamical signal.\n"
    "  rank: 100             The default rank for the randomized matrix SVD.\n"
    "                        The default value will give decent accuracy with\n"
    "                        small computational overhead and will be a good\n"
    "                        choice for most applications.\n"
    "  RedSym: true          Use the randomized symmetric eigenvalue solver\n"
    "                        RedSym rather rather than RedSVD for the co-\n"
    "                        variance matrix SVD (Traj: false). The main use\n"
    "                        for this is checking the accuracy of the default\n"
    "                        randomized matrix methods.\n"
    "  allchan: true         Perform k-means clustering analysis using all\n"
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
        The MSSA analysis class

        Parameters
        ----------
        config : mssaConfig
             the input database of components; see Notes. 
        window : int
             the length of the cross-correlation interval. Start with
             half of the time series length.
        numpc : int
             the default number of eigenvalues to compute
        flags : str, default=""
             YAML stanza of parameter values

        Returns
        -------
        None

        Notes
        -----
        The configuration should be in the format:

        {'example': (coefs, keylst, []), ...}

        where keylst is a list of selected PCs.  For example, a SphericalBasis 
        would have keys in the format: 

        [[l1, m1, n1], [l2, m2, n2], ...]

        Each sublist represents a PC, where l, m, and n are the spherical 
        harmonics basis parameters.
        )",
	py::arg("config"),
	py::arg("window"),
	py::arg("numpc"),
	py::arg("flags") = "");


  f.def("eigenvalues", &expMSSA::eigenvalues,
    R"(
    Return the vector of eigenvalues from the MSSA analysis

    Returns
    -------
    list(float)
        vector of eigenvalues obtained from the MSSA analysis
    )");

  f.def("cumulative", &expMSSA::cumulative,
    R"(
    Cumulatively summed vector of eigenvalues from the MSSA analysis

    Returns
    -------
    list (float) 
        cumulatively summed vector
    )");

  
  f.def("getU", &expMSSA::getU,
	R"(
        right-singular vectors from the MSSA analysis

        These vectors describe the contribution of each channel to each PC

        Returns
        -------
        numpy.ndarray
            right-singular vectors
        )");


  f.def("getPC", &expMSSA::getPC,
	R"(
        left-singular vectors from the MSSA analysis

        These are the 'principal component vectors' that 
        describe the correlated temporal variation

        Returns
        -------
        numpy.ndarray
            principal component vectors
    )");

  f.def("pcDFT", &expMSSA::pcDFT,
	R"(
        DFT of the principal component vectors

        Returns
        -------
        tuple(numpy.ndarray, numpy.ndarray)
            The frequency vector and power distribution for each pc

        Notes
        -----
        frequency values are given in angular frequency with 2 * np.pi / T
        )");

  f.def("channelDFT", &expMSSA::channelDFT,
	R"(
        DFT of the selected data channels

        Returns
        -------
        tuple(numpy.ndarray, numpy.ndarray)
            The frequency vector and power distribution for each channel

        Notes
        -----
        Useful for comparing with the power of the principal components (PCs)
        )");
    
  f.def("singleDFT", &expMSSA::singleDFT, py::arg("key"),
	R"(
        DFT of the selected data channel with partial power for each PC

        Parameters
        ----------
        key : list(int)
            identifier indices of the selected data channel

        Returns
        -------
        tuple(numpy.ndarray, numpy.ndarray)
            The frequency vector and power distribution for each channel
        )");
  

  f.def("reconstruct", &expMSSA::reconstruct, py::arg("evlist"),
	R"(
        Reconstruct the data channels with the provided group list of eigenvalue indices

        Parameters
        ----------
        evlist : list(int)
            list of eigenvalue indices specifying the group

        Returns
        -------
        None

        Notes
        -----
        Use getReconstructed() to copy the reconstruction back to the coefficient db

        See also
        --------
        getReconstructed
        )");


  f.def("getReconstructed", &expMSSA::getReconstructed, py::arg("reconstructmean")=true,
	R"(
        provide reconstructed time series

        Parameters
        ----------
        reconstructmean : bool, default=True
             include the detrended mean in the data channel reconstruction

        Returns
        -------
        dict({id: Coefs},...)
             reconstructed time series in the original coefficient form

        Notes
        -----
        The reconstructed data will overwrite the memory of the original coefficient 
        data.  We suggest a deepcopy() if you wish to preserve the input coefficient db.

        Setting 'reconstructmean=False' may be helpful to compare the variance between 
        data channels.
        )");

  f.def("background", &expMSSA::background,
	R"(
        Copy the background data streams back to the working coefficient database

        Notes
        -----
        This can be used after a zerodata() call to include the background in the 
        reconstruction
        )");

  f.def("wCorr", &expMSSA::wCorr, py::arg("name"), py::arg("key"),
	py::arg("nPC") = std::numeric_limits<int>::max(),
	R"(
        The w-correlation matrix for the selected component and channel key

        Parameters
        ----------
        name : str
            The name of the selected component.
        key : tuple
            The identifier key of the selected channel.
        nPC : int
            The maximum rank for reconstruction

        Returns
        -------
        numpy.ndarray
            The w-correlation matrix for the selected component and channel

        Notes
        -----
        The w-correlation matrix needs the reconstructed trajectory matrices for each
        of the eigenvalue, PC pairs.  Calling this method will recompute the reconstruction
        for all eigenvalues up to 'nPC' and return an nPC x nPC matrix.  If the 'nPC'
        parameter is not specified, it will be set to `numpc` used to construct the 
        instance.  Any previous reconstruction will be overwritten.

        Returns the combined cosine+sine correlation for complex types for viewing 
        (e.g., with 'imshow').

        The w-correlation values range from 0 to 1, where a higher value corresponds 
        to a stronger correlation.
        )");

  f.def("wCorrKey", &expMSSA::wCorrKey, py::arg("key"),
	py::arg("nPC") = std::numeric_limits<int>::max(),
	R"(
        Get the w-correlation matrix for the selected component and channel key

        Parameters
        ----------
        key : list(int,...)
            identifier of the selected channel

        Returns
        -------
        numpy.ndarray
            w-correlation matrix for the selected component and channel key.

        Notes
        -----
        The index key here is 'extended' by the prefixed component index.
  
        Computation of the w-correlation matrix needs the reconstructed
        trajectory matrices for each of the (eigenvalue, PC) pairs.  Calling
        this method will recompute the reconstruction for all eigenvalues up to
        order 'npc' and return an (nPC x nPC) matrix.  If the 'nPC' parameter is
        not specified, it will be set to the `numpc` used in the original
        construction.  Any prior reconstruction will be overwritten.

        The rows and columns contain distinct cosine and sine indicies if the channel 
        is complex valued.

        This matrix can be visualized using 'imshow' for plotting.

        The w-correlation values range from 0 to 1, where a higher value corresponds to 
        a stronger correlation.
    )");

  f.def("wCorrAll", &expMSSA::wCorrAll,
	py::arg("nPC") = std::numeric_limits<int>::max(),

	R"(
        the w-correlation matrix for all channels in the reconstruction

        Parameters
        ----------
        nPC : int
            The maximum rank for reconstruction

        Returns
        -------
        numpy.ndarray
            w-correlation matrix for all channels in the reconstruction

        Notes
        -----
        The w-correlation values range from 0 to 1, where a higher value
        corresponds to a stronger correlation.

        Computation of the w-correlation matrix needs the reconstructed
        trajectory matrices for each of the (eigenvalue, PC) pairs.  Calling
        this method will recompute the reconstruction for all eigenvalues up to
        order 'npc' and return an (nPC x nPC) matrix.  If the 'nPC' parameter is
        not specified, it will be set to the `numpc` used in the original
        construction.  Any prior reconstruction will be overwritten.

        See also
        --------
        wCorr : return the w-correlation matrix by component name and key
        wCorrPNG : output w-correlation matrices in PNG images
        wCorrKey : return the w-correlation matrix by extended key
        )");

  f.def("wcorrPNG", &expMSSA::wcorrPNG,
	py::arg("nPC") = std::numeric_limits<int>::max(),
	R"(
        w-correlation matrices and output PNG image representations

        Parameters
        ----------
        nPC : int
            The maximum rank for reconstruction

        Notes
        -----
        The w-correlation values range from 0 to 1, where a higher value
        corresponds to a stronger correlation.

        Computation of the w-correlation matrix needs the reconstructed
        trajectory matrices for each of the (eigenvalue, PC) pairs.  Calling
        this method will recompute the reconstruction for all eigenvalues up to
        order 'npc' and return an (nPC x nPC) matrix.  If the 'nPC' parameter is
        not specified, it will be set to the `numpc` used in the original
        construction.  Any prior reconstruction will be overwritten.

        See also
        --------
        wCorr : return the w-correlation matrix by component name and key
        wCorrAll : return the combined correlation matrix for all components and keys
        wCorrKey : return the w-correlation matrix by extended key
        )");

  f.def("kmeans", &expMSSA::kmeans,
	py::arg("clusters") = 4,
	py::arg("stride") = 2,
	R"(
        Do a k-means analysis on the reconstructed trajectory matrices for a
        single channel (specified key value) to provide grouping insight.  A
        vector of channel indices that identify clusters is return in a vector
        ordered by PC index.

        Parameters
        ----------
        clusters : int, default=4
            number of clusters for the k-means analysis
        stride : int, default=2
            if positive, the initial cluster centers are stride selected from the PC list.
            If zero, the centers are selected randomly from the PC list

        Returns
        -------
        tuple : (numpy.niarray, numpy.ndarray, double)
            The PC indices of the k-means clusters, distance from the centroid, and the
            final update error.  A zero update error implies that the k-means algorithm
            converged.

        Notes
        -----
        The k-means partitions n vector observations into k clusters in which
        each observation belongs to the cluster with the nearest centers while
        minimizing the variance within each cluster.  In this case, the vectors
        are the full trajectory matrices and the distance is the distance
        between the trajectory matricies reconstructed from each eigentriple
        from mSSA.  The distance used here is the Frobenius distance or matrix
        norm distance: the square root of the sum of squares of all elements in
        the difference between two matrices.

        This version does the analysis for all channels together, the most
        useful for estimating groups.  For individual contributions by channel,
        use kmeansChannel.
        )");

  f.def("kmeansChannel", &expMSSA::kmeansChannel,
	py::arg("key"),
	py::arg("clusters") = 4,
	py::arg("stride") = 2,
	R"(
        Do a k-means analysis on the reconstructed trajectory matrices for a
        single channel (specified key value) to provide grouping insight.  In
        most cases, you will want to use the kmeans() version which analyzes all
        channels together.

        Parameters
        ----------
        clusters : int, default=4
            number of clusters for the k-means analysis
        key : list(int)
            identifier indices of the selected data channel

        Returns
        -------
        tuple : (numpy.niarray, numpy.ndarray, double)
            The PC indices of the k-means clusters, distance from the centroid, and the
            final update error.  A zero update error implies that the k-means algorithm
            converged.

        Notes
        -----
        This version does the analysis channel-by-channel.  You may wish to see all channels together using
        kmeansTotal.  See kemans() for more details.
        )");

  f.def("contrib", &expMSSA::contributions,
	R"(
        the relative contribution of each principal component (PC)

        Returns
        -------
        tuple(numpy.ndarray, numpy.ndarray)

        Notes
        -----
        The reutrn tuple (F, G) represents two views of the contributions as follows
            - F: Each PC's contribution to each channel. The columns are L2 normed.
            - G: Each channel's contribution to each PC. The rows are L2 normed.

        By default, channels for non-zero 'm' are split into cosine and sine
        components from the real+imaginary values.

        The L2 norm, or Euclidean norm, computes the length of a vector in a
        multi-dimensional space.  For a vector v = [v1, v2, ..., vn], the L2
        norm is calculated as sqrt(v1^2 + v2^2 + ... + vn^2).

        The L2 normed views provide a measure of the relative contribution of
        each PC to each channel and the relative contribution of each channel to
        each PC. These contributions can be plotted using 'imshow'.
        )");

  f.def("saveState", &expMSSA::saveState,
	R"(
        Save the current MSSA state to an HDF5 file

        Parameters
        ----------
        prefix : str
            prefix used for the HDF5 file

        Returns
        -------
        None
        )", py::arg("prefix"));

  f.def("restoreState", &expMSSA::restoreState,
	R"(
        Restore the current MSSA state from an HDF5 file

        Parameters
        ----------
        prefix : str
             prefix used for the HDF5 file

        Returns
        -------
        None

        Notes
        -----
        The expMSSA instance must be constructed with the same data and parameters as 
        the saved state. The restoreState routine will check for the same data dimension 
        and trend state but cannot ensure complete consistency.
        )", py::arg("prefix"));


  f.def("getTotVar", &expMSSA::getTotVar,
	R"(
	variance value used for normalizing the coefficient series

	Returns
        ------
	float
            variance value
	)");

  f.def("getTotPow", &expMSSA::getTotPow,
	R"(
        power value used for normalizing the coefficient series

	Returns
        -------
        float
	    power value
	)");


  f.def("getKoopmanModes", &expMSSA::getKoopmanModes,
	R"(
        Compute the Koopman mode estimate from the right-singular vectors

        Uses eDMD to estimate the modes

        Parameters
        ----------
        tol : double
            singular value truncation level
        window: int
            Smoothing between serialized channels (0 for no smoothing)
        debug : bool
            flag indicating whether to print debug information

        Notes
        -----
        Use getReconstructedKoopman() to copy the reconstruction for a
        particular mode back to the coefficient db

        Returns
        -------
        tuple(numpy.ndarray, numpy.ndarray)
            vector of eigenvalues and modes
    )", py::arg("tol")=1.0e-12, py::arg("window")=0, py::arg("debug")=false);

  f.def("getReconstructedKoopman", &expMSSA::getReconstructedKoopman,
	R"(
        Reconstruct the coefficients for a particular Koopman mode

        Parameters
        ----------
        mode: int
            The index of the mode to be reconstructed

        Returns
        -------
        dict({id: Coefs},...)
             reconstructed time series in the original coefficient form

    )", py::arg("mode"));

  f.def("getRC", &expMSSA::getRC,
	R"(
        Access the detrended reconstructed channel series by internal key

        Parameters
        ----------
        key : list(int,...)
            internal key for the desired channel

        Returns
        -------
        numpy.ndarray: 
            detrended reconstructed channel series

        Notes
        -----
        keys are a lists of integer values
        )", py::arg("key"));

  f.def("getRCkeys", &expMSSA::getRCkeys,
	R"(
	Provides a list of internal keys for accessing the detrended channel series using getRC().

	Returns
        -------
	list(list(int,...)) 
            list of internal keys representing the detrended channel series

        Notes
        -----
        keys are a lists of integer values
	)");

  f.def("getAllKeys", &expMSSA::getAllKeys,
	R"(
	Provides a list of all internal channel keys for reference

	Returns
        -------
	list(list)
            list of all internal channel keys which are list(int,...)

        Notes
        -----
        keys are a lists of integer values
	)");
}
