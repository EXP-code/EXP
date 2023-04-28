#ifndef EXP_MSSA_H
#define EXP_MSSA_H

#include <yaml-cpp/yaml.h>
#include "CoefContainer.H"

namespace MSSA
{
  /**
     Class for eDMD processing of EXP coefficients along with
     user-specified auxiliary channels
  */
  class Koopman
  {

  protected:

    //@{
    //! Repacked stream data for eDMD
    std::map<Key, std::vector<double>, mSSAkeyCompare > data;
    //@}

    //! Coefficient container
    CoefContainer coefDB;

    //! Working updated copy of DB
    std::shared_ptr<CoefContainer> newDB;

    //! Parameter database
    YAML::Node params;

    //! Primary Koopman analysis
    void koopman_analysis();

    bool computed, reconstructed;

    //! EDMD modes
    Eigen::MatrixXcd Phi;

    //! DMD state matricies
    Eigen::MatrixXd X0, X1;

    //! Singular values
    Eigen::VectorXd S;

    //! Right singular vectors
    Eigen::MatrixXd U;

    //! Left singular vectors
    Eigen::MatrixXd V;

    //! Koopman matrix approximation
    Eigen::MatrixXd A;

    //! Eigenvalues of Koopman
    Eigen::VectorXcd L;

    //! Eigenvectors of Koopman
    Eigen::MatrixXcd W;

    //! Reconstructed time series
    Eigen::MatrixXd Y;

    //! Parameters
    //@{
    bool verbose, powerf, project;
    std::string prefix, config;
    int nev;
    //@}

    //! Construct YAML node from string
    void assignParameters(const std::string pars);

    //! Number of channels
    int nkeys;

    //! Number of points in the time series
    int numT;

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

  public:

    /** Constructor

	@param spec map/dictionary of tuples listing the Coefs object
	and a list of keys

	@param flags is a string of YAML with changes for the default '
	flag values

	@param window is the the window length
	@param maxEV is the maximum number of eigenvectors

	The map/dictionary has the following structure:
	{
	"mnemonic1": (Coefs1, [ [key11], [key12], [...], ...]),
	"mnemonic2": (Coefs2, [ [key21], [key22], [...], ...]),
	.
	.
	}

	where the mnemonic is choosen for convenience the set of keys
	for each coefficient set, Coefs, specify the indices in the
	dimensionaly specific the the Coefs instance itself.
	E.g. harmonic and radial indicies for spherical and
	cylindrical bases.
    */
    Koopman(const mssaConfig& spec, int maxEV, const std::string flags="");

    //! Destructor
    virtual ~Koopman() {}

    //! Get the eigenvalues
    Eigen::VectorXcd eigenvalues()
    {
      if (not computed) koopman_analysis();
      return L;
    }

    //! Return the EDMD modes
    Eigen::MatrixXcd getModes()
    {
      if (not computed) koopman_analysis();
      return Phi;
    }

    /** Return the frequency and power spectrum of the channel vectors for
	the reconstruction

	See @pcDFT for intended C++ call signature and return.
    */
    std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
    channelDFT();

    /** Reconstruction for desired eigenvalues

	@param evlist is a vector of indices correponding to the
	eigenvalue and PC ordering (i.e. a group). An empty list
	implies no PCs will be used for the reconstruction. That is,
	you will get the mean values or zero values for the 'totVar'
	and the 'totPow' options.
    */
    void reconstruct(const std::vector<int>& evlist);

    /** Get the reconstructed coefficients in an updated Coefs stuctures

	The values are returned as map/dictionary with the mnemonic
	name as the key.  You may call resetReconstructed() to get a
	fresh copy of the data.
    */
    std::map<std::string, CoefClasses::CoefsPtr> getReconstructed();

    /** Compute the diagnostic computation of contributions per channel
	In C++, consider using std::tie to get the unpacked values, e.g.:
	 //
	 Eigen::VectorXd F;
	 Eigen::MatrixXd G;
	 std::tie(F, G) = contributions();
	 // Do something with matrices F and G
    */
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> contributions();

    //! Compute the diagnostic computation of contributions per
    //! channel and make image
    void contributionsPNG();

    //! Save current MSSA state to an HDF5 file with the given prefix
    void saveState(const std::string& prefix);

    //! Restore current MSSA state to an HDF5 file with the given prefix
    void restoreState(const std::string& prefix);

    //! Provides a list of all channel keys
    std::vector<Key> getAllKeys()
    {
      std::vector<Key> ret;
      for (auto v : data) ret.push_back(v.first);
      return ret;
    }
  };


}
// END namespace MSSA

#endif
