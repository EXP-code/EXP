#ifndef EXP_MSSA_H
#define EXP_MSSA_H

#include <yaml-cpp/yaml.h>
#include "CoefContainer.H"

namespace MSSA
{
  /**
     Class for MSSA processing of EXP coefficients along with
     user-specified auxiliary channels
  */
  class expMSSA
  {

  protected:
    
    //@{
    //! Repacked stream data
    std::map<Key, std::vector<double>, mSSAkeyCompare > data;
    std::map<Key, double, mSSAkeyCompare> mean, var;
    //@}

    //! Coefficient container
    CoefContainer coefDB;

    //! Working updated copy of DB
    std::shared_ptr<CoefContainer> newDB;

    //! Parameter database
    YAML::Node params;

    //! Primary MSSA analysis
    void mssa_analysis();

    bool computed, reconstructed;

    //! The reconstructed coefficients for each PC
    std::map<Key, Eigen::MatrixXd, mSSAkeyCompare> RC;

    //! The left singular vectors (PC)
    Eigen::MatrixXd PC;

    //! MSSA variables
    Eigen::MatrixXd Y;

    //! Singular values
    Eigen::VectorXd S;

    //! Right singular vectors
    Eigen::MatrixXd U;

    //! Parameters
    //@{
    bool flip, verbose, powerf;
    std::string prefix, config, spec;
    int numW, nmin, nmax, npc;
    double evtol;
    //@}

    //! Construct YAML node from string
    void assignParameters(const std::string pars);

    //! Number of channels
    int nkeys;

    //! Number of points in the time series
    int numT;

    //! Columns in trajectory matrix
    int numK;

    //! Number of components in the reconstruction
    int ncomp;

    //! Normalization values
    double totVar, totPow;

    //! Toggle for detrending
    bool useMean;
    
    //! Detrending type
    TrendType type;

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
    expMSSA(const mssaConfig& spec, int window, int maxEV,
	 const std::string flags="");

    //! Destructor
    virtual ~expMSSA() {}

    //! Get the eigenvalues
    Eigen::VectorXd eigenvalues()
    {
      if (not computed) mssa_analysis();
      return S;
    }

    //! Cumulative sum of the eigenvalues
    Eigen::VectorXd cumulative()
    {
      if (not computed) mssa_analysis();
      auto C = S;
      for (int i=1; i<C.size(); i++) C(i) += C(i-1);
      return C;
    }

    //! Right singular vectors
    Eigen::MatrixXd getU()
    {
      if (not computed) mssa_analysis();
      return U;
    }

    //! Return the PC vectors
    Eigen::MatrixXd getPC()
    {
      if (not computed) mssa_analysis();
      return PC;
    }

    //! Return the number of PC vectors
    int numPC()
    {
      if (not computed) mssa_analysis();
      return PC.cols();
    }

    /**  Return the frequency and power spectrum of the PC vectors
	 In C++, one might use this as follows:
	 //------------------------------------------------------------
	 Eigen::VectorXd F;
	 Eigen::MatrixXd P;
	 std::tie(F, P) = pcDFT();
	 // Do something with F and P.  
	 //------------------------------------------------------------
	 This example is the C++ equivalent of the Python: (F, P) = pcDFT()
    */
    std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
    pcDFT();

    /** Return the frequency and power spectrum of the channel vectors for
	the reconstruction
	
	See @pcDFT for intended C++ call signature and return.
    */
    std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
    channelDFT();

    /** Return the frequency and power spectrum by PC for the selected channel
	
	See @pcDFT for intended C++ call signature and return.
    */
    std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
    singleDFT(const Key& key);

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

    //@{
    /** Compute w-correlation matrix for a single trajectory matrix

	@param name is the nmemonic name (see constructor)
	@param index is one of the specified keys (see constructor)
    */
    Eigen::MatrixXd wCorrKey(const Key& key);
    Eigen::MatrixXd wCorr(const std::string& name, const Key& ckey);
    //@}

    //! Compute w-correlation matrix for all channels
    Eigen::MatrixXd wCorrAll();

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

    //! Create wcorrlation matricies and output PNG
    void wcorrPNG();

    /**
       Kmean analysis of the trajectories by PC with fixed cluster size
       
       @param toTerm write to stdout if true
       @param toFile write to file if true
    */
    void kmeans(int clusters, bool toTerm=true, bool toFile=false);

    //! Save current MSSA state to an HDF5 file with the given prefix
    void saveState(const std::string& prefix);

    //! Restore current MSSA state to an HDF5 file with the given prefix
    void restoreState(const std::string& prefix);

    //! Return total variance value used for normalizing coefficient series
    double getTotVar() { return totVar; }

    //! Return total power value used for normalizing coefficient series
    double getTotPow() { return totPow; }

    //! Provide access to the detrended reconstruction
    Eigen::MatrixXd getRC(const Key& key) {
      auto it = RC.find(key);
      if (it == RC.end())
	throw std::runtime_error("expMSSA:getRC: key not found");
      return it->second;
    }

    //! Provides a list of keys for accessing the detrended
    //! reconstruction using getRC()
    std::vector<Key> getRCkeys()
    {
      std::vector<Key> ret;
      for (auto v : RC) ret.push_back(v.first);
      return ret;
    }

    //! Provides a list of all channel keys
    std::vector<Key> getAllKeys()
    {
      std::vector<Key> ret;
      for (auto v : mean) ret.push_back(v.first);
      return ret;
    }
  };


}
// END namespace MSSA

#endif
