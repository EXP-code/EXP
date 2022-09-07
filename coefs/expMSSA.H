#ifndef EXP_MSSA_H
#define EXP_MSSA_H

#include <yaml-cpp/yaml.h>
#include "CoefDB.H"

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
    std::map<Key, std::vector<double> > data;
    std::map<Key, double> mean, var;
    //@}

    //! Coefficient container
    CoefContainer coefDB;

    //! Parameter database
    YAML::Node params;

    //! Primary MSSA analysis
    void mssa_analysis();

    bool computed, reconstructed;

    //! The reconstructed coefficients for each PC
    std::map<Key, Eigen::MatrixXd> RC;

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
    bool chatty, flip, verbose, zeropad, dfiles;
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

    //! Return the power spectrum of the PC vectors
    Eigen::MatrixXd pcDFT(Eigen::VectorXd& F, Eigen::VectorXd& P);

    //! Return the power spectrum of the channel vectors
    Eigen::MatrixXd channelDFT(Eigen::VectorXd& F, Eigen::VectorXd& P);

    /** Reconstruction for desired eigenvalues 

	@param evlist is a vector of indices correponding to the
	eigenvalue and PC ordering (i.e. a group).  All eigenvalues
	will be used for an empty list.
    */
    void reconstruct(const std::vector<int>& evlist);

    /** Get the reconstructed coefficients in an updated Coefs stuctures

	The values are returned as map/dictionary with the mnemonic
	name as the key.

	@param zero controls whether the unused channels are set to (zero=true)
	or the original data is updated with the used channels (zero=false)
    */
    std::map<std::string, Coefs::CoefsPtr> getReconstructed(bool zero);

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

    //! Compute the diagnostic computation of contributions per channel
    void contributions();

    //! Create wcorrlation matricies and output PNG
    void wcorrPNG();

    /**
       Kmean analysis of the trajectories by PC with fixed cluster size
       
       @param toTerm write to stdout if true
       @param toFile write to file if true
    */
    void kmeans(int clusters, bool toTerm=true, bool toFile=false);
  };


}
// END namespace MSSA

#endif
