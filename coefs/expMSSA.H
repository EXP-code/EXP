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
    std::map<MSSA::Key, std::vector<double> > data;
    std::map<MSSA::Key, double> mean, var;
    //@}

    //! Coefficient container
    CoefContainer coefDB;

    //! Parameter database
    YAML::Node params;

    //! Primary MSSA analysis
    void mssa_analysis();

    bool computed, reconstructed;

    //! The reconstructed coefficients for each PC
    std::map<MSSA::Key, Eigen::MatrixXd> RC;

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
    std::string prefix, grpfile, config, spec;
    bool quiet, flip, verbose, zeropad, dfiles;
    int numW, nmin, nmax, npc;
    int stride, skip, clusters;
    std::vector<int> MM;
    std::vector<std::vector<int>> groups;
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
	eigenvalue and PC ordering
    */
    void reconstruct(const std::vector<int>& evlist);

    /** Reconstruction for desired eigenvalues

	@param evlist is a vector of indices correponding to the
	eigenvalue and PC ordering
	@param type is the reconstruction style
    */
    std::map<std::string, Coefs::CoefsPtr> getReconstructed();

    //@{
    /** Compute w-correlation matrix for a single trajectory matrix

	@param name is the nmemonic name (see constructor)
	@param index is one of the specified keys (see constructor)
    */
    Eigen::MatrixXd wCorrKey(const std::vector<unsigned>& key);
    Eigen::MatrixXd wCorr(const std::string& name,
			  const std::vector<unsigned>& ckey)
    {
      auto key = ckey;
      key.push_back(coefDB.index(name));
      return wCorrKey(key);
    }
    //@}

    //! Compute w-correlation matrix for all channels
    Eigen::MatrixXd wCorrAll();

    //! Compute the diagnostic computation of contributions per channel
    void contributions();

    //! Create wcorrlation matricies and output PNG
    void wcorrPNG();

    //! Kmean analysis of the trajectories by PC
    void kmeans();

    //! Detrend channels
    void channelReconstruct();
  };


}

#endif
