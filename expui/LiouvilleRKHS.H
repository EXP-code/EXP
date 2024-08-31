#ifndef EXP_LIOUVILLE_RKHS_H
#define EXP_LIOUVILLE_RKHS_H

#include <yaml-cpp/yaml.h>
#include "CoefContainer.H"

#include <map>

namespace MSSA
{
  /**
     Class for eDMD processing of trajectories using Liouville opertor
     RKHS.  This implements both Algorithms 6.1 and 7.1 from the paper
     by Joel A. Rosenfeld and Rushikesh Kamalapurkar, "Singular
     Dynamic Mode Decomposition", 2023, SIAM J. Appl. Dynamical
     Systems, Vol. 22, No. 3, pp 2357-2381.
  */
  class LiouvilleRKHS
  {

  public:

    //! RKHS types
    enum class RKHS
    {
      Polynomial,
      Exponential,
      Gaussian
    };
   
    //! Analysis method
    enum class Method
    {
      Singular,
      Eigenfunction
    };
   
    //! For RKHS type reflection and parsing
    static std::map<RKHS, std::string> RKHS_names;
    static std::map<std::string, RKHS> RKHS_values;

    //! For method type reflection and parsing
    static std::map<Method, std::string> Method_names;
    static std::map<std::string, Method> Method_values;

  protected:

    //@{
    //! Repacked stream data for eDMD
    std::map<Key, Eigen::MatrixXd, mSSAkeyCompare > data;
    //@}

    //! Coefficient container
    CoefContainer coefDB;

    //! Working updated copy of DB
    std::shared_ptr<CoefContainer> newDB;

    //! Parameter database
    YAML::Node params;

    //! Method
    Method method = Method::Eigenfunction;

    //! Primary Liouville analysis
    void analysis()
    {
      switch (method)
	{
	case Method::Singular:
	  singular_analysis();
	  break;
	case Method::Eigenfunction:
	  eigenfunction_analysis();
	  break;
	}
    }

    //! Singular Louville analysis
    void singular_analysis();

    //! Eigenfunction Liouville analysis
    void eigenfunction_analysis();

    //! Compute G matrix for RKHS space defined by mu value
    Eigen::MatrixXd computeGrammian(double mu, double rat=1.0);

    //! Compute occupation-kernel difference matrix
    Eigen::MatrixXd computeGammaDiff(double mu);

    //! Compute occupation-kernel vector
    Eigen::VectorXd computeGamma(const Eigen::VectorXd& x, double mu);

    //! Compute singular Liouville RKHS matrix
    Eigen::MatrixXd computeGramGamma(double mu);

    //! Compute singular Liouville trajectory matrix
    Eigen::MatrixXd trajectory();

    bool computed, reconstructed;

    //@{
    //! Parameters
    double mu1=1.0, mu2=1.1, eps=1.0e-4;
    bool use_red = true;
    //@}
    
    //! DMD state matrices
    Eigen::MatrixXd G1, G2, G3;

    //! Occupation RKHS matrix
    Eigen::MatrixXd A;

    //! Occupation kernel matrix
    Eigen::MatrixXd O;

    //! Liouville modes for eigenfunction method
    Eigen::MatrixXcd Xi;

    //! Liouville modes for singular method
    Eigen::MatrixXd Xh;

    //! Eigenbasis
    Eigen::MatrixXcd V;

    //! Normalized eigenbasis
    Eigen::MatrixXcd Vbar;

    //! Eigenvectors
    Eigen::VectorXcd L;

    //! For projection of initial state into eigenfunctions
    Eigen::MatrixXcd Ginv;

    //! EDMD modes
    Eigen::MatrixXcd Phi;

    //! Eigenvalue tolerance
    double tol = 1.0e-6;

    //! Max eigenvalue count
    int evCount = 100;

    //! Eigenvalues values
    Eigen::VectorXd S1, S2, S3;

    //! Eigenvectors
    Eigen::MatrixXd Q1, Q2, Q3;

    //@{
    //! SVD vectors and normalizations
    Eigen::MatrixXd UU, VV, V0, Vt;
    Eigen::VectorXd SS, Dq, Dp;
    //@}

    //! Parameters
    //@{
    bool verbose, powerf, project;
    std::string prefix, config;
    int nev;
    //@}

    //! Construct YAML node from string
    void assignParameters(const std::string pars);

    //! Number of trajectories
    int traj;

    //! Number of active trajectories
    int nkeys;

    //! Phase-space rank
    int rank;

    //! Number of points in the time series
    int numT;

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! RKHS parameters
    double d, alpha;

    //! RKHS type
    RKHS rkhs = RKHS::Polynomial;

    //! RKHS kernels
    double kernel(const Eigen::VectorXd& x,
		  const Eigen::VectorXd& y,
		  double mu);

    //! Occupation kernel matrix
    Eigen::MatrixXd occupation();

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
    LiouvilleRKHS(const mssaConfig& spec, double tol, int count, const std::string flags="");

    //! Destructor
    virtual ~LiouvilleRKHS() {}

    //! Get the eigenvalues
    Eigen::VectorXcd eigenvalues()
    {
      if (not computed) analysis();
      if (method == Method::Eigenfunction)
	return L;
      else
	return Eigen::VectorXcd(SS);
    }

    //! Return the Liouville modes, the coefficients to the eigenfunctions
    Eigen::MatrixXcd modeEval()
    {
      if (method != Method::Eigenfunction)
	throw std::runtime_error("LiouvilleRKHS::modeEval is only for the Eigenfunction method");
      
      if (not computed) analysis();

      return Xi;
    }

    //! Return the eigenfunctions
    Eigen::VectorXcd evecEval(const Eigen::VectorXd& x);

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

    //! Switch method: Singular or Eigenfunction
    void setMethod(const std::string& method_name)
    {
      auto it = Method_values.find(method_name);
      if (it != Method_values.end()) method = it->second;
      else throw std::runtime_error("LiouvilleRKHS::setMethod: unknown method <" + method_name + ">.  Expecting <Singular> or <Eigenfunction>");
    }

    //! Flow vector evaluated at a point
    Eigen::VectorXd flow(const Eigen::VectorXd& x);

    //! Compute trajectory for an initial state using RK4 (for singular method only)
    Eigen::MatrixXd computeTrajectory(const Eigen::VectorXd& x, double h);

  };

}
// END namespace MSSA

#endif
