#ifndef EXP_LIOUVILLE_RKHS_H
#define EXP_LIOUVILLE_RKHS_H

#include <yaml-cpp/yaml.h>
#include "CoefContainer.H"

#include <map>

namespace MSSA
{
  /**
     Class for eDMD processing of trajectories using Liouville opertor
     RKHS
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
   
    //! For type reflection and parsing
    static std::map<RKHS, std::string> RKHS_names;
    static std::map<std::string, RKHS> RKHS_values;

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

    //! Primary Koopman/Liouville analysis
    void analysis();

    //! Compute G matrix for RKHS space defined by mu value
    Eigen::MatrixXd computeGrammian(double mu, double rat=1.0);

    //! Compute occupation-kernel difference matrix
    Eigen::MatrixXd computeGammaDiff(double mu);

    bool computed, reconstructed;

    //@{
    //! Parameters
    double mu1=1.0, mu2=1.1, eps=1.0e-4;
    bool use_red = true;
    //@}
    
    //! DMD state matrices
    Eigen::MatrixXd G1, G2, G3;

    //! Occupation matrix
    Eigen::MatrixXd A;

    //! Eigenvectors
    Eigen::VectorXcd L;

    //! Eigenfunctions
    Eigen::MatrixXcd Xi;

    //! EDMD modes
    Eigen::MatrixXcd Phi;

    //! Eigenvalue tolerance
    double tol = 1.0e-6;

    //! Max eigenvalue count
    int evCount = 100;

    //! Eigenvalues values
    Eigen::VectorXd S2, S3;

    //! Eigenvectors
    Eigen::MatrixXd Q2, Q3;

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
      return L;
    }

    //! Return the Liouville modes, the coefficients to the eigenfunctions
    Eigen::MatrixXcd getModes()
    {
      if (not computed) analysis();
      return Xi;
    }

    //! Evaluate the contribution from the indexed triple
    Eigen::VectorXcd modeEval(int index, const Eigen::VectorXd& x);

    //! Return the eigenfunction
    std::complex<double> evecEval(int index, const Eigen::VectorXd& x);

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
