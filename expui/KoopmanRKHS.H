#ifndef EXP_KOOPMAN_RKHS_H
#define EXP_KOOPMAN_RKHS_H

#include <yaml-cpp/yaml.h>
#include "CoefContainer.H"

#include <map>

namespace MSSA
{
  /**
     Class for eDMD processing of EXP coefficients along with
     user-specified auxiliary channels
  */
  class KoopmanRKHS
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

    //! Primary Koopman analysis
    void analysis();

    bool computed, reconstructed;

    //@{
    //! Test stuff
    double lam=0.9, mu=0.5, c=0.0, tscale=1.0;
    bool plummer = true, radial = false, oscil = false, testF = true;
    //@}
    

    //! EDMD modes
    Eigen::MatrixXcd Xi;

    //! DMD state matrices
    Eigen::MatrixXd G, A, X;

    //! Eigenvalue tolerance
    double tol = 1.0e-6;

    //! Eigenvalues values
    Eigen::VectorXd S2, SP;

    //! Koopman eigenvalues
    Eigen::VectorXcd SR, SL;

    //! Koopman operator estimate
    Eigen::MatrixXd K;

    //! Eigenvectors
    Eigen::MatrixXd Q;

    //! Left and right eigenvectors
    Eigen::MatrixXcd U, V;

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

    //! RKHS parameters
    double d, alpha;

    //! RKHS type
    RKHS rkhs = RKHS::Polynomial;

    //! RKHS kernel
    double kernel(const Eigen::VectorXd& x,
		  const Eigen::VectorXd& y);

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
    KoopmanRKHS(const mssaConfig& spec, double tol, const std::string flags="");

    //! Destructor
    virtual ~KoopmanRKHS() {}

    //! Get the eigenvalues
    Eigen::VectorXcd eigenvalues()
    {
      if (not computed) analysis();
      return SR;
    }

    //! Return the Koopman modes, the coefficients to the eigenfunctions
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
