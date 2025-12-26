#ifndef _AxisymmetricBasis_H
#define _AxisymmetricBasis_H

#include "Basis.H"
#include <Eigen/Eigen>

/** @class AxisymmetricBasis
    @brief Defines a basis-based potential and acceleration class

    @details **YAML configuration**

    @param Lmax  is the maximum spherical harmonic order

    @param Mmax  is the maximum cylindrical harmonic order

    @param nmax  is the maximum radial order

    @param dof   is the dimensionality of the basis

    @param npca  is the frequency for pca analyses

    @param npca0 is the number of steps to skip before the first pca analysis

    @param pcavar boolean turns on pca analyses

    @param pcaeof boolean turns on emperical function recomputation

    @param pcadiag boolean turns on diagnostic pca output

    @param pcavtk boolean prints out visualizable representations of pca info

    @param hexp is the Hall smoothing exponent

    @param samplesz is the number of subsamples for pca variance
    analysis (default: 1). Value 0 means use sqrt number of particles

    @param vtkfreq is the frequency for printing pca variance analysis info

    @param tksmooth is the S/N for Hall truncation for VarianceCut method

    @param tkcum is the cumulative fraction eigenvalue power for truncation for CmululativeCut method

    @param tk_type is the smoothing type, one of: Hall, VarianceCut, CumulativeCut, VarianceWeighted
    
    @param subsamp true sets partition variance computation (default: false)

    @param nint is the step cadence for subsample variance computation
    (default: 0, no subsampling)
*/
class AxisymmetricBasis : public Basis
{

protected:

  //! Maximum angular harmonic
  int Lmax, Mmax;

  //! Maximum radial basis index
  int nmax;

  //! Use principal component bootstrap analysis of variance to select
  //! coefficients
  bool pcavar;

  //! Output PCA diagnostics
  bool pcadiag;

  //! Output PCA diagnostics in vtk format
  bool pcavtk;

  //! Compute eof decomposition
  bool pcaeof;

  //! VTK diagnostic frequency
  unsigned vtkfreq = std::numeric_limits<unsigned>::max();

  //! Frequency of PCA computation
  int npca = std::numeric_limits<int>::max();

  //! First step for PCA computation
  int npca0 = std::numeric_limits<int>::max();

  // Frequency of subsample computation
  int nint = 0;

  //! Hall smoothing exponent (default: 1.0)
  double hexp;

  //! Signal-to-noise scaling parameter (default: 1.0)
  double snr;

  /**
     Expansion coefficients (summed over all nodes)

     The vector entries store coefficients for a multistep
     expansion.  Move pointers rather than objects for
     efficiency.
  */
  //@{
  typedef std::shared_ptr<Eigen::VectorXd > VectorP;
  typedef std::shared_ptr<Eigen::MatrixXd > MatrixP;
  typedef std::shared_ptr<Eigen::VectorXcd> VectorQ;
  typedef std::shared_ptr<Eigen::MatrixXcd> MatrixQ;

  std::vector<VectorP> expcoef, expcoef1, expcoefP;

  //! Swap coefficients
  void swap_coefs(std::vector<VectorP>& from, std::vector<VectorP>& to)
  {
    if (from.size() != to.size()) {
      std::ostringstream sout;
      sout << "swap_coefs: size(from)=" << from.size() << " != "
	   << "size(to)=" << to.size();
      throw std::runtime_error(sout.str());
    }

    std::vector<VectorP> tmp(from);
    from = to;
    to   = tmp;
  }

  //@{
  //! Interpolation arrays
  std::vector<std::vector<VectorP>> expcoefN;
  std::vector<std::vector<VectorP>> expcoefL;
  //@}

  //@{
  //! Covariance arrays
  std::vector<std::vector<VectorQ>> expcoefT, expcoefT1;
  std::vector<std::vector<MatrixQ>> expcoefM, expcoefM1;
  //@}

  //@{
  //! Mass and counts for subsample covariance
  std::vector<unsigned>              countT, countT1;
  std::vector<double>                massT, massT1;
  unsigned                           sampT, defSampT;
  //@}

  //! Test for eof trim
  std::vector<MatrixP>              tvar;
  //@}

  //! Normalization for covariance matrix based on the biorthogonal basis norm
  Eigen::MatrixXd normM;

  //! Sqr root of the normalization
  Eigen::MatrixXd sqnorm;

  //! Accumulated mass for variance computation
  double muse;

  /** MPI routines to distribution the expansion coefficients
      (<code>expcoef</code>) from all nodes. */
  virtual void parallel_distribute_coefficients(void);


  /** MPI routines to sum the expansion coefficients
      (<code>expcoef</code>) from all nodes. */
  virtual void parallel_gather_coefficients(void);

  /** MPI routines to sum the expansion coefficients variance moments
      (<code>cc</code>) from all nodes. */
  virtual void parallel_gather_coef2(void);

  //! Partitioned variance computation
  bool subsamp;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  /*! Enum listing the possible selection algorithms for coefficient
    selection */
  enum TKType {
    Hall,             /*!< Tapered signal-to-noise power defined by Hall   */
    VarianceCut,      /*!< Cut all coefficients below some S/N level       */
    CumulativeCut,    /*!< Cut coefficients below some cumulative fraction */
    VarianceWeighted, /*!< Weight coefficients be S/N for S/N<1            */
    None              /*!< Compute the S/N but do not modify coefficients  */
  };

  //! Constructor
  AxisymmetricBasis(Component* c0, const YAML::Node& line);

  //! Destructor
  virtual ~AxisymmetricBasis();

  //@{
  //! Compute component noise suppression using PCA analysis
  void pca_hall(bool compute);	  // Choose between 2d and 3d
  void pca_hall_2d(bool compute); // 2d-polar implementation
  void pca_hall_3d(bool compute); // 3d spherical harmonic implementation
  //@}

  //! Reset used particle counter
  virtual void multistep_reset() { used=0; }

  //! Set tk_type from string
  TKType setTK(const std::string& tk);

  //! Access coefficients
  VectorP getCoefs(int l) const { return expcoef[l]; }

  //@{
  //! Get angular order
  int getLmax() { return Lmax; }
  int getMmax() { return Mmax; }
  //@{

  //! Get radial order
  int getNmax() { return nmax; }

private:

  int    Ldim, L0;
  double fac02, tksmooth, tkcum;
  TKType tk_type;

  Eigen::VectorXd smth;
  std::vector<Eigen::VectorXd> weight;
  std::vector<Eigen::VectorXd> s_Hall;
  std::vector<Eigen::VectorXd> b_Hall;
  Eigen::VectorXd inv;
  Eigen::VectorXd eval;
  Eigen::VectorXd cuml;
  std::vector<Eigen::MatrixXd> evec;
  std::vector<Eigen::MatrixXd> Tevec;
  Eigen::MatrixXd covar;

};


#endif // _AxisymmetricBasis_H


