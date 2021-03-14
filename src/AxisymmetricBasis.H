#ifndef _AxisymmetricBasis_H
#define _AxisymmetricBasis_H

#include <Basis.H>
#include <Vector.h>

//! Defines a basis-based potential and acceleration class
/*! Configurable parameters are:

  @param Lmax  is the maximum spherical harmonic order

  @param nmax  is the maximum radial order

  @param dof   is the dimensionality of the basis

  @param npca  is the frequency for pca analyses

  @param npca0 is the number of steps to skip before the first pca analysis

  @param pcavar boolean turns on pca analyses

  @param pcaeof boolean turns on emperical function recomputation

  @param pcadiag boolean turns on diagnostic pca output

  @param pcavtk boolean prints out visualizable representations of pca info

  @param samplesz is the number of subsamples for pca variance analysis (default: 0 means use sqrt number o particles

  @param vtkfreq is the frequency for printing pca variance analysis info

  @param tksmooth is the S/N for Hall truncation for VarianceCut method

  @param tkcum is the cumulative fraction eigenvalue power for truncation for CmululativeCut method

  @param tk_type is the smoothing type, one of: Hall, VarianceCut, CumulativeCut, VarianceWeighted
*/
class AxisymmetricBasis : public Basis
{

protected:

  //! Maximum angular harmonic
  int Lmax;

  //! Maximum radial basis index
  int nmax;

  //! Use principal component boostrap analysis of variance to select
  //! coefficients
  bool pcavar;

  //! Output PCA diagnostics
  bool pcadiag;

  //! Output PCA diagnostics in vtk format
  bool pcavtk;

  //! Compute eof decomposition
  bool pcaeof;

  //! VTK diagnostic frequency
  unsigned vtkfreq;

  //! Frequency of PCA computation
  int npca;

  //! First step for PCA computation
  int npca0;

  /**
     Expansion coefficients (summed over all nodes)

     The vector entries store coefficients for a multistep
     expansion.  Move pointers rather than objects for
     efficiency.
  */
  //@{
  typedef boost::shared_ptr<Vector> VectorP;
  typedef boost::shared_ptr<Matrix> MatrixP;

  std::vector<VectorP> expcoef, expcoef1;

  //@{
  //! Interpolation arrays
  std::vector<std::vector<VectorP>> expcoefN;
  std::vector<std::vector<VectorP>> expcoefL;
  //@}

  //@{
  //! Covariance arrays
  std::vector<std::vector<VectorP>> expcoefT, expcoefT1;
  std::vector<std::vector<MatrixP>> expcoefM, expcoefM1;
  //@}

  //@{
  //! Mass and counts for subsample covariance
  std::vector<double>               massT, massT1;
  unsigned                          sampT, defSampT;
  //@}

  //! Test for eof trim
  std::vector<MatrixP>              tvar;
  //@}

  //! Normalization for covariance matrix based on the biorthogonal basis norm
  Matrix normM;

  //! Sqr root of the normalization
  Matrix sqnorm;

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
  AxisymmetricBasis(const YAML::Node& line);

  //! Destructor
  virtual ~AxisymmetricBasis();

  //! Compute component noise suppression using PCA analysis
  void pca_hall(bool compute);

  //! Reset used particle counter
  virtual void multistep_reset() { used=0; }

  //! Set tk_type from string
  TKType setTK(const std::string& tk);

private:

  int    Ldim, L0;
  double fac02, tksmooth, tkcum;
  TKType tk_type;

  Vector smth;
  Vector *weight;
  Vector *b_Hall;
  Vector inv;
  Vector eval;
  Vector cuml;
  Matrix *evec;
  Matrix *Tevec;
  Matrix covar;

};


#endif // _AxisymmetricBasis_H


