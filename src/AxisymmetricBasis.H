#ifndef _AxisymmetricBasis_H
#define _AxisymmetricBasis_H

#include <Basis.H>
#include <Vector.h>

//! Defines a basis-based potential and acceleration class
class AxisymmetricBasis : public Basis
{

protected:

  //! Maximum angular harmonic
  int Lmax;

  //! Maximum radial basis index
  int nmax;

  //! Use principal component analysis to select coefficients
  bool pca;

  //! Output PCA diagnostics
  bool pcadiag;

  //! Output PCA diagnostics in vtk format
  bool pcavtk;

  //! vtk diagnostic frequency
  unsigned vtkfreq;

  //! PCA jackknife variance
  bool pcajknf;

  //! Frequency of PCA recomputation
  int npca;

  /**
     Expansion coefficients (summed over all nodes)

     The vector entries store coefficients for a multistep
     expansion.  Move pointers rather than objects for
     efficiency.
  */
  //@{
  Matrix expcoef, expcoef1;
  vector<Matrix*> expcoefN;
  vector<Matrix*> expcoefL;

  typedef boost::shared_ptr<Matrix> MatrixP;
  vector<MatrixP> expcoefT, expcoefT1;
  vector<double>  massT, massT1;
  unsigned        sampT;
  //@}

  //! Normalization for covariance matrix based on the biorthogonal basis norm
  Matrix normM;

  //! Sqr root of the normalization
  Matrix sqnorm;

  //@{
  //! The variance for all nodes and per node (respectively)
  Matrix *cc, *cc1;
  //@}

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
    Null              /*!< Compute the S/N but do not modify coefficients  */
  };

  //! Constructor
  AxisymmetricBasis(const YAML::Node& line);

  //! Destructor
  virtual ~AxisymmetricBasis();

  //! Compute component noise suppression using PCA analysis
  void pca_hall(int compute);

  //! Reset used particle counter
  virtual void multistep_reset() { used=0; }

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
  Matrix Tevec;
  Matrix covar;

};


#endif // _AxisymmetricBasis_H


