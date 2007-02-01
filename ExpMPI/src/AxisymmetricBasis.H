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
  //! Number of degrees of freedom (e.g. 3 for sphere)
  int dof;
  //! Use principal component analysis to select coefficients
  bool pca;
  //! Frequency of PCA recomputation
  int npca;

  //! Expansion coefficients (summed over all nodes)
  Matrix expcoef;
  //! Expansion coefficients (per node)
  Matrix expcoef1;
  //! Normalization for covariance matrix based on the biorthogonal basis norm
  Matrix normM;
  //! The variance for all nodes and per node (respectively)
  Matrix *cc, *cc1;

  /** MPI routines to distribution the expansion coefficients
      (<code>expcoef</code>) from all nodes. */
  void parallel_distribute_coefficients(void);
  /** MPI routines to sum the expansion coefficients
      (<code>expcoef</code>) from all nodes. */
  void parallel_gather_coefficients(void);

public:

  /*! Enum listing the possible selection algorithms for coefficient
    selection */
  enum TKType {
    Hall,             /*!< Tapered signal-to-noise power defined by Hall */
    VarianceCut,      /*!< Cut all coefficients below some S/N level */
    CumulativeCut,    /*!< Cut coefficients below some cumulative fraction */
    VarianceWeighted  /*!< Weight coefficients be S/N for S/N<1 */
  };

  //! Constructor
  AxisymmetricBasis(string& line);
  //! Destructor
  virtual ~AxisymmetricBasis();

  //! Compute component noise suppression using PCA analysis
  void pca_hall(int compute);

private:

  int Ldim, L0;
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
  Matrix sqnorm;

  Matrix covar;

};


#endif // _AxisymmetricBasis_H


