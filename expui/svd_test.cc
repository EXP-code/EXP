#include <iostream>
#include <Eigen/Dense>

namespace MSSA {
  void SvdSignChoice
  (const Eigen::MatrixXd& X,
   Eigen::MatrixXd& U, const Eigen::VectorXd& S, Eigen::MatrixXd& V,
   bool sample);
}

int main()
{
  srand((unsigned int) time(0));

  // Instantiate a matrix filled with random numbers
  const int NROWS = 32;    // N
  const int MCOLS = 32;    // M

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(NROWS, MCOLS);
  std::cout << "Here is the data matrix A:" << std::endl << A
	    << std::endl << std::endl;

  // Carry out the SVD decomposition of the matrix, A such that,
  // A = U S V^T where A is NxM, and 
  //             U is NxM, S=diag(w_{0}, w_{1} ...) is MxM, and V is MxM.
  // U is column-orthogonal
  // D is the diagonal matrix with the singular values
  // V is orthogonal
  
  Eigen::JacobiSVD<Eigen::MatrixXd>
    svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

  std::cout << "Its singular values are:" << std::endl 
	    << svd.singularValues() << std::endl << std::endl;

  std::cout << "Its left singular vectors are the columns of the thin U matrix:" 
	    << std::endl << svd.matrixU() << std::endl << std::endl;

  std::cout << "Its right singular vectors are the columns of the thin V matrix:" 
	    << std::endl << svd.matrixV() << std::endl << std::endl;

  // Test that this SVD is doing what we expect.
  //
  Eigen::MatrixXd S = svd.singularValues().asDiagonal();

  std::cout << std::endl << "Singular value matrix:" << std::endl;
  std::cout << S << std::endl << std::endl;
  
  Eigen::MatrixXd B (NROWS, MCOLS); // N x M
  Eigen::MatrixXd V (MCOLS, MCOLS); // M x M
  Eigen::MatrixXd VT(MCOLS, MCOLS); // Transpose of V
  Eigen::MatrixXd U (NROWS, MCOLS); // N x M
  
  V=svd.matrixV();
  U=svd.matrixU();
  VT=(svd.matrixV()).transpose();
  
  B = U*S*VT; // Recompose the matrix from the decomposition
  
  std::cout << "Reconstructed data matrix" << std::endl;
  std::cout << B << std::endl << std::endl;
  
  Eigen::MatrixXd C(NROWS,MCOLS);
  C = A - B;
  std::cout << "Difference with input data" << std::endl;
  std::cout << C << std::endl;
  std::cout << "Norm: " << C.norm() << std::endl << std::endl;
  
  Eigen::MatrixXd D(MCOLS, MCOLS);
  D = ((svd.matrixU()).transpose())*U;
  std::cout << "Left vector orthogonality" << std::endl;
  std::cout << D << std::endl << std::endl;

  D = ( (svd.matrixV()).transpose() ) * svd.matrixV();
  std::cout << "Right vector orthogonality" << std::endl;
  std::cout << D << std::endl << std::endl;

  // Sign correction test
  {
    auto U = svd.matrixU();
    auto S = svd.singularValues();
    auto V = svd.matrixV();

    MSSA::SvdSignChoice(A, U, S, V, true);

    std::cout << "Corrected left singular vectors:"
	      << std::endl << U << std::endl << std::endl;
    
    std::cout << "Corrected right singular vectors:"
	      << std::endl << V << std::endl << std::endl;

    Eigen::MatrixXd C = A - U * S.asDiagonal() * V.transpose();
    std::cout << "Difference with input data" << std::endl;
    std::cout << C << std::endl;
    std::cout << "Norm: " << C.norm() << std::endl;

    if (C.norm() > 1.0e-10) {

      for (int k=0; k<svd.singularValues().size(); k++) {
	std::cout << "Singular value " << k << " = " << svd.singularValues()(k)
		  << std::endl << std::endl;
	  
	std::cout << svd.matrixU().col(k) * svd.matrixV().col(k).transpose()
		  << std::endl << std::endl;

	std::cout << U.col(k) * V.col(k).transpose()
		  << std::endl << std::endl;
      }
    }
  }

}
