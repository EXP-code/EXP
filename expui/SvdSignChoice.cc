#include <Eigen/Dense>

namespace MSSA {

  // SVD with corrected signs
  //
  // Input
  // -----
  // X is data matrix (constant)
  // U is the matrix of left singular vectors
  // S is the array of singular values (constant)
  // V is the matrix of right singular vectors
  //
  // Output
  // ------
  // U, V returned corrected, disambiguated signs
  //
  // Reference
  // ---------
  // Bro, R., Acar, E., & Kolda, T. G. (2008). Resolving the sign
  // ambiguity in the singular value decomposition.  Journal of
  // Chemometrics: A Journal of the Chemometrics Society, 22(2),
  // 135-140.
  //
  // URL:
  // https://prod-ng.sandia.gov/techlib-noauth/access-control.cgi/2007/076422.pdf
  //
  void SvdSignChoice
  (const Eigen::MatrixXd& X,
   Eigen::MatrixXd& U, const Eigen::VectorXd& S, Eigen::MatrixXd& V)
  {
    // SDV dimensions
    int I = U.rows();
    int J = V.rows();
    int K = S.size();
    
    // Dimensions
    // ----------
    // X is a [I x J] data matrix
    // U is a [I x K] matrix
    // S is a [K x 1] vector (diagonal matrix)
    // V is a [J x K] matrix
    
    // Sanity checks
    if (U.cols() != K)
      throw std::invalid_argument("SvdSignChoice: U has wrong dimensions");
    
    if (V.cols() != K)
      throw std::invalid_argument("SvdSignChoice: V has wrong dimensions");
    
    if (X.rows() != I || X.cols() != J)
      throw std::invalid_argument("SvdSignChoice: X dimensions do not match SVD input");
    
    // Sign determination loop
    //
    Eigen::VectorXd sL(K), sR(K);
    sL.setZero();
    sR.setZero();
    
    // Working, non-const instance
    auto S1 = S;
    
    // Get projections from left and right singular vectors onto data
    // matrix
    //
    for (int k=0; k<K; k++) {
      // Remove all but target dimension for numerical stability
      S1(k) = 0.0;
      auto Y = X - U * S1.asDiagonal() * V.transpose();

      // Restore the value
      S1(k) = S(k);

      // d_j = U^T_k * Y_j
      Eigen::VectorXd dL = Y.transpose() * U.col(k);

      // sum of sgn(dL_j)*dL_j^2
      sL(k) += dL.dot(dL.cwiseAbs());

      // d_i = V^T_k * (Y^T)_i
      Eigen::VectorXd dR = Y * V.col(k);

      // sum of sgn(dR_i)*dR_i^2
      sR(k) += dR.dot(dR.cwiseAbs());
    }
    
    auto sgn = [](double val) -> int
    {
      return (0.0 < val) - (val < 0.0);
    };
    
    // Determine and apply the sign correction
    //
    for (int k=0; k<K; k++) {
      // If signs are opposite, flip the one with the smaller absolute
      // value
      //
      if (sL(k)*sR(k) < 0.0) {
	if (std::abs(sL(k)) < std::abs(sR(k)))
	  sL(k) = -sL(k);
	else
	  sR(k) = -sR(k);
      }
      
      // Apply
      U.col(k) *= sgn(sL(k));
      V.col(k) *= sgn(sR(k));
    }
    
    // Done
  }

  void SvdSignChoice
  (const Eigen::MatrixXd& X,
   Eigen::MatrixXd& U, const Eigen::VectorXd& S, const Eigen::MatrixXd& V)
  {
    // SDV dimensions
    int I = U.rows();
    int J = V.rows();
    int K = S.size();
    
    // Dimensions
    // ----------
    // X is a [I x J] data matrix
    // U is a [I x K] matrix
    // S is a [K x 1] vector (diagonal matrix)
    // V is a [J x K] matrix
    
    // Sanity checks
    if (U.cols() != K)
      throw std::invalid_argument("SvdSignChoice: U has wrong dimensions");
    
    if (V.cols() != K)
      throw std::invalid_argument("SvdSignChoice: V has wrong dimensions");
    
    if (X.rows() != I || X.cols() != J)
      throw std::invalid_argument("SvdSignChoice: X dimensions do not match SVD input");
    
    // Sign determination loop
    //
    Eigen::VectorXd sL(K), sR(K);
    sL.setZero();
    sR.setZero();
    
    // Working, non-const instance
    auto S1 = S;
    
    // Get projections from left and right singular vectors onto data
    // matrix
    //
    for (int k=0; k<K; k++) {
      // Remove all but target dimension for numerical stability
      S1(k) = 0.0;
      auto Y = X - U * S1.asDiagonal() * V.transpose();

      // Restore the value
      S1(k) = S(k);

      // d_j = U^T_k * Y_j
      Eigen::VectorXd dL = Y.transpose() * U.col(k);

      // sum of sgn(dL_j)*dL_j^2
      sL(k) += dL.dot(dL.cwiseAbs());

      // d_i = V^T_k * (Y^T)_i
      Eigen::VectorXd dR = Y * V.col(k);

      // sum of sgn(dR_i)*dR_i^2
      sR(k) += dR.dot(dR.cwiseAbs());
    }
    
    auto sgn = [](double val) -> int
    {
      return (0.0 < val) - (val < 0.0);
    };
    
    // Determine and apply the sign correction
    //
    for (int k=0; k<K; k++) {
      // If signs are opposite, flip the one with the smaller absolute
      // value
      //
      if (sL(k)*sR(k) < 0.0) {
	if (std::abs(sL(k)) < std::abs(sR(k)))
	  sL(k) = -sL(k);
	else
	  sR(k) = -sR(k);
      }
      
      // Apply
      U.col(k) *= sgn(sL(k));
    }
    
    // Done
  }

  void SvdSignChoice
  (const Eigen::MatrixXd& X,
   const Eigen::MatrixXd& U, const Eigen::VectorXd& S, Eigen::MatrixXd& V)
  {
    // SDV dimensions
    int I = U.rows();
    int J = V.rows();
    int K = S.size();
    
    // Dimensions
    // ----------
    // X is a [I x J] data matrix
    // U is a [I x K] matrix
    // S is a [K x 1] vector (diagonal matrix)
    // V is a [J x K] matrix
    
    // Sanity checks
    if (U.cols() != K)
      throw std::invalid_argument("SvdSignChoice: U has wrong dimensions");
    
    if (V.cols() != K)
      throw std::invalid_argument("SvdSignChoice: V has wrong dimensions");
    
    if (X.rows() != I || X.cols() != J)
      throw std::invalid_argument("SvdSignChoice: X dimensions do not match SVD input");
    
    // Sign determination loop
    //
    Eigen::VectorXd sL(K), sR(K);
    sL.setZero();
    sR.setZero();
    
    // Working, non-const instance
    auto S1 = S;
    
    // Get projections from left and right singular vectors onto data
    // matrix
    //
    for (int k=0; k<K; k++) {
      // Remove all but target dimension for numerical stability
      S1(k) = 0.0;
      auto Y = X - U * S1.asDiagonal() * V.transpose();

      // Restore the value
      S1(k) = S(k);

      // d_j = U^T_k * Y_j
      Eigen::VectorXd dL = Y.transpose() * U.col(k);

      // sum of sgn(dL_j)*dL_j^2
      sL(k) += dL.dot(dL.cwiseAbs());

      // d_i = V^T_k * (Y^T)_i
      Eigen::VectorXd dR = Y * V.col(k);

      // sum of sgn(dR_i)*dR_i^2
      sR(k) += dR.dot(dR.cwiseAbs());
    }
    
    auto sgn = [](double val) -> int
    {
      return (0.0 < val) - (val < 0.0);
    };
    
    // Determine and apply the sign correction
    //
    for (int k=0; k<K; k++) {
      // If signs are opposite, flip the one with the smaller absolute
      // value
      //
      if (sL(k)*sR(k) < 0.0) {
	if (std::abs(sL(k)) < std::abs(sR(k)))
	  sL(k) = -sL(k);
	else
	  sR(k) = -sR(k);
      }
      
      // Apply
      V.col(k) *= sgn(sR(k));
    }
    
    // Done
  }

}
