#define TESTING
//
// EDMD (Liouville occupuation kernel version) for EXP coefficients
//
// Uses fixed-rank approximation for the SVD to save time and space.
// Uses the approximate SVD randomized algorithm from Halko,
// Martinsson, and Tropp by default.  Use BDCSVD or Jacobi flags for
// the standard methods.  Using the randomized method together with
// trajectory matrix rather than covariance matrix analysis may lead
// to large errors in eigenvectors.
//
// The implementation here is based on the
//
// Joel A. Rosenfeld and Rushikesh Kamalapurkar, "Singular Dynamic
// Mode Decomposition", 2023, SIAM J. Appl. Dynamical Systems,
// Vol. 22, No. 3, pp 2357-2381
//

#include <filesystem>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cmath>
#include <map>

#include <config_exp.h>

#include <Eigen/Dense>

#include <highfive/highfive.hpp>
#include <highfive/eigen.hpp>

#include <omp.h>

/* For debugging
   #undef eigen_assert
   #define eigen_assert(x)						\
   if (!(x)) { throw (std::runtime_error("Eigen assert error")); }
*/

#include <omp.h>

#include <LiouvilleRKHS.H>

#include <RedSVD.H>
#include <YamlConfig.H>
#include <YamlCheck.H>
#include <EXPException.H>
#include <libvars.H>

namespace MSSA {

  // Helper ostream manipulator for debug coefficient key info
  //
  std::ostream& operator<< (std::ostream& out, const std::vector<unsigned>& t);

  // Get RKHS name from enum
  //
  std::map<LiouvilleRKHS::RKHS, std::string> LiouvilleRKHS::RKHS_names =
    {
      {LiouvilleRKHS::RKHS::Polynomial,  "Polynomial" },
      {LiouvilleRKHS::RKHS::Exponential, "Exponential"},
      {LiouvilleRKHS::RKHS::Gaussian,    "Gaussian"   }
    };

  // Get RKHS enum from name
  //
  std::map<std::string, LiouvilleRKHS::RKHS> LiouvilleRKHS::RKHS_values =
    {
      {"Polynomial",  LiouvilleRKHS::RKHS::Polynomial },
      {"Exponential", LiouvilleRKHS::RKHS::Exponential},
      {"Gaussian",    LiouvilleRKHS::RKHS::Gaussian   }
    };


  // Get Method name from enum
  //
  std::map<LiouvilleRKHS::Method, std::string> LiouvilleRKHS::Method_names =
    {
      {LiouvilleRKHS::Method::Singular,      "Singular" },
      {LiouvilleRKHS::Method::Eigenfunction, "Eigenfunction"}
    };

  // Get RKHS enum from name
  //
  std::map<std::string, LiouvilleRKHS::Method> LiouvilleRKHS::Method_values =
    {
      {"Singular",      LiouvilleRKHS::Method::Singular   },
      {"Eigenfunction", LiouvilleRKHS::Method::Eigenfunction}
    };


  // RKHS kernel
  //
  double LiouvilleRKHS::kernel(const Eigen::VectorXd& x,
			       const Eigen::VectorXd& y,
			       double mu)
  {
    double D2 = d*d/(mu*mu);

    if (rkhs == RKHS::Polynomial) {
      double prod = x.adjoint()*y;
      return pow(1.0 + prod/D2, alpha);
    } else if (rkhs == RKHS::Exponential) {
      double prod = x.adjoint()*y;
      return exp(prod/D2);
    } else if (rkhs == RKHS::Gaussian) {
      Eigen::VectorXd diff = x - y;
      double prod = diff.adjoint()*diff;
      return exp(-prod/D2);
    } else {
      throw std::runtime_error("LiouvilleRKHS: Unknown kernel type");
    }
  }

  // Structure for sorting complex numbers and retrieving a permutation index
  //
  using complex_elem = std::pair<std::complex<double>, int>;
  static
  bool complex_comparator(const complex_elem &lhs, const complex_elem &rhs)
  {
    auto a = lhs.first;
    auto b = rhs.first;

    double abs_a = std::abs(a);
    double abs_b = std::abs(b);

    // Three step sort:
    // (1) By modulus
    // (2) By real value
    // (3) By imag value
    //
    if (abs_a == abs_b) {
      if (std::real(a) == std::real(b))
	return std::imag(a) < std::imag(b);
      else
	return std::real(a) < std::real(b);
    } else {
      return abs_a < abs_b;
    }
  }

  // Fixed-width matrix element printing
  //
  static
  void matrix_print(std::ostream& out, const Eigen::MatrixXcd& M)
  {
    for (int i=0; i<M.rows(); i++) {
      for (int j=0; j<M.cols(); j++) {
	out << "(" << std::setw(8) << std::fixed << std::setprecision(3) << M(i, j).real()
	    << "," << std::setw(8) << std::fixed << std::setprecision(3) << M(i, j).imag() << ") ";
      }
      out << std::endl;
    }
  }

  // Check for similarity to the identity matrix
  //
  static
  bool matrix_idtest(std::ostream& out, const Eigen::MatrixXcd& M)
  {
    double max_diag = 0.0, max_offd = 0.0;
    for (int i=0; i<M.rows(); i++) {
      for (int j=0; j<M.cols(); j++) {
	if (i==j) max_diag = std::max(max_diag, std::abs(M(i, j)-1.0));
	else      max_offd = std::max(max_offd, std::abs(M(i, j)    ));
      }
    }
    out << std::string(70, '-') << std::endl << std::setprecision(8);
    out << "diag diff=" << max_diag << " offd diff=" << max_offd << std::endl;
    out << std::string(70, '-') << std::endl;

    if (max_diag>1.0e-3 or max_offd>1.0e-3) return true;
    return false;
  }

  // Compute G matrix for RKHS space defined by mu value
  Eigen::MatrixXd LiouvilleRKHS::computeGrammian(double mu, double rat)
  {
    // Use Simpson's rule quadrature
    
    // Number of trajectories
    //
    auto keys = getAllKeys();

    if (keys.size() != traj)
      throw std::runtime_error("LiouvilleRKHS::computeGrammian: "
			       "number of trajectories does not match");
    
    // Check for odd number of times for Simpson's 1/3 rule
    //
    int nt = numT, f1, f2;
    if (nt % 2 == 0) nt -= 1;

    // Allocate Grammian matrix and set to zero
    //
    Eigen::MatrixXd G(traj, traj);

#pragma omp parallel for collapse(2)
    for (int i=0; i<traj; i++) {

      for (int j=0; j<traj; j++) {

	G(i, j) = 0.0;

	for (int t1=0; t1<nt; t1++) {

	  int f1 = 4;			   // Even interior point
	  if (t1 == 0 or t1 == nt-1) f1 = 1; // End points
	  else if (t1 % 2) f1 = 2;	   // Odd interior point
	
	  double inner = 0.0;
	  for (int t2=0; t2<nt; t2++) {
	  
	    int f2 = 4;
	    if (t2 == 0 or t2 == nt-1) f2 = 1;
	    else if (t2 % 2) f2 = 2;
	  
	    auto x1 = data[keys[i]].row(t1)*rat;
	    auto x2 = data[keys[j]].row(t2)*rat;
	    inner += kernel(x1, x2, mu) * f2;
	  }
	
	G(i, j) += inner * f1;
	}
      }
    }

    /*
    for (int i=0; i<traj; i++) {

      for (int j=0; j<traj; j++) {

	G(i, j) = 0.0;

	for (int t1=0; t1<nt; t1++) {

	  if (t1 == 0 or t1 == nt-1) f1 = 1; // End points
	  else if (t1 % 2) f1 = 2;	     // Odd interior point
	  else f1 = 4;			     // Even interior point
	  
	  double inner = 0.0;
	  for (int t2=0; t2<nt; t2++) {
	    
	    if (t2 == 0 or t2 == nt-1) f2 = 1;
	    else if (t2 % 2) f2 = 2;
	    else f2 = 4;
	    
	    auto x1 = data[keys[i]].row(t1)*rat;
	    auto x2 = data[keys[j]].row(t2)*rat;
	    inner += kernel(x1, x2, mu) * f2;
	  }
	  
	  G(i, j) += inner * f1;
	}
      }
    }
    */
      
    // Time step factor
    //
    double SF = (coefDB.times[1] - coefDB.times[0])/3.0;

    return G*SF*SF;
  }

  // Compute G matrix for RKHS space defined by mu value
  Eigen::MatrixXd LiouvilleRKHS::computeGammaDiff(double mu)
  {
    // Use Simpson's rule quadrature
    
    // Number of trajectories
    //
    auto keys = getAllKeys();
    if (traj != keys.size())
      throw std::runtime_error("LiouvilleRKHS::computeGammaDiff: "
			       "number of trajectories does not match");
    
    // Check for odd number of times for Simpson's 1/3 rule
    //
    int nt = numT, f;
    if (nt % 2 == 0) nt -= 1;

    // Allocate Grammian matrix and set to zero
    //
    Eigen::MatrixXd A(traj, traj);

#pragma omp parallel for collapse(2)
    for (int i=0; i<traj; i++) {

      for (int j=0; j<traj; j++) {


	A(i, j) = 0.0;

	double sumT = 0.0, sum0 = 0.0;
	
	for (int t=0; t<nt; t++) {
	
	  int  f = 4;			  // Even interior point
	  if (t == 0 or t == nt-1) f = 1; // End points
	  else if (t % 2) f = 2;	  // Odd interior point
	
	  sum0 += kernel(data[keys[i]].row(0   ), data[keys[j]].row(t), mu) * f;
	  sumT += kernel(data[keys[i]].row(nt-1), data[keys[j]].row(t), mu) * f;
	}

	A(i, j) += sumT - sum0;
      }
    }

    /*
    for (int i=0; i<traj; i++) {

      for (int j=0; j<traj; j++) {
	      
	A(i, j) = 0.0;

	double sumT = 0.0, sum0 = 0.0;
	
	for (int t=0; t<nt; t++) {
	  
	  if (t == 0 or t == nt-1) f = 1; // End points
	  else if (t % 2) f = 2;		// Odd interior point
	  else f = 4;			// Even interior point
	  
	  sum0 += kernel(data[keys[i]].row(0   ), data[keys[j]].row(t), mu) * f;
	  sumT += kernel(data[keys[i]].row(nt-1), data[keys[j]].row(t), mu) * f;
	}
	
	A(i, j) += sumT - sum0;
      }
    }
    */

    // Time step factor
    //
    double SF = (coefDB.times[1] - coefDB.times[0])/3.0;

    return A * SF;
  }

  // Compute G matrix for RKHS space defined by mu value
  Eigen::MatrixXd LiouvilleRKHS::computeGramGamma(double mu)
  {
    // Number of trajectories
    //
    auto keys = getAllKeys();
    if (traj != keys.size())
      throw std::runtime_error("LiouvilleRKHS::computeGammaDiff: "
			       "number of trajectories does not match");
    
    // Allocate Grammian matrix and set to zero
    //
    Eigen::MatrixXd A(traj, traj);

    int lT = numT - 1;

#pragma omp parallel for
    for (int k=0; k<traj*(traj+1)/2; k++) {
    
      // This computes indices for i=[0, traj) j=[i, traj)
      //
      int i = floor( (sqrt(1.0 + 8.0*k) - 1.0)/2.0 );
      int j = k - i*(i+1)/2;

      A(i, j) =
	kernel(data[keys[i]].row(lT), data[keys[j]].row(lT), mu) +
	kernel(data[keys[i]].row(0 ), data[keys[j]].row(0 ), mu) -
	kernel(data[keys[i]].row(0 ), data[keys[j]].row(lT), mu) -
	kernel(data[keys[i]].row(lT), data[keys[j]].row(0 ), mu) ;
      
      if (i != j) A(j, i) = A(i, j); // Matrix is symmetric
    }

    /*
    for (int i=0; i<traj; i++) {
      for (int j=i; j<traj; j++) {
	A(i, j) =
	  kernel(data[keys[i]].row(lT), data[keys[j]].row(lT), mu) +
	  kernel(data[keys[i]].row(0 ), data[keys[j]].row(0 ), mu) -
	  kernel(data[keys[i]].row(0 ), data[keys[j]].row(lT), mu) -
	  kernel(data[keys[i]].row(lT), data[keys[j]].row(0 ), mu) ;

	if (i != j) A(j, i) = A(i, j); // Matrix is symmetric
      }
    }
    */

    return A;
  }


  // Compute gamma vector for RKHS space defined by mu value
  Eigen::VectorXd LiouvilleRKHS::computeGamma
  (const Eigen::VectorXd& x, double mu)
  {
    // Use Simpson's rule quadrature
    
    auto keys = getAllKeys();
    if (traj != keys.size())
      throw std::runtime_error("LiouvilleRKHS::computeGamma: "
			       "number of trajectories does not match");
    
    // Check for odd number of times for Simpson's 1/3 rule
    //
    int nt = numT;
    if (nt % 2 == 0) nt -= 1;

    // Allocate return vector
    //
    Eigen::VectorXd g(traj);

#pragma omp parallel for
    for (int i=0; i<traj; i++) {

      g(i) = 0.0;

      for (int t=0; t<nt; t++) {
	
	int f = 4;			// Even interior point
	if (t == 0 or t == nt-1) f = 1; // End points
	else if (t % 2) f = 2;		// Odd interior point
	
	g(i) += kernel(x, data[keys[i]].row(t), mu) * f;
      }
    }

    // Time step factor
    //
    double SF = (coefDB.times[1] - coefDB.times[0])/3.0;

    return g * SF;
  }


  // Compute occupation kernel matrix
  Eigen::MatrixXd LiouvilleRKHS::occupation()
  {
    // Use Simpson's rule quadrature
    
    // Number of trajectories
    //
    auto keys = getAllKeys();

    if (traj != keys.size())
      throw std::runtime_error("LiouvilleRKHS::occupation: "
			       "number of trajectories does not match");

    // Check for odd number of times for Simpson's 1/3 rule
    //
    int nt = numT, f;
    if (nt % 2 == 0) nt -= 1;

    // Allocate return vector and set to zero
    //
    Eigen::MatrixXd g(traj, rank);
    g.setZero();
    
    for (int t=0; t<nt; t++) {
	
      if (t == 0 or t == nt-1) f = 1; // End points
      else if (t % 2) f = 2;	      // Odd interior point
      else f = 4;		      // Even interior point
      
#pragma omp parallel for
      for (int i=0; i<traj; i++) {
	g.row(i) += data[keys[i]].row(t) * f;
      }
    }

    // Time step factor
    //
    double SF = (coefDB.times[1] - coefDB.times[0])/3.0;

    return g * SF;
  }

  // Compute trajectory matrix
  Eigen::MatrixXd LiouvilleRKHS::trajectory()
  {
    // Number of trajectories
    //
    auto keys = getAllKeys();

    if (traj != keys.size())
      throw std::runtime_error("LiouvilleRKHS::trajectory: "
			       "number of trajectories does not match");

    Eigen::MatrixXd g(traj, rank);
    g.setZero();
    
    int nT = numT - 1;

    for (int i=0; i<traj; i++)
      g.row(i) = data[keys[i]].row(nT) - data[keys[i]].row(0);

    return g;
  }


  // Algorithm 6.1 from Rosenfeld and Kamalapurkhar
  //
  void LiouvilleRKHS::singular_analysis()
  {
    // Number of channels
    //
    auto keys = getAllKeys();
    nkeys = keys.size();

    // The number of time points
    //
    numT = data[keys[0]].rows();
    rank = data[keys[0]].cols();

    // Get Grammian matrices
    //
    G1 = computeGramGamma(mu1);
    G2 = computeGrammian(mu2);

    // For regularization
    //
    auto R = Eigen::MatrixXd::Identity(traj, traj) * eps;

    // Perform eigenanalysis of Grammian matrix (pos def)
    //
    if (use_red) {
      RedSVD::RedSymEigen<Eigen::MatrixXd> eigensolver1(G1 + R, evCount);
      S1 = eigensolver1.eigenvalues();
      Q1 = eigensolver1.eigenvectors();

      RedSVD::RedSymEigen<Eigen::MatrixXd> eigensolver2(G2 + R, evCount);
      S2 = eigensolver2.eigenvalues();
      Q2 = eigensolver2.eigenvectors();
    } else {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver1(G1 + R);
      if (eigensolver1.info() != Eigen::Success) {
	throw std::runtime_error("LiouvilleRKHS: Eigensolver for Grammian 1 failed");
      }
      S1 = eigensolver1.eigenvalues();
      Q1 = eigensolver1.eigenvectors();

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver2(G2 + R);
      if (eigensolver2.info() != Eigen::Success) {
	throw std::runtime_error("LiouvilleRKHS: Eigensolver for Grammian 2 failed");
      }
      S2 = eigensolver2.eigenvalues();
      Q2 = eigensolver2.eigenvectors();
    }
    
    // Compute the normalized bases
    //
    Eigen::VectorXd D(traj), Dt(traj);

    for (int i=0; i<Q1.cols(); i++) {
      auto dp = Q1.col(i).dot(G1 * Q1.col(i));
      D(i) = 1.0/sqrt(dp);
    }

    for (int i=0; i<Q2.cols(); i++) {
      auto dp = Q2.col(i).dot(G2 * Q2.col(i));
      Dt(i) = 1.0/sqrt(dp);
    }

    V0 = Q1 * D .asDiagonal();
    Vt = Q2 * Dt.asDiagonal();
    

    Eigen::MatrixXd Gq = Vt.transpose() * G2 * Vt;
    Eigen::MatrixXd Gp = V0.transpose() * G1 * V0;

    // Compute SVD of PP
    //
    auto PP = V0.inverse() * Vt;

    if (use_red) {
      RedSVD::RedSVD<Eigen::MatrixXd> svd(PP, evCount);
      UU = svd.matrixU();
      SS = svd.singularValues();
      VV = svd.matrixV();
    } else {
      Eigen::BDCSVD<Eigen::MatrixXd>
	svd(PP, Eigen::ComputeThinU | Eigen::ComputeThinV);
      if (svd.info() != Eigen::Success) {
	throw std::runtime_error("LiouvilleRKHS::singular_compute: SVD failure");
      }
      UU = svd.matrixU();
      SS = svd.singularValues();
      VV = svd.matrixV();
    }

    // Compute the projection of the state into the eigenbasis
    //
    Dp.resize(traj);
    Dq.resize(traj);

#pragma omp parallel for
    for (int i=0; i<traj; i++) {
      auto dp = UU.col(i).dot(Gp * UU.col(i));
      Dp(i) = 1.0/sqrt(dp);
    }

#pragma omp parallel for
    for (int i=0; i<traj; i++) {
      auto dp = VV.col(i).dot(Gq * VV.col(i));
      Dq(i) = 1.0/sqrt(dp);
    }
    
    Eigen::MatrixXd T = trajectory();

    Xh   = Dp.asDiagonal() * UU.transpose() * V0.transpose() * T;

    computed      = true;
    reconstructed = false;

#ifdef TESTING
    std::ofstream fout("LiouvilleTest.S");
    if (fout) {
      auto header = [&fout] (const std::string& label) {
	fout << std::string(80, '-') << std::endl
	     << "---- " << label << std::endl
	     << std::string(80, '-') << std::endl;
      };

      auto symcheck = [&fout] (const Eigen::MatrixXd& A) {
	Eigen::MatrixXd diff = A - A.transpose();
	double cmin =  std::numeric_limits<double>::infinity();
	double cmax = -std::numeric_limits<double>::infinity();
	auto d = diff.data();
	for (int i=0; i<diff.size(); i++) {
	  cmin = std::min(cmin, d[i]);
	  cmax = std::max(cmax, d[i]);
	}
	fout << "Min: " << cmin << " Max: " << cmax << std::endl;
	if (fabs(cmin) < 1e-10 and fabs(cmax) < 1e-10)
	  fout << "SYMMETRIC" << std::endl;
	else
	  fout << "NOT symmetric" << std::endl;
      };

      auto zerocheck = [&fout](const Eigen::MatrixXcd& M)
      {
	double maxval = 0.0;
	for (int i=0; i<M.rows(); i++) {
	  for (int j=0; j<M.cols(); j++) {
	    maxval = std::max(maxval, std::abs(M(i, j)));
	  }
	}
	fout << "max value=" << maxval << std::endl;
	if (maxval < 1e-10)
	  fout << "ZERO" << std::endl;
	else
	  fout << "NOT zero" << std::endl;
      };

      header("Parameters");    fout << "Traj=" << traj
				    << " Rank=" << rank
				    << " Ntimes=" << numT
				    << " method=" << Method_names[method]
				    << std::endl;

      header("Grammian 1");    fout << G1   << std::endl;
      header("Grammian 2");    fout << G2   << std::endl;
      header("Grammian SV 1"); fout << S1   << std::endl;
      header("Grammian SV 2"); fout << S2   << std::endl;
      header("Traj kernel");   fout << T    << std::endl;
      header("PP");            fout << PP   << std::endl;
      header("S");             fout << SS   << std::endl;
      header("U");             fout << UU   << std::endl;
      header("V");             fout << VV   << std::endl;
      header("Dq");            fout << Dq   << std::endl;
      header("Dp");            fout << Dp   << std::endl;
      header("Xh"  );          fout << Xh   << std::endl;

      header("Is G1 symmetric? Expect YES" ); symcheck(G1  );
      header("Is G2 symmetric? Expect YES" ); symcheck(G2  );
      header("Is PP symmetric? Expect NO"  ); symcheck(PP  );
    } else {
      std::cerr << "LiouvilleRKHS: Failed to open output file" << std::endl;
    }
#endif
  }

  // Algorithm 7.1 from Rosenfeld and Kamalapurkhar
  //
  void LiouvilleRKHS::eigenfunction_analysis()
  {
    // Number of channels
    //
    auto keys = getAllKeys();
    nkeys = keys.size();

    // The number of time points
    //
    numT = data[keys[0]].rows();
    rank = data[keys[0]].cols();

    // Get Grammian matrices
    //
    G1 = computeGrammian(mu1);
    G2 = computeGrammian(mu2);
    G3 = computeGrammian(mu2, mu1/mu2);
    A  = computeGammaDiff(mu1);

    // For regularization
    //
    auto R = Eigen::MatrixXd::Identity(traj, traj) * eps;

    // Perform eigenanalysis of Grammian matrix (pos def)
    //
    if (use_red) {
      RedSVD::RedSymEigen<Eigen::MatrixXd> eigensolver1(G2 + R, evCount);
      S2 = eigensolver1.eigenvalues();
      Q2 = eigensolver1.eigenvectors();

      RedSVD::RedSymEigen<Eigen::MatrixXd> eigensolver2(G3 + R, evCount);
      S3 = eigensolver2.eigenvalues();
      Q3 = eigensolver2.eigenvectors();
    } else {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver1(G2 + R);
      if (eigensolver1.info() != Eigen::Success) {
	throw std::runtime_error("LiouvilleRKHS: Eigensolver for Grammian 2 failed");
      }
      S2 = eigensolver1.eigenvalues();
      Q2 = eigensolver1.eigenvectors();

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver2(G3 + R);
      if (eigensolver2.info() != Eigen::Success) {
	throw std::runtime_error("LiouvilleRKHS: Eigensolver for Grammian 2 failed");
      }
      S3 = eigensolver2.eigenvalues();
      Q3 = eigensolver2.eigenvectors();
    }
    
    // Compute pseudoinvese of G2
    //
    Eigen::MatrixXd G2inv(traj, traj);
    Eigen::VectorXd S2inv(traj);

    double TOL = S2(S2.size()-1) * tol;

    for (int n=std::max<int>(0, S2.size()-evCount); n<S2.size(); n++) {
      if (S2(n) < TOL) {
	S2inv(n) = 0.0;
      } else {
	S2inv(n) = 1.0 / S2(n);
      }
    }

    // Compute pseudoinvese of G3
    //
    Eigen::MatrixXd G3inv(traj, traj);
    Eigen::VectorXd S3inv(traj);

    TOL = S3(S3.size()-1) * tol;

    for (int n=std::max<int>(0, S3.size()-evCount); n<S3.size(); n++) {
      if (S3(n) < TOL) {
	S3inv(n) = 0.0;
      } else {
	S3inv(n) = 1.0 / S3(n);
      }
    }

    Eigen::MatrixXd PP = (Q3*S3inv.asDiagonal()*Q3.transpose()) * G1 *
      (Q2*S2inv.asDiagonal()*Q2.transpose()) * A;


    // Perform eigenanalysis for PP
    //
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(PP);

    // Eigenvalues and eigenvectors of PP
    //
    L = eigensolver.eigenvalues();
    V = eigensolver.eigenvectors();

    // Compute the normalized eigenbasis
    //
    Vbar = V;
    for (int i=0; i<V.cols(); i++) {
      auto dp = V.col(i).dot(Eigen::MatrixXcd(G3) * V.col(i));
      Vbar.col(i) /= sqrt(dp.real());
    }

    // Compute the projection of the state into the eigenbasis
    //
    O    = occupation();
    Ginv = (Vbar.transpose() * G3 * Vbar).inverse();
    Xi   = Ginv * Vbar.transpose() * O;
    // Ginv = (Vbar.adjoint() * G3 * Vbar).inverse();
    // Xi   = Ginv * Vbar.adjoint() * O;

    computed      = true;
    reconstructed = false;

#ifdef TESTING
    std::ofstream fout("LiouvilleTest.E");
    if (fout) {
      auto header = [&fout] (const std::string& label) {
	fout << std::string(80, '-') << std::endl
	     << "---- " << label << std::endl
	     << std::string(80, '-') << std::endl;
      };

      auto symcheck = [&fout] (const Eigen::MatrixXd& A) {
	Eigen::MatrixXd diff = A - A.transpose();
	double cmin =  std::numeric_limits<double>::infinity();
	double cmax = -std::numeric_limits<double>::infinity();
	auto d = diff.data();
	for (int i=0; i<diff.size(); i++) {
	  cmin = std::min(cmin, d[i]);
	  cmax = std::max(cmax, d[i]);
	}
	fout << "Min: " << cmin << " Max: " << cmax << std::endl;
	if (fabs(cmin) < 1e-10 and fabs(cmax) < 1e-10)
	  fout << "SYMMETRIC" << std::endl;
	else
	  fout << "NOT symmetric" << std::endl;
      };

      auto zerocheck = [&fout](const Eigen::MatrixXcd& M)
      {
	double maxval = 0.0;
	for (int i=0; i<M.rows(); i++) {
	  for (int j=0; j<M.cols(); j++) {
	    maxval = std::max(maxval, std::abs(M(i, j)));
	  }
	}
	fout << "max value=" << maxval << std::endl;
	if (maxval < 1e-10)
	  fout << "ZERO" << std::endl;
	else
	  fout << "NOT zero" << std::endl;
      };

      header("Parameters");    fout << "Traj=" << traj
				    << " Rank=" << rank
				    << " Ntimes=" << numT
				    << " method=" << Method_names[method]
				    << std::endl;

      header("Grammian 1");    fout << G1   << std::endl;
      header("Grammian 2");    fout << G2   << std::endl;
      header("Grammian 3");    fout << G3   << std::endl;
      header("Occ kernel");    fout << A    << std::endl;
      header("Grammian 2 SV"); fout << S2   << std::endl;
      header("Grammian 3 SV"); fout << S3   << std::endl;
      header("PP");            fout << PP   << std::endl;
      header("Lambda");        fout << L    << std::endl;
      header("V");             fout << V    << std::endl;
      header("Vbar");          fout << Vbar << std::endl;
      header("Ginv");          fout << Ginv << std::endl;
      header("O"   );          fout << O    << std::endl;
      header("Xi"  );          fout << Xi   << std::endl;

      auto tst1 = Vbar.transpose() * G3 * Vbar;
      auto tst2 = Vbar.adjoint()   * G3 * Vbar;

      header("Gram Vbar");     fout << tst1 << std::endl;
      header("Gram Vbar 2");   fout << tst2 << std::endl;

      header("Is G1 symmetric? Expect YES" ); symcheck(G1 );
      header("Is G2 symmetric? Expect YES" ); symcheck(G2 );
      header("Is G3 symmetric? Expect YES" ); symcheck(G3 );
      header("Is A symmetric?  Expect NO"  ); symcheck(A  );
      header("Is PP symmetric? Expect NO"  ); symcheck(PP );

      header("Good eigenanalysis? Expect YES");
      zerocheck(V*L.asDiagonal()*V.inverse() - PP);
    } else {
      std::cerr << "LiouvilleRKHS: Failed to open output file" << std::endl;
    }
#endif
  }

  Eigen::VectorXd LiouvilleRKHS::flow(const Eigen::VectorXd& x)
  {
    if (method != Method::Singular)
      throw std::runtime_error("LiouvilleRKHS::flow is defined for the Singular method only");

    Eigen::VectorXd Psi = Dq.asDiagonal()*VV.transpose()*Vt.transpose() *
      computeGamma(x, mu1);
    
    return Xh.transpose() * SS.asDiagonal() * Psi;
  }

  Eigen::MatrixXd LiouvilleRKHS::computeTrajectory
  (const Eigen::VectorXd& x, double h)
  {
    if (method != Method::Singular)
      throw std::runtime_error("LiouvilleRKHS::computeTrajectory only makes sense for the Singular method");

    // Define an ODE integrator for a time-independent flow (RK4 for now)
    auto ODE = [this](const Eigen::VectorXd& x, double t, double h, int substeps)
    {
      double hh = h/substeps;
      Eigen::VectorXd y  = x;

      for (int i=0; i<substeps; i++) {

	Eigen::VectorXd k1 = hh * this->flow(y         );
	Eigen::VectorXd k2 = hh * this->flow(y + 0.5*k1);
	Eigen::VectorXd k3 = hh * this->flow(y + 0.5*k2);
	Eigen::VectorXd k4 = hh * this->flow(y + k3    );
	
	y += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
      }

      return std::tuple<Eigen::VectorXd, double>(y, t + h);
    };

#ifdef TESTING
    std::cout << "Computing trajectory with " << numT << " X " << rank
	      << " array" << std::endl;
#endif
    
    // An array to contain the returned trajectory
    Eigen::MatrixXd ret(numT, rank);

    // Compute the time step
    double dt = coefDB.times[1] - coefDB.times[0];

    // Compute the substep number
    int substeps = std::max<int>(1, floor(dt/h));

    // Initial time
    double t = coefDB.times[0];

    // Initial state vector
    Eigen::VectorXd y = x;
    ret.row(0) = y;

    // Integrate the flow
    for (int i=1; i<numT; i++) {
      std::tie(y, t) = ODE(y, t, dt, substeps);
      ret.row(i) = y;
    }
    
    return ret;
  }


  const std::set<std::string>
  LiouvilleRKHS::valid_keys = {
    "mu1",
    "mu2",
    "eps",
    "d",
    "alpha",
    "use_red",
    "kernel",
    "method",
    "verbose",
    "output"
  };

  template <typename Derived>
  std::string get_shape(const Eigen::EigenBase<Derived>& x)
  {
    std::ostringstream oss;
    oss  << "(" << x.rows() << ", " << x.cols() << ")";
    return oss.str();
  }

  Eigen::VectorXcd LiouvilleRKHS::evecEval(const Eigen::VectorXd& x)
  {
    if (method != Method::Eigenfunction)
      throw std::runtime_error("LiouvilleRKHS::evecEval make sense for the Eigenfunction method only");

    if (not computed) analysis();

    return Vbar.transpose() * computeGamma(x, mu1);
  }

  void LiouvilleRKHS::assignParameters(const std::string flags)
  {
    // Default parameter values
    //
    d        = 1.0;
    alpha    = 10.0;
    verbose  = false;

    // Parse the parameters database
    //
    try {

      // Load parameters from string
      //
      params = YAML::Load(flags);

      // Check for unmatched keys
      //
      auto unmatched = YamlCheck(params, valid_keys);
      if (unmatched.size())
	throw YamlConfigError("MSSA::LiouvilleRKHS", "parameter", unmatched, __FILE__, __LINE__);

      // Compute flags
      //
      computed      = false;
      reconstructed = false;

      // Top level parameter flags
      //
      if (params["d"])
	d  = double(params["d"].as<double>());

      if (params["eps"])
	eps = double(params["eps"].as<double>());

      if (params["mu1"])
	mu1 = double(params["mu1"].as<double>());

      if (params["mu2"])
	mu2 = double(params["mu2"].as<double>());

      if (params["use_red"])
	use_red = bool(params["use_red"].as<bool>());

      if (params["alpha"])
	alpha = double(params["alpha"].as<double>());

      if (params["kernel"]) {
	std::string type = params["kernel"].as<std::string>();
	if (RKHS_values.find(type) != RKHS_values.end())
	  rkhs = RKHS_values[type];
	else
	  throw std::runtime_error("LiouvilleRKHS: kernel type [" + type +
				   "] not found");
      }
      
      if (params["method"]) {
	std::string type = params["method"].as<std::string>();
	if (Method_values.find(type) != Method_values.end())
	  method = Method_values[type];
	else
	  throw std::runtime_error("LiouvilleRKHS: method type [" + type +
				   "] not found");
      }
      
      if (params["verbose"])
	verbose  = bool(params["verbose"]);

      if (params["output"]) prefix = params["output"].as<std::string>();
      else                  prefix = "exp_rkhs";

    }
    catch (const YAML::ParserException& e) {
      std::cout << "LiouvilleRKHS::assignParameters, parsing error=" << e.what()
		<< std::endl;
      throw;
    }
  }


  // Save current KOOPMAN state to an HDF5 file with the given prefix
  void LiouvilleRKHS::saveState(const std::string& prefix)
  {
    if (not computed) return;	// No point in saving anything

    if (std::filesystem::exists(prefix + "_rkhs.h5")) {
      std::ostringstream sout;
      sout << "LiouvilleRKHS::saveState: the file <" <<  prefix + "_edmd.h5"
	   << "> already exists.\nPlease delete this file or choose a "
	   << "different file name";
      throw std::runtime_error(sout.str());
    }

    try {
      // Create a new hdf5 file
      //
      HighFive::File file(prefix + "_rkhs.h5",
			  HighFive::File::ReadWrite |
			  HighFive::File::Create);

      // Write the time dimension
      //
      file.createAttribute<int>("numT", HighFive::DataSpace::From(numT)).write(numT);

      // Write the number of channels
      //
      file.createAttribute<int>("nKeys", HighFive::DataSpace::From(nkeys)).write(nkeys);

      // Write the number of eigenvalues
      //
      file.createAttribute<int>("nEV", HighFive::DataSpace::From(nev)).write(nev);

      // Hilbert space domain parameter
      //
      file.createAttribute<double>("mu1", HighFive::DataSpace::From(mu1)).write(mu1);

      // Hilbert space range parameter
      //
      file.createAttribute<double>("mu2", HighFive::DataSpace::From(mu2)).write(mu2);

      // Regularization parameter
      //
      file.createAttribute<double>("eps", HighFive::DataSpace::From(eps)).write(eps);

      // The algorithm method
      //
      std::string mthd = Method_names[method];
      file.createAttribute<std::string>("method", HighFive::DataSpace::From(mthd)).write(mthd);

      // Save the key list
      //
      std::vector<Key> keylist;
      for (auto k : data) keylist.push_back(k.first);

      // Pad keylist entries to maximum key length for HighFive
      //
      size_t maxSZ = 0;
      for (auto v : keylist) maxSZ = std::max<size_t>(maxSZ, v.size());

      // The padding value
      //
      auto padVal = std::numeric_limits<unsigned>::max();
      for (auto & v : keylist) {
	if (v.size() < maxSZ) {
	  for (auto k=v.size(); k<maxSZ; k++) v.push_back(padVal);
	}
      }

      // Finally, create the dataset
      //
      file.createDataSet("keylist", keylist);

      // Save koopman_analysis state
      //
      HighFive::Group analysis = file.createGroup("koopmanRKHS_analysis");

      if (method == Method::Eigenfunction) {

	analysis.createDataSet("rkhs",  RKHS_names[rkhs]);
	analysis.createDataSet("G1",    G1  );
	analysis.createDataSet("G2",    G2  );
	analysis.createDataSet("G3",    G3  );
	analysis.createDataSet("A",     A   );
	analysis.createDataSet("Vbar",  Vbar);
	analysis.createDataSet("L",     L   );
	analysis.createDataSet("Xi",    Xi  );
	analysis.createDataSet("Phi",   Phi );
	analysis.createDataSet("S2",    S2  );
	analysis.createDataSet("S3",    S3  );
	analysis.createDataSet("Q2",    Q2  );
	analysis.createDataSet("Q3",    Q3  );
      } else {
	analysis.createDataSet("rkhs",  RKHS_names[rkhs]);
	analysis.createDataSet("G1",    G1  );
	analysis.createDataSet("G2",    G2  );
	analysis.createDataSet("S1",    S1  );
	analysis.createDataSet("S2",    S2  );
	analysis.createDataSet("Q1",    Q1  );
	analysis.createDataSet("Q2",    Q2  );
	analysis.createDataSet("V0",    V0  );
	analysis.createDataSet("Vt",    Vt  );
	analysis.createDataSet("U",     UU  );
	analysis.createDataSet("V",     VV  );
	analysis.createDataSet("S",     SS  );
	analysis.createDataSet("Dp",    Dp  );
	analysis.createDataSet("Dq",    Dq  );
	analysis.createDataSet("Xh",    Xh  );
      }

    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }
  }

  // Restore current KOOPMAN state to an HDF5 file with the given prefix
  void LiouvilleRKHS::restoreState(const std::string& prefix)
  {
    try {
      // Silence the HDF5 error stack
      //
      HighFive::SilenceHDF5 quiet;

      // Try opening the file as HDF5
      //
      HighFive::File h5file(prefix + "_rkhs.h5", HighFive::File::ReadOnly);

      // Read and check parameters
      //
      int nTime, nKeys;
      double MU1, MU2, EPS;
      std::string mthd;

      h5file.getAttribute("numT"  ).read(nTime);
      h5file.getAttribute("nKeys" ).read(nKeys);
      h5file.getAttribute("nEV"   ).read(nev  );
      h5file.getAttribute("mu1"   ).read(MU1  );
      h5file.getAttribute("mu2"   ).read(MU2  );
      h5file.getAttribute("eps"   ).read(EPS  );
      h5file.getAttribute("method").read(mthd);


      // Number of channels
      //
      nkeys = data.size();

      // Test recovered parameters
      //
      if (nTime != numT) {
	std::ostringstream sout;
	sout << "LiouvilleRKHS::restoreState: saved state has numT="
	     << nTime << " but LiouvilleRKHS expects numT=" << numT
	     << ".\nCan't restore LiouvilleRKHS state!";
	throw std::runtime_error(sout.str());
      }

      if (nKeys != nkeys) {
	std::ostringstream sout;
	sout << "LiouvilleRKHS::restoreState: saved state has nkeys="
	     << nKeys << " but LiouvilleRKHS expects nkeys=" << nkeys
	     << ".\nCan't restore LiouvilleRKHS state!";
	throw std::runtime_error(sout.str());
      }

      if (nKeys != nkeys) {
	std::ostringstream sout;
	sout << "LiouvilleRKHS::restoreState: saved state has nkeys="
	     << nKeys << " but LiouvilleRKHS expects nkeys=" << nkeys
	     << ".\nCan't restore LiouvilleRKHS state!";
	throw std::runtime_error(sout.str());
      }

      if (fabs(mu1-MU1) > 1.0e-12) {
	std::ostringstream sout;
	sout << "LiouvilleRKHS::restoreState: saved state has mu1="
	     << MU1 << " but LiouvilleRKHS expects mu1=" << mu1
	     << ".\nCan't restore LiouvilleRKHS state!";
	throw std::runtime_error(sout.str());
      }

      if (fabs(mu2-MU2) > 1.0e-12) {
	std::ostringstream sout;
	sout << "LiouvilleRKHS::restoreState: saved state has mu2="
	     << MU2 << " but LiouvilleRKHS expects mu2=" << mu2
	     << ".\nCan't restore LiouvilleRKHS state!";
	throw std::runtime_error(sout.str());
      }

      if (fabs(eps-EPS) > 1.0e-18) {
	std::ostringstream sout;
	sout << "LiouvilleRKHS::restoreState: saved state has eps="
	     << EPS << " but LiouvilleRKHS expects eps=" << eps
	     << ".\nCan't restore LiouvilleRKHS state!";
	throw std::runtime_error(sout.str());
      }

      if (Method_values[mthd] != method) {
	std::ostringstream sout;
	sout << "LiouvilleRKHS::restoreState: saved state has method="
	     << mthd << " but LiouvilleRKHS expects method=" << Method_names[method]
	     << ".\nCan't restore LiouvilleRKHS state!";
	throw std::runtime_error(sout.str());
      }

      std::vector<Key> keylist;
      h5file.getDataSet("keylist").read(keylist);

      // Remove padded values from K5 store
      //
      auto padVal = std::numeric_limits<unsigned>::max();
      for (auto & v : keylist) {
	std::vector<unsigned int>::iterator it;
	while ((it = std::find(v.begin(), v.end(), padVal)) != v.end())
	  v.erase(it);
      }

      // Check key list
      //
      bool bad = false;
      for (int n=0; n<keylist.size(); n++) {
	auto it = data.find(keylist[n]);
	if (it == data.end()) bad = true;
	else if (it->first.size() != keylist[n].size()) bad = true;
	else {
	  for (int i=0; i<keylist[n].size(); i++) {
	    if (keylist[n][i] != it->first[i]) bad = true;
	  }
	}
      }

      if (bad) {
	std::ostringstream sout;
	sout << "LiouvilleRKHS::restoreState: keylist mismatch." << std::endl
	     << "Can't restore LiouvilleRKHS state! Wanted keylist: ";
	for (auto v : data) {
	  sout << "[";
	  for (auto u : v.first) sout << u << ' ';
	  sout << "] ";
	}
	sout << std::endl << "Found keylist: ";
	for (auto v : keylist) {
	  sout << "[";
	  for (auto u : v) sout << u << ' ';
	  sout << "] ";
	}

	throw std::runtime_error(sout.str());
      }

      auto analysis = h5file.getGroup("koopmanRKHS_analysis");

      std::string type = analysis.getDataSet("rkhs").read<std::string>();
      rkhs = RKHS_values[type];

      if (method == Method::Eigenfunction) {

	G1   = analysis.getDataSet("G1"  ).read<Eigen::MatrixXd >();   
	G2   = analysis.getDataSet("G2"  ).read<Eigen::MatrixXd >();   
	G3   = analysis.getDataSet("G3"  ).read<Eigen::MatrixXd >();   
	A    = analysis.getDataSet("A"   ).read<Eigen::MatrixXd >();   
	Vbar = analysis.getDataSet("Vbar").read<Eigen::MatrixXcd>(); 
	L    = analysis.getDataSet("L"   ).read<Eigen::VectorXcd>();  
	Xi   = analysis.getDataSet("Xi"  ).read<Eigen::MatrixXcd>();   
	Phi  = analysis.getDataSet("Phi" ).read<Eigen::MatrixXcd>(); 
	S2   = analysis.getDataSet("S2"  ).read<Eigen::VectorXd >();   
	S3   = analysis.getDataSet("S3"  ).read<Eigen::VectorXd >();   
	Q2   = analysis.getDataSet("Q2"  ).read<Eigen::MatrixXd >();   
	Q3   = analysis.getDataSet("Q3"  ).read<Eigen::MatrixXd >();   

      } else {

	G1   = analysis.getDataSet("G1"  ).read<Eigen::MatrixXd>();
	G2   = analysis.getDataSet("G2"  ).read<Eigen::MatrixXd>();
	S1   = analysis.getDataSet("S1"  ).read<Eigen::VectorXd>();
	S2   = analysis.getDataSet("S2"  ).read<Eigen::VectorXd>();
	Q1   = analysis.getDataSet("Q1"  ).read<Eigen::MatrixXd>();
	Q2   = analysis.getDataSet("Q2"  ).read<Eigen::MatrixXd>();
	V0   = analysis.getDataSet("V0"  ).read<Eigen::MatrixXd>();
	Vt   = analysis.getDataSet("Vt"  ).read<Eigen::MatrixXd>();
	UU   = analysis.getDataSet("U"   ).read<Eigen::MatrixXd>();
	VV   = analysis.getDataSet("V"   ).read<Eigen::MatrixXd>();
	SS   = analysis.getDataSet("S"   ).read<Eigen::VectorXd>();
	Dp   = analysis.getDataSet("Dp"  ).read<Eigen::VectorXd>();
	Dq   = analysis.getDataSet("Dq"  ).read<Eigen::VectorXd>();
	Xh   = analysis.getDataSet("Xh"  ).read<Eigen::MatrixXd>();

      }

      computed = true;

    } catch (HighFive::Exception& err) {
      std::cerr << "**** Error opening or reading H5 file ****" << std::endl;
      throw;
    }

  }

  LiouvilleRKHS::LiouvilleRKHS(const mssaConfig& config, double tol, int count,
			   const std::string flags) : tol(tol), evCount(count)
  {
    // Parse the YAML string
    //
    assignParameters(flags);

    // Eigen OpenMP reporting
    //
    static bool firstTime = true;
    if (firstTime) {
      std::cout << "---- Eigen is using " << Eigen::nbThreads()
		<< " threads" << std::endl;
      firstTime = false;
    }

    // Now open and parse the coefficient files
    //
    coefDB = CoefContainer(config, flags);

    // Number of snapshots
    //
    numT = coefDB.times.size();

    // Database names; there should only be one
    //
    auto names = coefDB.getNames();

    // Get coefficients for this database
    //
    auto coefs = coefDB.getCoefs(names[0]);

    // Get the first coefficient set for parameter reflection
    //
    auto cf = std::dynamic_pointer_cast<CoefClasses::TrajStruct>
      (coefs->getCoefStruct(coefDB.times[0]));

    traj = cf->traj;
    rank = cf->rank;

    // Generate data store
    //
    auto keys = coefDB.getKeys(); // For future generalization...

    // Allocate storage and reorder coefficients
    //
    for (auto key : coefDB.getKeys()) {
      data[key].resize(numT, rank);
      auto v = coefDB.getData(key);
      for (int t=0; t<numT; t++) {
	for (int n=0; n<rank; n++) data[key](t, n) = v[t*rank + n];
      }
    }
    
    computed      = false;
    reconstructed = false;
  }
  // END LiouvilleRKHS constructor

}
// END namespace MSSA
