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
    traj = keys.size();
    
    // Check for odd number of times for Simpson's 1/3 rule
    //
    int nt = numT;
    if (nt/2*2 == nt) nt -= 1;

    // Allocate Grammian matrix and set to zero
    //
    Eigen::MatrixXd G(traj, traj);
    int f1, f2;

    for (int i=0; i<traj; i++) {
      for (int j=0; j<traj; j++) {

	G(i, j) = 0.0;

	for (int t1=0; t1<nt; t1++) {

	  if (t1 == 0 or t1 == nt-1) f1 = 1;
	  else {
	    if (t1 == 1) f1 = 4;
	    if (t1 == 4) f1 = 2;
	  }

	  double inner = 0.0;
	  for (int t2=0; t2<nt; t2++) {
	    
	    if (t2 == 0 or t2 == nt-1) f2 = 1;
	    else {
	      if (t2 == 1) f2 = 4;
	      if (t2 == 4) f2 = 2;
	    }

	    auto x1 = data[t1].row(keys[i])*rat;
	    auto x2 = data[t2].row(keys[j])*rat;
	    inner += kernel(data[t1].row(keys[i])*rat, data[t2].row(keys[j])*rat, mu) * f2;
	  }

	  G(i, j) += inner*f1;
	}
      }
    }

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
    traj = keys.size();
    
    // Check for odd number of times for Simpson's 1/3 rule
    //
    int nt = numT;
    if (nt/2*2 == nt) nt -= 1;

    // Allocate Grammian matrix and set to zero
    //
    Eigen::MatrixXd A(traj, traj);
    int f;

    for (int i=0; i<traj; i++) {
      for (int j=0; j<traj; j++) {

	A(i, j) = 0.0;

	double sumT = 0.0, sum0 = 0.0;

	for (int t=0; t<nt; t++) {

	  if (t == 0 or t == nt-1) f = 1;
	  else {
	    if (t == 1) f = 4;
	    if (t == 4) f = 2;
	  }

	  sum0 += kernel(data[0   ].row(keys[i]), data[t].row(keys[j]), mu) * f;
	  sumT += kernel(data[nt-1].row(keys[i]), data[t].row(keys[j]), mu) * f;
	}

	A(i, j) += sumT - sum0;
      }
    }

    // Time step factor
    //
    double SF = (coefDB.times[1] - coefDB.times[0])/3.0;

    return A * SF;
  }


  // Algorithm 7.1 from Rosenfeld and Kamalapurkhar
  //
  void LiouvilleRKHS::analysis()
  {
    // Number of channels
    //
    auto keys = getAllKeys();
    nkeys = keys.size();

    // The number of time points
    //
    numT = data[keys[0]].size();

    // Get Grammian matrices
    //
    G1 = computeGrammian(mu1);
    G2 = computeGrammian(mu2);
    G3 = computeGrammian(mu2, mu1/mu2);
    A  = computeGammaDiff(mu1);

    // Perform eigenanalysis of Grammian matrix (pos def)
    //
    if (use_red) {
      RedSVD::RedSymEigen<Eigen::MatrixXd> eigensolver1(G2, evCount);
      S2 = eigensolver1.eigenvalues();
      Q3 = eigensolver1.eigenvectors();

      RedSVD::RedSymEigen<Eigen::MatrixXd> eigensolver2(G3, evCount);
      S3 = eigensolver2.eigenvalues();
      Q3 = eigensolver2.eigenvectors();
    } else {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver1(G2);
      if (eigensolver1.info() != Eigen::Success) {
	throw std::runtime_error("LiouvilleRKHS: Eigensolver for Grammian 2 failed");
      }
      S2 = eigensolver1.eigenvalues();
      Q2 = eigensolver1.eigenvectors();

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver2(G3);
      if (eigensolver2.info() != Eigen::Success) {
	throw std::runtime_error("LiouvilleRKHS: Eigensolver for Grammian 2 failed");
      }
      S3 = eigensolver2.eigenvalues();
      Q3 = eigensolver2.eigenvectors();
    }
    
    // Compute pseudoinvese of G2
    //
    Eigen::MatrixXd G2inv(traj, traj);
    Eigen::Vectorxd S2inv(traj);

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
    Eigen::Vectorxd S3inv(traj);

    for (int n=std::max<int>(0, S3.size()-evCount); n<S3.size(); n++) {
      if (S3(n) < TOL) {
	S3inv(n) = 0.0;
      } else {
	S3inv(n) = 1.0 / S3(n);
      }
    }

    Eigen::MatrixXd PPG = Q2.transpose()*S2inv.asDiagonal()*Q2 *G1 *
      Q3.tranpose()*S3inv.asDiagonal*Q3 * A;


    // Perform eigenanalysis for PPG
    //
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(PPG);

    // Eigenvalues and eigenvectors of PPG
    //
    Eigen::VectorXcd lam = eigensolver.eigenvalues();
    Eigen::MatrixXcd V   = eigensolver.eigenvectors();

    // Use equation 7.2 to compute eigenfunctions
    //


    // Use equation 7.3 to compute modes
    //

#ifdef TESTING
    std::ofstream file("testing.dat");
    if (file) {
      file << std::string(70, '-') << std::endl;
      file << "---- Left/right eigenvalues difference" << std::endl;
      file << std::string(70, '-') << std::endl;
      for (int i=0; i<SL.size(); i++)
	file << std::setw(5) << i
	     << " " << std::setw(18) << std::abs(SL(i) - SR(i))
	     << " " << SL(i)
	     << " " << SR(i)
	     << std::endl;
      file << std::string(70, '-') << std::endl;
      file << std::string(70, '-') << std::endl;
      file << "---- L/R eigenvector test" << std::endl;
      file << std::string(70, '-') << std::endl;
      matrix_print(file, U.adjoint()*V);
      file << std::string(70, '-') << std::endl;
    }
#endif

#ifdef TESTING
    if (file) {
      file << "---- Swapped left eigenvalues" << std::endl;
      file << std::string(70, '-') << std::endl;
      for (int i=0; i<SL.size(); i++)
	file << std::setw(5) << i
	     << " " << std::setw(18) << std::abs(SL(i))
	     << " " << std::setw(18) << std::abs(SR(i))
	     << std::endl;
      file << std::string(70, '-') << std::endl;
      file << "---- Left normalization" << std::endl;
      file << std::string(70, '-') << std::endl;
    } else {
      std::cerr << "ERROR: Could not open testing.dat" << std::endl;
    }
#endif

    // Normalize left eigenvectors
    for (int i=0; i<U.rows(); i++) {
      auto nrm = U.adjoint().row(i) * V.col(i);
#ifdef TESTING
      if (file) file << std::setw(5) << i << ": " << nrm << std::endl;
#endif
      U.col(i) /= std::conj(nrm(0));
    }

#ifdef TESTING
    if (file) {
      file << std::string(70, '-') << std::endl;
      file << "---- Action matrix" << std::endl;
      file << "---- min/max = " << A.minCoeff() << "/" << A.maxCoeff() << std::endl;
      file << std::string(70, '-') << std::endl;
      file << A << std::endl;
      file << std::string(70, '-') << std::endl;
      file << "---- Grammian eigenvalues" << std::endl;
      file << std::string(70, '-') << std::endl;
      for (int i=0; i<S2.size(); i++) {
	file << std::setw(5) << i
		  << std::setw(18) << S2(i)
		  << std::setw(18) << SP(i) << std::endl;
      }
      file << std::string(70, '-') << std::endl;
      file << "---- Postnorm right/left eigenvector test" << std::endl;
      file << std::string(70, '-') << std::endl;
      if (matrix_idtest(file, U.adjoint()*V)) {
	file << std::string(70, '-') << std::endl;
	matrix_print(file, U.adjoint()*V);
      }
      file << std::string(70, '-') << std::endl;
    }
#endif

    // The Liouville modes from eq. 22
    //
    Xi = (U.adjoint()*SP.asDiagonal()*Q.transpose()*X).transpose();

    computed      = true;
    reconstructed = false;
    }

  const std::set<std::string>
  LiouvilleRKHS::valid_keys = {
    "mu1",
    "mu2"
    "eps",
    "d",
    "alpha",
    "use_red",
    "kernel",
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

  Eigen::VectorXcd LiouvilleRKHS::modeEval(int index, const Eigen::VectorXd& x)
  {
    if (not computed) analysis();
    
    Eigen::VectorXcd ret(nkeys);
    ret.setZero();

    if (index>=0 and index < numT) {
      Eigen::VectorXd Y(numT);
      for (int j=0; j<numT; j++) Y(j) = kernel(x, X.row(j));

      std::complex<double> Psi = Y.transpose()*Q*SP.asDiagonal()*V.col(index);

      ret = Xi.col(index)*Psi;
    }

    return ret;
  }

  std::complex<double> LiouvilleRKHS::evecEval(int index, const Eigen::VectorXd& x)
  {
    if (not computed) analysis();
    
    static bool first = true;
    if (first) {
      std::ofstream test("test.data");
      test << "Q"  << std::endl << Q  << std::endl;
      test << "SP" << std::endl << SP << std::endl;
      test << "V"  << std::endl << V  << std::endl;
      first = false;
    }

    std::complex<double> ret(0.0);

    if (index>=0 and index < numT) {
      Eigen::VectorXd Y(numT);
      for (int j=0; j<numT; j++) {
	Y(j) = kernel(x, X.row(j));
      }

      ret = Y.transpose()*Q*SP.asDiagonal()*V.col(index);
    }

    return ret;
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

      analysis.createDataSet("rkhs",  RKHS_names[rkhs]);
      analysis.createDataSet("Xi",    Xi );
      analysis.createDataSet("G",     G  );
      analysis.createDataSet("A",     A  );
      analysis.createDataSet("X",     X  );
      analysis.createDataSet("S2",    S2 );
      analysis.createDataSet("SP",    SP );
      analysis.createDataSet("SR",    SR );
      analysis.createDataSet("SL",    SL );
      analysis.createDataSet("Q",     Q  );
      analysis.createDataSet("K",     K  );
      analysis.createDataSet("U",     U  );
      analysis.createDataSet("V",     V  );
      analysis.createDataSet("nEV",   nev);

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

      h5file.getAttribute("numT" ).read(nTime);
      h5file.getAttribute("nKeys").read(nKeys);

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

      Xi  = analysis.getDataSet("Xi" ).read<Eigen::MatrixXd>();
      G   = analysis.getDataSet("G"  ).read<Eigen::MatrixXd >();
      A   = analysis.getDataSet("A"  ).read<Eigen::MatrixXd >();
      X   = analysis.getDataSet("X"  ).read<Eigen::MatrixXd >();
      S2  = analysis.getDataSet("S2" ).read<Eigen::VectorXd >();
      SP  = analysis.getDataSet("SP" ).read<Eigen::VectorXd >();
      SR  = analysis.getDataSet("SR" ).read<Eigen::VectorXd >();
      SL  = analysis.getDataSet("SL" ).read<Eigen::VectorXd >();
      Q   = analysis.getDataSet("Q"  ).read<Eigen::MatrixXd >();
      K   = analysis.getDataSet("K"  ).read<Eigen::MatrixXd >();
      U   = analysis.getDataSet("U"  ).read<Eigen::MatrixXcd>();
      V   = analysis.getDataSet("V"  ).read<Eigen::MatrixXcd>();
      nev = analysis.getDataSet("nev").read<int>();

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

    numT = coefDB.times.size();

    // Generate all the channels
    //
    for (auto key : coefDB.getKeys()) {
      data[key] = coefDB.getData(key);
    }

    computed      = false;
    reconstructed = false;
  }
  // END LiouvilleRKHS constructor

}
// END namespace MSSA
