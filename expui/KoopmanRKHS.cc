#define TESTING

//
// EDMD (Koopman theory) for EXP coefficients
//
// Uses fixed-rank approximation for the SVD to save time and space.
// Uses the approximate SVD randomized algorithm from Halko,
// Martinsson, and Tropp by default.  Use BDCSVD or Jacobi flags for
// the standard methods.  Using the randomized method together with
// trajectory matrix rather than covariance matrix analysis may lead
// to large errors in eigenvectors.
//
// The implementation here follows:
//
// M. O. Williams, C. W. Rowley, I. G. Kevrekidis, 2015, "A
// kernel-based method for data-driven Koopman spectral analysis",
// Journal of Computational Dynamics, 2 (2), 247-265
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

#include <KoopmanRKHS.H>

#include <RedSVD.H>
#include <YamlConfig.H>
#include <YamlCheck.H>
#include <EXPException.H>
#include <libvars.H>

namespace MSSA {

  // Helper ostream manipulator for debug coefficient key info
  std::ostream& operator<< (std::ostream& out, const std::vector<unsigned>& t);

  std::map<KoopmanRKHS::RKHS, std::string> KoopmanRKHS::RKHS_names =
    {
      {KoopmanRKHS::RKHS::Polynomial,  "Polynomial" },
      {KoopmanRKHS::RKHS::Exponential, "Exponential"},
      {KoopmanRKHS::RKHS::Gaussian,    "Gaussian"   }
    };

  std::map<std::string, KoopmanRKHS::RKHS> KoopmanRKHS::RKHS_values =
    {
      {"Polynomial",  KoopmanRKHS::RKHS::Polynomial },
      {"Exponential", KoopmanRKHS::RKHS::Exponential},
      {"Gaussian",    KoopmanRKHS::RKHS::Gaussian   }
    };


  //! RKHS kernel
  double KoopmanRKHS::kernel(const Eigen::VectorXd& x,
			     const Eigen::VectorXd& y)
  {
    if (rkhs == RKHS::Polynomial) {
      double prod = x.adjoint()*y;
      return pow(1.0 + prod/(d*d), alpha);
    } else if (rkhs == RKHS::Exponential) {
      double prod = x.adjoint()*y;
      return exp(prod/(d*d));
    } else if (rkhs == RKHS::Gaussian) {
      Eigen::VectorXd diff = x - y;
      double prod = diff.adjoint()*diff;
      return exp(-prod/(d*d));
    } else {
      throw std::runtime_error("KoopmanRKHS: Unknown kernel type");
    }
  }

  // Structure for sorting complex numbers and retrieving a permutation index
  using complex_elem = std::pair<std::complex<double>, int>;
  bool complex_comparator(const complex_elem &lhs, const complex_elem &rhs)
  {
    auto a = lhs.first;
    auto b = rhs.first;

    double abs_a = std::abs(a);
    double abs_b = std::abs(b);

    if (abs_a == abs_b) {
      if (std::real(a) == std::real(b))
	return std::imag(a) < std::imag(b);
      else
	return std::real(a) < std::real(b);
    } else {
      return abs_a < abs_b;
    }
  }

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

  bool matrix_idtest(std::ostream& out, const Eigen::MatrixXcd& M)
  {
    double max_diag = 0.0, max_offd = 0.0;
    for (int i=0; i<M.rows(); i++) {
      for (int j=0; j<M.cols(); j++) {
	if (i==j) max_diag = std::max(max_diag, std::abs(M(i, j)-1.0));
	else      max_offd = std::max(max_offd, std::abs(M(i, j)    ));
      }
    }
    out << std::string(70, '-') << std::endl << std::setprecision(6);
    out << "diag diff=" << max_diag << " offd diff=" << max_offd << std::endl;
    out << std::string(70, '-') << std::endl;

    if (max_diag>1.0e-3 or max_offd>1.0e-3) return true;
    return false;
  }


  // Algorithm 3 from Williams, Rowley, Kevrekidis
  //
  void KoopmanRKHS::analysis()
  {
    // Number of channels
    //
    auto keys = getAllKeys();
    nkeys = keys.size();

    numT = data[keys[0]].size();

    // Grammian and action matrices
    G.resize(numT, numT);
    A.resize(numT, numT);

    // Data matrix (rows = time, columns = keys)
    X.resize(numT, nkeys);
    
    for (int k=0; k<nkeys; k++) {
      for (int i=0; i<numT; i++) X(i, k) = data[keys[k]][i];
    }

    // Flow data
    Eigen::VectorXd Y(nkeys);

    // Time step
    double DT = coefDB.times[1] - coefDB.times[0];

    // Compute Grammian and action matrices
    for (int i=0; i<numT; i++) {

      for (int k=0; k<nkeys; k++) {
	if (i==0) {
	  Y(k) = (data[keys[k]][1] - data[keys[k]][0])/DT;
	} else if (i==numT-1) {
	  Y(k) = (data[keys[k]][numT-1] - data[keys[k]][numT-2])/DT;
	} else {
	  Y(k) = (data[keys[k]][i+1] - data[keys[k]][i-1])/(2.0*DT);
	}
      }

      if (testF) {
	if (kepler) {
	  for (int l=0; l<nkeys/4; l++) {
	    double x = data[keys[0+4*l]][i];
	    double y = data[keys[1+4*l]][i];
	    double u = data[keys[2+4*l]][i];
	    double v = data[keys[3+4*l]][i];
	      
	    double r2 = x*x + y*y;
	    double r  = sqrt(r2);
	    double dp = -1.0/r2;
	      
	    Y(0+4*l) = u;
	    Y(1+4*l) = v;
	    Y(2+4*l) = dp*x/r;
	    Y(3+4*l) = dp*y/r;
	  }
	}
	else if (plummer) {
	  if (radial) {
	    for (int l=0; l<nkeys/4; l++) {
	      double r  = data[keys[0+4*l]][i];
	      double pr = data[keys[2+4*l]][i];
	      double L  = data[keys[3+4*l]][i];

	      Y(0+4*l) = pr;
	      Y(1+4*l) = L/(r*r);
	      Y(2+4*l) = -r*pow(1.0+r*r, -1.5)/tscale;
	      Y(3+4*l) = 0.0;
	    }
	  } else {
	    for (int l=0; l<nkeys/4; l++) {
	      double x = data[keys[0+4*l]][i];
	      double y = data[keys[1+4*l]][i];
	      double u = data[keys[2+4*l]][i];
	      double v = data[keys[3+4*l]][i];
	      
	      double r = sqrt(x*x + y*y);
	      double dphi = -r*pow(1.0+r*r, -1.5)/tscale;
	      
	      Y(0+4*l) = u;
	      Y(1+4*l) = v;
	      Y(2+4*l) = dphi*x/r;
	      Y(3+4*l) = dphi*y/r;
	    }
	  }
	} else if (oscil) {
	  for (int l=0; l<nkeys/4; l++) {
	    double x0 = data[keys[0+4*l]][i];
	    double y0 = data[keys[1+4*l]][i];
	    double x1 = data[keys[2+4*l]][i];
	    double y1 = data[keys[3+4*l]][i];
	    Y(0+4*l) = -lam * y0;
	    Y(1+4*l) =  lam * x0;
	    Y(2+4*l) = -mu * y1 + c*(x0*x0 - y0*y0);
	    Y(3+4*l) =  mu * x1 + 2.0*c*x0*y0;
	  }
	} else {
	  for (int l=0; l<nkeys/2; l++) {
	    double x0 = data[keys[0+2*l]][i];
	    double x1 = data[keys[1+2*l]][i];
	    Y(0+2*l) = lam * x0;
	    Y(1+2*l) = mu  * x1 + c * x0 * x0;
	  }
	}
      }
	
      for (int j=0; j<numT; j++) {
	G(i, j) = kernel(X.row(i), X.row(j));
	A(i, j) = kernel(Y, X.row(j));
      }
    }

    // Perform eigenanalysis of Grammian matrix (pos def)
    //
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(G);
    if (eigensolver.info() != Eigen::Success) {
      throw std::runtime_error("KoopmanRKHS: Eigensolver for Grammian failed");
    }
    S2 = eigensolver.eigenvalues();
    Q  = eigensolver.eigenvectors();
    
    // Find inversion tolerance condition
    //
    SP.resize(S2.size());
    double TOL = S2(S2.size()-1) * tol;
    TOL = tol;
    int nev = 0;
    for (int n=0; n<S2.size(); n++) {
      if (S2(n) < TOL) {
	SP(n) = 0.0;
      } else {
	SP(n) = 1.0 / sqrt(S2(n));
	nev++;
      }
    }

    // Get Koopman operator estimate
    //
    K = (SP.asDiagonal() * Q.transpose()) * A * (Q * SP.asDiagonal());

    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver2(K);
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver3(K.adjoint());

    // Right and left eigenvectors
    //
    Eigen::VectorXcd SR2 = eigensolver2.eigenvalues();
    Eigen::MatrixXcd V2  = eigensolver2.eigenvectors();

    Eigen::VectorXcd SL2 = eigensolver3.eigenvalues().conjugate();
    Eigen::MatrixXcd U2  = eigensolver3.eigenvectors();

    U  = U2;
    V  = V2;

    SL = SL2;
    SR = SR2;

    // Reorder to match eigenvalues in left and right eigenvectors
    // using a comparison functor which sorts by magnitude, real
    // value, and imaginary value by precidence
    {
      // Right eigenvectors
      std::vector<complex_elem> v(SR2.size());
      for (int i=0; i<SR2.size(); i++) v[i] = complex_elem(SR2(i), i);
      std::sort(v.begin(), v.end(), complex_comparator);
      for (int i=0; i<SR2.size(); i++) {
	V.col(i) = V2.col(v[i].second);
	SR(i)    = SR2(v[i].second);
      }

      // Left eigenvectors
      std::vector<complex_elem> u(SL2.size());
      for (int i=0; i<SL2.size(); i++) u[i] = complex_elem(SL2(i), i);
      std::sort(u.begin(), u.end(), complex_comparator);
      for (int i=0; i<SL2.size(); i++) {
	U.col(i) = U2.col(u[i].second);
	SL(i)    = SL2(u[i].second);
      }
    }

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

    Xi = (U.adjoint()*SP.asDiagonal()*Q.transpose()*X).transpose();

    computed      = true;
    reconstructed = false;
    }

  const std::set<std::string>
  KoopmanRKHS::valid_keys = {
    "d",
    "alpha",
    "kepler", "plummer", "radial", "oscil", "testF", "lam", "mu", "c", "tscale",
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

  Eigen::VectorXcd KoopmanRKHS::modeEval(int index, const Eigen::VectorXd& x)
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

  std::complex<double> KoopmanRKHS::evecEval(int index, const Eigen::VectorXd& x)
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

  void KoopmanRKHS::assignParameters(const std::string flags)
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
	throw YamlConfigError("MSSA::KoopmanRKHS", "parameter", unmatched, __FILE__, __LINE__);

      // Compute flags
      //
      computed      = false;
      reconstructed = false;

      // Top level parameter flags
      //
      if (params["d"])
	d  = double(params["d"].as<double>());

      if (params["lam"])
	lam = double(params["lam"].as<double>());

      if (params["mu"])
	mu = double(params["mu"].as<double>());

      if (params["c"])
	c = double(params["c"].as<double>());

      if (params["kepler"])
	kepler = bool(params["kepler"].as<bool>());

      if (params["plummer"])
	plummer = bool(params["plummer"].as<bool>());

      if (params["radial"])
	radial = bool(params["radial"].as<bool>());

      if (params["tscale"])
	tscale = double(params["tscale"].as<double>());

      if (params["oscil"])
	oscil = bool(params["oscil"].as<bool>());

      if (params["testF"])
	testF = bool(params["testF"].as<bool>());

      if (params["alpha"])
	alpha = double(params["alpha"].as<double>());

      if (params["kernel"]) {
	std::string type = params["kernel"].as<std::string>();
	if (RKHS_values.find(type) != RKHS_values.end())
	  rkhs = RKHS_values[type];
	else
	  throw std::runtime_error("KoopmanRKHS: kernel type [" + type +
				   "] not found");
      }
      
      if (params["verbose"])
	verbose  = bool(params["verbose"]);

      if (params["output"]) prefix = params["output"].as<std::string>();
      else                  prefix = "exp_rkhs";

    }
    catch (const YAML::ParserException& e) {
      std::cout << "KoopmanRKHS::assignParameters, parsing error=" << e.what()
		<< std::endl;
      throw;
    }
  }


  // Save current KOOPMAN state to an HDF5 file with the given prefix
  void KoopmanRKHS::saveState(const std::string& prefix)
  {
    if (not computed) return;	// No point in saving anything

    if (std::filesystem::exists(prefix + "_rkhs.h5")) {
      std::ostringstream sout;
      sout << "KoopmanRKHS::saveState: the file <" <<  prefix + "_edmd.h5"
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
  void KoopmanRKHS::restoreState(const std::string& prefix)
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
      int nTime, nKeys, nEV;

      h5file.getAttribute("numT" ).read(nTime);
      h5file.getAttribute("nKeys").read(nKeys);
      h5file.getAttribute("nEV"  ).read(nEV  );

      // Number of channels
      //
      nkeys = data.size();

      // Test recovered parameters
      //
      if (nTime != numT) {
	std::ostringstream sout;
	sout << "KoopmanRKHS::restoreState: saved state has numT="
	     << nTime << " but KoopmanRKHS expects numT=" << numT
	     << ".\nCan't restore KoopmanRKHS state!";
	throw std::runtime_error(sout.str());
      }

      if (nKeys != nkeys) {
	std::ostringstream sout;
	sout << "KoopmanRKHS::restoreState: saved state has nkeys="
	     << nKeys << " but KoopmanRKHS expects nkeys=" << nkeys
	     << ".\nCan't restore KoopmanRKHS state!";
	throw std::runtime_error(sout.str());
      }

      if (nEV != nev) {
	std::ostringstream sout;
	sout << "KoopmanRKHS::restoreState: saved state has nEV="
	     << nEV << " but KoopmanRKHS expects nEV=" << nev
	     << ".\nCan't restore KoopmanRKHS state!";
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
	sout << "KoopmanRKHS::restoreState: keylist mismatch." << std::endl
	     << "Can't restore KoopmanRKHS state! Wanted keylist: ";
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

  KoopmanRKHS::KoopmanRKHS(const mssaConfig& config, double tol, const std::string flags) : tol(tol)
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
  // END KoopmanRKHS constructor

}
// END namespace MSSA
