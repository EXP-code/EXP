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
  double KoopmanRKHS::kernel(const Eigen::VectorXd& x, const Eigen::VectorXd& y)
  {
    if (rkhs == RKHS::Polynomial) {
      double prod = x.adjoint()*y;
      return pow(1.0 + prod/(d*d), alpha);
    } else if (rkhs == RKHS::Exponential) {
      double prod = x.adjoint()*y;
      return exp(mu*prod);
    } else if (rkhs == RKHS::Gaussian) {
      auto diff = x - y;
      double prod = diff.adjoint()*diff;
      return exp(-prod/(sigma*sigma));
    } else {
      throw std::runtime_error("KoopmanRKHS: Unknown kernel type");
    }
  }

  // Algorithm 3 from Williams, Rowley, Kevrekidis
  //
  void KoopmanRKHS::koopman_analysis()
  {
    // Number of channels
    //
    auto keys = getAllKeys();
    nkeys = keys.size();

    numT = data[keys[0]].size();

    G.resize(numT, numT);
    A.resize(numT, numT);
    X.resize(numT, nkeys);
    
    for (int k=0; k<nkeys; k++) {
      for (int i=0; i<numT; i++) X(i, k) = data[keys[k]][i];
    }
      

    Eigen::VectorXd Y(nkeys);

    double DT = coefDB.times[1] - coefDB.times[0];

    // Compute Grammian and action matrices
    for (int i=0; i<numT; i++) {
      for (int j=0; j<numT; j++) {
	
	for (int k=0; k<nkeys; k++) {
	  if (j==0) {
	    Y(k) = (data[keys[k]][1] - data[keys[k]][0])/DT;
	  } else if (j==numT-1) {
	    Y(k) = (data[keys[k]][numT-1] - data[keys[k]][numT-2])/DT;
	  } else {
	    Y(k) = (data[keys[k]][j+1] - data[keys[k]][j-1])/(2.0*DT);
	  }
	}
	
	G(i, j) = kernel(X.row(i), X.row(j));
	A(i, j) = kernel(X.row(i), Y);
      }
    }

    // Perform eigenanalysis of Grammian matrix (pos def)
    //
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(G);
    if (eigensolver.info() != Eigen::Success) {
      throw std::runtime_error("KoopmanRKHS: Eigen solver failed");
    }
    S2 = eigensolver.eigenvalues();
    Q  = eigensolver.eigenvectors();
    
    // Find inversion tolerance condition
    //
    SP.resize(S2.size());
    double TOL = S2(S2.size()-1) * tol;
    int nev = 0;
    for (int n=0; n<S2.size(); n++) {
      if (S2(n) < TOL) {
	SP(n) = 0.0;
      } else {
	SP(n) = 1.0 / S2(n);
	nev++;
      }
    }

    // Get Koopman operator estimate
    //
    K = (SP.asDiagonal() * Q.transpose()) * A * (Q * SP.asDiagonal());

    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver2(K), eigensolver3(K.transpose());
    SR = eigensolver2.eigenvalues();
    V  = eigensolver2.eigenvectors();
    SL = eigensolver3.eigenvalues();
    U  = eigensolver3.eigenvectors();

    Xi = (U.transpose()*SP.asDiagonal()*Q.transpose()*X).transpose();

    computed = true;
    reconstructed = false;
  }

  const std::set<std::string>
  KoopmanRKHS::valid_keys = {
    "d",
    "alpha",
    "mu",
    "sigma",
    "verbose",
    "output"
  };

  void KoopmanRKHS::assignParameters(const std::string flags)
  {
    // Default parameter values
    //
    d        = 1.0;
    alpha    = 10.0;
    mu       = 1.0;
    sigma    = 1.0;
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
      d        = double(params["d"      ].as<double>());
      alpha    = double(params["alpha"  ].as<double>());
      mu       = double(params["mu"     ].as<double>());
      sigma    = double(params["sigma"  ].as<double>());
      verbose  = bool  (params["verbose"]);

      if (params["output"] ) prefix = params["output"].as<std::string>();
      else                   prefix = "exp_rkhs";

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

  KoopmanRKHS::KoopmanRKHS(const mssaConfig& config, const std::string flags)
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
