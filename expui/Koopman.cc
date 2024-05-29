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
// Jonathan H. Tu, Clarence W. Rowley, Dirk M. Luchtenburg, Steven
// L. Brunton and J. Nathan Kutz, 2014, "On Dynamic Mode
// Decomposition: Theory and Applications", Journal of Computational
// Dynamics, 1(2), 391-421
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

#include <KMeans.H>

#include <Koopman.H>

#include <RedSVD.H>
#include <YamlConfig.H>
#include <YamlCheck.H>
#include <EXPException.H>
#include <libvars.H>

#include <TransformFFT.H>

#ifdef HAVE_LIBPNGPP
#include <ColorGradient.H>
#endif

#include "Koopman.H"

namespace MSSA {

  // Helper ostream manipulator for debug coefficient key info
  std::ostream& operator<< (std::ostream& out, const std::vector<unsigned>& t);

  // Do the SVD, populate singular value vectors
  //
  void Koopman::koopman_analysis()
  {
    // Number of channels
    //
    nkeys = data.size();

    // Enforce nev to be <= rank
    //
    if (nev > nkeys) std::cout << "Koopman: setting nEV=" << nkeys << std::endl;
    nev = std::min<int>(nev, nkeys);

    // Allocate trajectory arrays
    //
    X0.resize(nkeys, numT-1);
    X1.resize(nkeys, numT-1);

    X0.fill(0.0);
    X1.fill(0.0);

    // Build input & output data series
    //
    int n=0;
    for (auto k : data) {
      for (int j=0; j<numT-1; j++) {
	X0(n, j) = data[k.first][j+0];
	X1(n, j) = data[k.first][j+1];
      }
      n++;
    }

    // Approximate the pseudoinverse using the rank-r SVD
    // approximation of the initial state matrix
    //
    // Use one of the built-in Eigen3 algorithms
    //
    if (params["Jacobi"]) {
      // -->Using Jacobi
      Eigen::JacobiSVD<Eigen::MatrixXd>
	svd(X0, Eigen::ComputeThinU | Eigen::ComputeThinV);
      S = svd.singularValues();
      U = svd.matrixU();
      V = svd.matrixV();
    } else if (params["BDCSVD"]) {
      // -->Using BDC
      Eigen::BDCSVD<Eigen::MatrixXd>
	svd(X0, Eigen::ComputeThinU | Eigen::ComputeThinV);
      S = svd.singularValues();
      U = svd.matrixU();
      V = svd.matrixV();
    } else {
      // -->Use Random approximation algorithm from Halko, Martinsson,
      //    and Tropp
      RedSVD::RedSVD<Eigen::MatrixXd> svd(X0, nev);
      S = svd.singularValues();
      U = svd.matrixU();
      V = svd.matrixV();
    }

    // Compute the approximation to the Koopman operator for the rank
    // reduced approximation
    //
    Eigen::MatrixXd D = S.asDiagonal();

    // E.g. Tu et al. 2014, equation 4 (parens to enforce effficient
    // order)
    //
    A = U.transpose() * (X1 * V) * D.inverse();

    // Now compute the eigenvalues and eigenvectors
    //
    Eigen::EigenSolver<Eigen::MatrixXd> sol(A, true);

    L = sol.eigenvalues();
    W = sol.eigenvectors();

    // Compute the EDMD modes
    //
    Eigen::VectorXcd Linv(L);
    for (int i=0; i<L.size(); i++) {
      if (Linv(i) != 0.0) Linv(i) = 1.0/Linv(i);
      else Linv(i) = 0.0;
    }

    // Projected mode for testing
    //
    if (project) {
      Phi = U * W;
    }
    // This is the exact mode from Tu et al. 2014, equation 9
    //
    else {
      Phi = Linv.asDiagonal() * X1 * V * D.inverse() * W;
    }
    
    computed = true;
    reconstructed = false;
  }

  void Koopman::reconstruct(const std::vector<int>& evlist)
  {
    // Prevent a belly-up situation
    //
    if (not computed) koopman_analysis();

    // Make a zero vector
    //
    Eigen::VectorXcd I(nev); I.setZero();

    auto lsz = evlist.size();

    Y.resize(numT, nkeys);
    Y.setZero();

    if (lsz) {

      int n = 0;
      Eigen::VectorXd xx(nkeys);
      xx.setZero();
      for (auto u : data) xx[n++] = data[u.first][0];
	
      for (auto v : evlist) {
	if (v<nev) I[v] = 1.0;
      }

      Eigen::VectorXcd B  = Phi.inverse() * xx;
      Eigen::MatrixXcd LL = I.asDiagonal();
	
      // Propate the solution with the operator
      //
      for (int i=0; i<numT; i++) {
	Y.row(i) = (Phi*LL*B).real();
	LL *= L.asDiagonal();
      }
    }

    reconstructed = true;
  }

  // This computes an image of the contributions
  //
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> Koopman::contributions()
  {
    Eigen::MatrixXd retF, retG;

    if (not reconstructed) {

      if (not computed)
	std::cout << "Koopman::contributions: "
		  << "call eigenvalues() or getPC() followed by "
		  << "reconstruct() before contributions()"
		  << std::endl;
      else
	std::cout << "Koopman::contributions: "
		  << "call reconstruct() before contributions()"
		  << std::endl;

      return std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>(retF, retG);
    }

    retF.resize(nev, nkeys);
    retG.resize(nev, nkeys);

    retF.setZero();

    int nn = 0;
    Eigen::VectorXd xx(nkeys);
    for (auto u : data) xx[nn++] = data[u.first][0];
    
    Eigen::VectorXcd B = Phi.inverse() * xx;
    Eigen::VectorXcd LL = Eigen::VectorXd::Ones(L.size());

    for (int i=0; i<numT; i++) {
      for (int j=0; j<nev; j++) {
	for (int n=0; n<nkeys; n++) {
	  retF(j, n) += std::norm( Phi(j, n)*LL(n)*B(n) );
	}
      }
      LL = LL.array() * L.array();
    }
    retF /= numT;
    retG = retF;

    // This is norm for each series over the entire reconstruction
    // (that is, fixed coefficient channel, summed over all PCs)

    Eigen::VectorXd norm(nev);
    norm.setZero();

    for (int j=0; j<nev; j++) {
      for (int k=0; k<nkeys; k++) {
	norm[j] += retF(j, k);
      }
    }

    for (int j=0; j<nev; j++) {
      for (int k=0; k<nkeys; k++) {
	if (norm[j]>0.0) retF(j, k) /= norm[j];
	retF(j, k) = sqrt(retF(j, k));
      }
    }

    norm.resize(nkeys);
    norm.setZero();

    for (int j=0; j<nev; j++) {
      for (int k=0; k<nkeys; k++) {
	norm[k] += retG(j, k);
      }
    }

    for (int j=0; j<nev; j++) {
      for (int k=0; k<nkeys; k++) {
	if (norm[k]>0.0) retG(j, k) /= norm[k];
	retG(j, k) = sqrt(retG(j, k));
      }
    }

    return std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>(retF, retG);
  }

  // This computes an image of the contributions
  //
  void Koopman::contributionsPNG()
  {
    if (not computed) {
      std::cout << "Koopman::contributions: "
		<< "call eigenvalues() or getModes() before contributions()"
		<< std::endl;
      return;
    }

#ifdef HAVE_LIBPNGPP

    auto ret  = contributions();
    auto retF = std::get<0>(ret);
    auto retG = std::get<1>(ret);

    std::string filename = prefix + ".f_contrib";
    std::ofstream out(filename);

    // Used in both computations
    //
    std::map<Key, std::vector<double>, mSSAkeyCompare> values;

    //
    const int minSize = 600;
    int ndupX = 1, ndupY = 1;
    if (nev < minSize) ndupX = minSize/nev + 1;
    if (data.size() < minSize) ndupY = minSize/data.size() + 1;

    std::cout << "Matrix size: " << nev << " X " << data.size() << std::endl;

    png::image< png::rgb_pixel > image1(nev*ndupX, data.size()*ndupY);
    png::image< png::rgb_pixel > image2(nev*ndupX, data.size()*ndupY);
    ColorGradient color;
    color.createFiveColorHeatMapGradient();

    if (out) {

      // Column header
      //
      out << std::setw(4) << "# m"
	  << std::setw(4) << "n"
	  << std::setw(4) << "cs";
      for (int j=0; j<nev; j++) {
	std::ostringstream sout; sout << "EV" << j;
	out << std::setw(18) << sout.str();
      }
      out << std::endl;

      int n = 0;
      double maxVal = 0.0;
      for (auto u : data) {
	auto key = u.first;
	for (auto q : key) out << std::setw(4) << q;

	for (int j=0; j<nev; j++) {
	  out << std::setw(18) << sqrt(retF(n, j));
	  maxVal = std::max<double>(maxVal, retF(n, j));
	}
	out << std::endl;
	n++;
      }
      out.close();

      int i = 0;
      for (auto u : data) {
	auto key = u.first;
	for (int j=0; j<nev; j++) {
	  png::rgb_pixel cval = color(sqrt(retF(i, j)/maxVal));
	  for (size_t yy = i*ndupY; yy < (i+1)*ndupY; yy++) {
	    for (size_t xx = j*ndupX; xx < (j+1)*ndupX; xx++) {
	      image1[yy][xx] = cval;
	    }
	  }
	}
	i++;
      }
      out.close();

      image1.write(filename + ".png");

    } else {
      std::cout << "Could not open <" << filename << ">" << std::endl;
      exit(-1);
    }

    filename = prefix + ".g_contrib";
    out.open(filename);
    if (out) {
      // Column header
      //
      out << std::setw(8) << "# coef";
      for (auto u : data) {
	std::ostringstream sout;
	sout << u.first;
	out << std::setw(18) << sout.str();
      }
      out << std::endl;

      std::vector<double> norm(nev, 0.0);

      double maxVal = 0.0;
      for (int j=0; j<nev; j++) {
	out << std::setw(8) << j;
	int i = 0;
	for (auto u : data) {
	  auto key = u.first;
	  out << std::setw(18) << sqrt(retG(i, j));
	  maxVal = std::max<double>(maxVal, retG(i, j));
	}
	out << std::endl;
      }
      out.close();

      for (int j=0; j<nev; j++) {
	int i = 0;
	for (auto u : data) {
	  auto key = u.first;
	  png::rgb_pixel cval = color(sqrt(retG(i, j)/maxVal));
	  for (size_t yy = i*ndupY; yy < (i+1)*ndupY; yy++) {
	    for (size_t xx = j*ndupX; xx < (j+1)*ndupX; xx++) {
	      image2[yy][xx] = cval;
	    }
	  }
	  i++;
	}
      }

      image2.write(filename + ".png");

    } else {
      std::cout << "Could not open <" << filename << ">" << std::endl;
      exit(-1);
    }
#else
    std::cout << "PNG is not available, so I can't make contribution images"
	      << std::endl;
#endif

  }

  // Return the DFT of the individual reconstructed data channels
  //
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
  Koopman::channelDFT()
  {
    Eigen::VectorXd F, P, fw;
    Eigen::MatrixXd pw;

    if (not reconstructed) {

      if (not computed)
	std::cout << "Koopman::channelDFT: "
		  << "call eigenvalues() or getModes() followed by "
		  << "reconstruct() before channelDFT()"
		  << std::endl;
      else
	std::cout << "Koopman::channelDFT: "
		  << "call reconstruct() before channelDFT()"
		  << std::endl;

      return {fw, pw};
    }

    {
      double DT = coefDB.times[1] - coefDB.times[0];

      int nfreq = numT/2 + 1;
      int nchan = data.size();

      fw.resize(nfreq);
      pw.resize(nfreq, nchan);

      Eigen::VectorXd p0(nfreq), in(numT);

      int nch = 0;
      for (auto u : data) {

	// Compute the summed data stream
	//
	for (int i=0; i<numT; i++) in(i) = Y(i, nch);

	// Compute the DFT for the summed data stream
	//
	TransformFFT fft(DT, in);
	fft.Power(F, p0);

	// Pack the results for return
	//
	for (int k=0; k<nfreq; k++) {
	  if (nch==0) fw(k) = F(k); // Only need to do this once
	  pw(k, nch) = p0(k);	    // Get the data for each channel
	}

	if (powerf) {

	  // Only need to do the partial DFTs if files are requested
	  //
	  for (int i=0; i<numT; i++) in(i) = Y(i, nch);

	  TransformFFT fft(DT, in);
	  fft.Power(F, P);

	  std::ostringstream filename;
	  filename << prefix << ".power_" << u.first;
	  std::ofstream out(filename.str());
	  if (out) {
	    out << "# " << u.first << std::endl;
	    out << "# " << std::setw(13) << "Freq"
		<< std::setw(15) << "Period"
		<< std::setw(15) << "Summed"
		<< std::setw(15) << "Full";
	    for (int j=0; j<nev; j++) {
	      std::ostringstream sout; sout << "Mode " << j;
	      out << std::setw(15) << sout.str();
	    }
	    out << "# " << std::setw(13) << "[1]"
		<< std::setw(15) << "[2]"
		<< std::setw(15) << "[3]";
	    for (int j=0; j<nev; j++) {
	      std::ostringstream sout; sout << '[' << j+4 << ']';
	      out << std::setw(15) << sout.str();
	    }
	    out << std::endl;

	    for (int j=0; j<nfreq; j++) {
	      out << std::setw(15) << std::setprecision(6) << F(j)
		  << std::setw(15) << std::setprecision(6) << 2.0*M_PI/F(j)
		  << std::setw(15) << std::setprecision(6) << p0(j)
		  << std::setw(15) << std::setprecision(6) << P(j)
		  << std::endl;
	    }
	    out.close();
	  } else {
	    std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	    exit(-1);
	  }
	}

	// Augment the channel counter
	//
	nch++;
      }
    }

    return {fw, pw};
  }

  std::map<std::string, CoefClasses::CoefsPtr> Koopman::getReconstructed()
  {
    if (not reconstructed) {
      std::ostringstream sout;
      if (not computed)
	sout << "Koopman::getReconstructed(): "
	     << "call eigenvalues() or getModes() followed by "
	     << "reconstruct() before getReconstructed()";
      else
	sout << "Koopman::getReconstructed(): "
	     << "call reconstruct() before getReconstructed()";

      std::runtime_error(sout.str());
    }


    // Copy the original map for return
    //
    auto newdata = data;

    for (int i=0; i<numT; i++) {
      int n = 0;
      for (auto u : newdata) {
  	newdata[u.first][i] = Y(i, n++);
      }
    }

    // Copy to the working data
    //
    for (auto v : newdata) {
      if (verbose) std::cout << "Updating for: " << v.first << std::endl;
      coefDB.setData(v.first, v.second);
    }

    // Copies working data back to the coefficient structures
    //
    return coefDB.endUpdate();
  }

  const std::set<std::string>
  Koopman::valid_keys = {
    "verbose",
    "power",
    "Jacobi",
    "BDCSVD",
    "project",
    "output"
  };

  void Koopman::assignParameters(const std::string flags)
  {
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
	throw YamlConfigError("MSSA::Koopman", "parameter", unmatched, __FILE__, __LINE__);

      // Compute flags
      //
      computed      = false;
      reconstructed = false;

      // Top level parameter flags
      //
      verbose  = bool(params["verbose"]);
      powerf   = bool(params["power"  ]);
      project  = bool(params["project"]);

      if (params["output"] ) prefix = params["output"].as<std::string>();
      else                   prefix = "exp_edmd";

    }
    catch (const YAML::ParserException& e) {
      std::cout << "Koopman::assignParameters, parsing error=" << e.what()
		<< std::endl;
      throw;
    }
  }


  // Save current KOOPMAN state to an HDF5 file with the given prefix
  void Koopman::saveState(const std::string& prefix)
  {
    if (not computed) return;	// No point in saving anything

    if (std::filesystem::exists(prefix + "_edmd.h5")) {
      std::ostringstream sout;
      sout << "Koopman::saveState: the file <" <<  prefix + "_edmd.h5"
	   << "> already exists.\nPlease delete this file or choose a "
	   << "different file name";
      throw std::runtime_error(sout.str());
    }

    try {
      // Create a new hdf5 file
      //
      HighFive::File file(prefix + "_edmd.h5",
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
      HighFive::Group analysis = file.createGroup("koopman_analysis");

      analysis.createDataSet("Phi",  Phi);
      analysis.createDataSet("X0",   X0 );
      analysis.createDataSet("X1",   X1 );
      analysis.createDataSet("U",    U  );
      analysis.createDataSet("V",    V  );
      analysis.createDataSet("A",    A  );
      analysis.createDataSet("L",    L  );
      analysis.createDataSet("W",    W  );
      analysis.createDataSet("Y",    Y  );

    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }
  }

  // Restore current KOOPMAN state to an HDF5 file with the given prefix
  void Koopman::restoreState(const std::string& prefix)
  {
    try {
      // Silence the HDF5 error stack
      //
      HighFive::SilenceHDF5 quiet;

      // Try opening the file as HDF5
      //
      HighFive::File h5file(prefix + "_edmd.h5", HighFive::File::ReadOnly);

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
	sout << "Koopman::restoreState: saved state has numT="
	     << nTime << " but Koopman expects numT=" << numT
	     << ".\nCan't restore Koopman state!";
	throw std::runtime_error(sout.str());
      }

      if (nKeys != nkeys) {
	std::ostringstream sout;
	sout << "Koopman::restoreState: saved state has nkeys="
	     << nKeys << " but Koopman expects nkeys=" << nkeys
	     << ".\nCan't restore Koopman state!";
	throw std::runtime_error(sout.str());
      }

      if (nEV != nev) {
	std::ostringstream sout;
	sout << "Koopman::restoreState: saved state has nEV="
	     << nEV << " but Koopman expects nEV=" << nev
	     << ".\nCan't restore Koopman state!";
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
	sout << "Koopman::restoreState: keylist mismatch." << std::endl
	     << "Can't restore Koopman state! Wanted keylist: ";
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

      auto analysis = h5file.getGroup("koopman_analysis");

      Phi = analysis.getDataSet("Phi").read<Eigen::MatrixXcd>();
      X0  = analysis.getDataSet("X0" ).read<Eigen::MatrixXd >();
      X1  = analysis.getDataSet("X1" ).read<Eigen::MatrixXd >();
      U   = analysis.getDataSet("U"  ).read<Eigen::MatrixXd >();
      V   = analysis.getDataSet("V"  ).read<Eigen::MatrixXd >();
      A   = analysis.getDataSet("A"  ).read<Eigen::MatrixXd >();
      L   = analysis.getDataSet("L"  ).read<Eigen::VectorXcd>();
      W   = analysis.getDataSet("W"  ).read<Eigen::MatrixXcd>();
      Y   = analysis.getDataSet("Y"  ).read<Eigen::MatrixXd >();

      computed = true;

    } catch (HighFive::Exception& err) {
      std::cerr << "**** Error opening or reading H5 file ****" << std::endl;
      throw;
    }

  }

  Koopman::Koopman(const mssaConfig& config, int nEV, const std::string flags) :nev(nEV)
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
  // END Koopman constructor

}
// END namespace MSSA
