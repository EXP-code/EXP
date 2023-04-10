//
// M-SSA for EXP coefficients
//
// Updated to use fixed-rank approximation for the SVD to save time
// and space.  Uses the approximate SVD randomized algorithm from
// Halko, Martinsson, and Tropp by default.  Use BDCSVD or Jacobi
// flags for the standard methods.  Using the randomized method
// together with trajectory matrix rather than covariance matrix
// analysis may lead to large errors in eigenvectors.
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

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

/* For debugging
   #undef eigen_assert
   #define eigen_assert(x)						\
   if (!(x)) { throw (std::runtime_error("Eigen assert error")); }
*/

#include <omp.h>

#include <KMeans.H>

#include <expMSSA.H>

#include <RedSVD.H>
#include <YamlConfig.H>
#include <YamlCheck.H>
#include <EXPException.H>
#include <libvars.H>

#include <TransformFFT.H>

#ifdef HAVE_LIBPNGPP
#include <ColorGradient.H>
#endif

#include "expMSSA.H"

namespace MSSA {

  Eigen::MatrixXd expMSSA::wCorrKey(const Key& key)
  {
    if (RC.find(key)==RC.end()) {
      throw std::runtime_error("expMSSA::wCorrKey: no such key");
    }

    auto R     = RC[key];

    int numT   = R.rows();
    int numW   = R.cols();
    int Lstar  = std::min<int>(numT - numW, numW);
    int Kstar  = std::max<int>(numT - numW, numW);

    // A Lambda for the weight function
    auto w = [&](int i) {
      if      (i < Lstar) return i;
      else if (i < Kstar) return Lstar;
      else                return numT - i + 1;
    };

    Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(numW, numW);
    for (int m=0; m<numW; m++) {
      for (int n=m; n<numW; n++) {
	for (int i=0; i<numT; i++) ret(m, n) += w(i) * R(i, m)*R(i, n);
      }
    }

    // Normalize
    //
    for (int m=0; m<numW; m++) {
      for (int n=m+1; n<numW; n++) {
	if (ret(m, m)>0.0 and ret(n, n)>0.0)
	  ret(m, n) /= sqrt(ret(m, m)*ret(n, n));
      }
    }

    // Unit diagonal
    //
    for (int m=0; m<numW; m++) ret(m, m) = 1.0;

    // Complete
    //
    for (int m=0; m<numW; m++) {
      for (int n=0; n<m; n++) ret(m, n) = ret(n, m);
    }

    return ret;
  }


  Eigen::MatrixXd expMSSA::wCorrAll()
  {
    int numT   = RC.begin()->second.rows();
    int numW   = RC.begin()->second.cols();
    int Lstar  = std::min<int>(numT - numW, numW);
    int Kstar  = std::max<int>(numT - numW, numW);

    // A Lambda for the weight function
    auto w = [&](int i) {
      if      (i < Lstar) return i;
      else if (i < Kstar) return Lstar;
      else                return numT - i + 1;
    };

    Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(numW, numW);
    for (auto R : RC) {
      for (int m=0; m<numW; m++) {
	for (int n=m; n<numW; n++) {
	  for (int i=0; i<numT; i++)
	    ret(m, n) += w(i) * R.second(i, m)*R.second(i, n);
	}
      }
    }

    // Normalize
    //
    for (int m=0; m<numW; m++) {
      for (int n=m+1; n<numW; n++) {
	if (ret(m, m)>0.0 and ret(n, n)>0.0)
	  ret(m, n) /= sqrt(ret(m, m)*ret(n, n));
      }
    }

    // Unit diagonal
    //
    for (int m=0; m<numW; m++) ret(m, m) = 1.0;

    // Complete
    //
    for (int m=0; m<numW; m++) {
      for (int n=0; n<m; n++) ret(m, n) = ret(n, m);
    }

    return ret;
  }

  Eigen::MatrixXd expMSSA::wCorr(const std::string& name, const Key& key)
  {
    int indx = coefDB.index(name);
    if (indx<0) {
      std::cout << "No such name <" << name << ">" << std::endl
		<< "Available names are:";
      for (auto v : coefDB.getNames()) std::cout << " " << v;
      std::cout << std::endl;
      throw std::runtime_error("expMSSA::wCorr: invalid component name");
    }

    if (coefDB.isComplex(name)) {
      auto ekey = key;
      int pos = ekey.size();

      ekey.push_back(0);
      ekey.push_back(indx);

      auto mat = wCorrKey(ekey);

      ekey[pos] = 1;
      if (RC.find(ekey)!=RC.end()) {
	mat += wCorrKey(ekey);
	mat *= 0.5;
      }
      return mat;
    }
    else {
      auto ekey = key;
      ekey.push_back(indx);
      return wCorrKey(ekey);
    }

  }

  // Helper ostream manipulator for debug coefficient key info
  std::ostream& operator<< (std::ostream& out, const std::vector<unsigned>& t)
  {
    const char lab[] = {'c', 's'};
    std::ostringstream sout;
    sout << '(';
    for (int i=0; i<t.size()-2; i++) sout << t[i] << ',';
    sout << lab[t[t.size()-2]] << ")_" << t.back();
    out << sout.str();
    return out;
  }


  // Do the SVD, populate singular value vectors
  //
  void expMSSA::mssa_analysis()
  {
    // Number of channels
    //
    nkeys = mean.size();

    // Make sure parameters are sane
    //
    if (numW<=0) numW = numT/2;
    if (numW > numT/2) numW = numT/2;

    numK = numT - numW + 1;

    Y.resize(numK, numW*nkeys);
    Y.fill(0.0);

    // Build embedded time series.  Augmented vectors are the rows.
    //
    {
      size_t n=0;
      for (auto k : mean) {
	for (int j=0; j<numW; j++) {
	  for (int i=0; i<numK; i++)
	    if (i + j < numT) Y(i, numW*n+j) = data[k.first][i + j];
	}
	n++;
      }
    }

    Eigen::MatrixXd cov;
    double Scale;
    int rank;

    // Covariance is the default
    //
    if (params["Traj"]) {

      // Deduce the maximum rank of the trajectory matrix Y
      //
      rank = std::min<int>(
			   {static_cast<int>(Y.cols()),
			    static_cast<int>(Y.rows()), npc}
			   );
      Scale = Y.norm();

      if (Scale<=0.0) {
	std::cout << "Frobenius norm of trajectory matrix is <= 0!" << std::endl;
	exit(-1);
      }
    } else {
      cov   = Y.transpose() * Y/numK;
      rank  = std::min<int>(cov.cols(), npc);
      Scale = cov.norm();
    }

    if (Scale<=0.0) {
      std::cout << "Frobenius norm of trajectory or covariance matrix is <= 0!" << std::endl;
      exit(-1);
    } else {
      cov /= Scale;
    }

    // Only write covariance matrix on request
    //
    if (not params["Traj"] and params["writeCov"]) {
      std::string filename = prefix + ".cov";
      std::ofstream out(filename);
      out << cov;
      out.close();
    }

    // Use one of the built-in Eigen3 algorithms
    //
    if (params["Jacobi"]) {
      // -->Using Jacobi
      if (params["Traj"]) {	// Trajectory matrix
	auto YY = Y/Scale;
	Eigen::JacobiSVD<Eigen::MatrixXd>
	  svd(YY, Eigen::ComputeThinU | Eigen::ComputeThinV);
	S = svd.singularValues();
	U = svd.matrixV();
      }
      else {			// Covariance matrix
	Eigen::JacobiSVD<Eigen::MatrixXd>
	  svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
	S = svd.singularValues();
	U = svd.matrixU();
      }
    } else if (params["BDCSVD"]) {
      // -->Using BDC
      if (params["Traj"]) {	// Trajectory matrix
	auto YY = Y/Scale;
	Eigen::BDCSVD<Eigen::MatrixXd>
	  svd(YY, Eigen::ComputeThinU | Eigen::ComputeThinV);
	S = svd.singularValues();
	U = svd.matrixV();
      }
      else {			// Covariance matrix
	Eigen::BDCSVD<Eigen::MatrixXd>
	  svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
	S = svd.singularValues();
	U = svd.matrixU();
      }
    } else {
      // -->Use Random approximation algorithm from Halko, Martinsson,
      //    and Tropp
      if (params["Traj"]) {	// Trajectory matrix
	auto YY = Y/Scale;
	RedSVD::RedSVD<Eigen::MatrixXd> svd(YY, rank);
	S = svd.singularValues();
	U = svd.matrixV();
      }
      else {			// Covariance matrix
	if (params["RedSym"]) {
	  RedSVD::RedSymEigen<Eigen::MatrixXd> eigen(cov, rank);
	  S = eigen.eigenvalues().reverse();
	  U = eigen.eigenvectors().rowwise().reverse();
	} else {
	  RedSVD::RedSVD<Eigen::MatrixXd> svd(cov, rank);
	  S = svd.singularValues();
	  U = svd.matrixU();
	}
      }
    }

    std::cout << "shape U = " << U.rows() << " x "
	      << U.cols() << std::endl;

    // Rescale the SVD factorization by the Frobenius norm
    //
    S *= Scale;
    if (params["Traj"]) {
      for (int i=0; i<S.size(); i++) S(i) = S(i)*S(i)/numK;
    }

    double tot = 0.0;
    for (int j=0; j<S.size(); j++) tot += S(j);

    npc = std::min<int>(npc, numW*nkeys);

    // Compute the PCs by projecting the data
    //
    PC = Y * U;

    computed = true;
    reconstructed = false;
  }

  void expMSSA::reconstruct(const std::vector<int>& evlist)
  {
    // Prevent a belly-up situation
    //
    if (not computed) mssa_analysis();

    // Use the OpenMP implementation
    const bool useOpenMP  = true;

    // Some speed up by reordering memory access at the expense of
    // some memory
    const bool useReverse = false;

    // Reconstructed time series
    //
    ncomp = std::min<int>({numW, npc, static_cast<int>(PC.cols())});

    // Make a bool vector
    //
    std::vector<bool> I(ncomp, false);

    auto lsz = evlist.size();

    if (lsz) {
      for (auto v : evlist) if (v<ncomp) I[v] = true;
    }

    for (auto u : mean) {
      RC[u.first].resize(numT, ncomp);
      RC[u.first].setZero();
    }

    if (lsz) {

      // Embedded time series matrix
      //
      Eigen::MatrixXd Z;
      if (useReverse) Z.resize(numT, numW);

      // Split eigenvectors into channel sequences
      //
      std::map<Key, Eigen::MatrixXd, mSSAkeyCompare> rho;

      int n = 0;
      for (auto u : mean) {
	rho[u.first].resize(numW, ncomp);
	rho[u.first].fill(0.0);
	for (int i=0; i<numW; i++) {
	  for (int j=0; j<ncomp; j++)
	    if (I[j]) rho[u.first](i, j) = U(numW*n + i, j);
	}
	n++;
      }

#pragma omp parallel
      if (useOpenMP) {
	// Parallelize the map iteration by wrapping it in a standard loop.
	// Ungainly for sure but it works.
	//
	int thread_count = omp_get_num_threads();
	int thread_num   = omp_get_thread_num();
	size_t map_size  = mean.size();

	auto u = mean.begin();
	std::advance(u, thread_num);

	for (int q=thread_num; q<map_size; q+=thread_count) {

	  for (int w=0; w<ncomp; w++) {

	    // Choose limits and weight
	    //
	    for (int i=0; i<numT; i++) {
	      double W;
	      int L, U;

	      if (i<numW) {		// Lower
		W = 1.0/(1.0 + i);
		L = 0;
		U = i + 1;
	      } else if (i<numT-numW) { // Middle
		W = 1.0/numW;
		L = 0;
		U = numW;
	      } else {		// Upper
		W = 1.0/(numT - i);
		L = i - numT + numW;
		U = numW;
	      }

	      RC[u->first](i, w) = 0.0;
	      for (int j=L; j<U; j++) {
		RC[u->first](i, w) += PC(i - j, w) * rho[u->first](j, w) * W;
	      }
	    }
	  }

	  if (q+thread_count < map_size) std::advance(u, thread_count);
	}
      }
      // The original serial implementation
      //
      else {

	for (auto u : mean) {

	  for (int w=0; w<ncomp; w++) {

	    // Build reverse embedded time series.
	    //
	    if (useReverse) {
	      Z.fill(0.0);
	      for (int j=0; j<numW; j++) {
		for (int i=0; i<numT; i++)
		  if (i - j >= 0 and i - j < numK) Z(i, j) = PC(i - j, w);
	      }
	    }

	    // Choose limits and weight
	    //
	    for (int i=0; i<numT; i++) {
	      double W;
	      int L, U;

	      if (i<numW) {	// Lower
		W = 1.0/(1.0 + i);
		L = 0;
		U = i + 1;
	      } else if (i<numT-numW) { // Middle
		W = 1.0/numW;
		L = 0;
		U = numW;
	      } else {		// Upper
		W = 1.0/(numT - i);
		L = i - numT + numW;
		U = numW;
	      }

	      RC[u.first](i, w) = 0.0;
	      for (int j=L; j<U; j++) {
		if (useReverse)
		  RC[u.first](i, w) += Z(i, j) * rho[u.first](j, w) * W;
		else
		  RC[u.first](i, w) += PC(i - j, w) * rho[u.first](j, w) * W;
	      }
	    }
	  }
	}
      }
    }

    reconstructed = true;
  }

  // This computes an image of the contributions
  //
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> expMSSA::contributions()
  {
    Eigen::MatrixXd retF, retG;

    if (not reconstructed) {

      if (not computed)
	std::cout << "expMSSA::contributions: "
		  << "call eigenvalues() or getPC() followed by "
		  << "reconstruct() before contributions()"
		  << std::endl;
      else
	std::cout << "expMSSA::contributions: "
		  << "call reconstruct() before contributions()"
		  << std::endl;

      return std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>(retF, retG);
    }

    retF.resize(ncomp, mean.size());
    retG.resize(ncomp, mean.size());

    int nk = 0;
    retF.setZero();
    for (auto u : mean) {
      auto key = u.first;
      for (int j=0; j<ncomp; j++) {
	for (int i=0; i<numT; i++) {
	  retF(j, nk) += RC[u.first](i, j)*RC[u.first](i, j);
	}
      }
      nk++;
    }

    retG = retF;

    // This is norm for each series over the entire reconstruction
    // (that is, fixed coefficient channel, summed over all PCs)

    Eigen::VectorXd norm(nk);
    norm.setZero();

    for (int n=0; n<nk; n++) {
      for (int j=0; j<ncomp; j++) {
	norm[n] += retF(j, n);
      }
    }

    for (int n=0; n<nk; n++) {
      for (int j=0; j<ncomp; j++) {
	if (norm[n]>0.0) retF(j, n) /= norm[n];
	retF(j, n) = sqrt(retF(j, n));
      }
    }

    norm.resize(ncomp);
    norm.setZero();

    for (int n=0; n<nk; n++) {
      for (int j=0; j<ncomp; j++) {
	norm[j] += retG(j, n);
      }
    }

    for (int n=0; n<nk; n++) {
      for (int j=0; j<ncomp; j++) {
	if (norm[j]>0.0) retG(j, n) /= norm[j];
	retG(j, n) = sqrt(retG(j, n));
      }
    }


    return std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>(retF, retG);
  }

  // This computes an image of the contributions
  //
  void expMSSA::contributionsPNG()
  {
    if (not computed) {
      std::cout << "expMSSA::contributions: "
		<< "call eigenvalues() or getPC() before contributions()"
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
    if (npc < minSize) ndupX = minSize/npc + 1;
    if (mean.size() < minSize) ndupY = minSize/mean.size() + 1;

    std::cout << "Matrix size: " << npc << " X " << mean.size() << std::endl;

    png::image< png::rgb_pixel > image1(npc*ndupX, mean.size()*ndupY);
    png::image< png::rgb_pixel > image2(npc*ndupX, mean.size()*ndupY);
    ColorGradient color;
    color.createFiveColorHeatMapGradient();

    if (out) {

      // Column header
      //
      out << std::setw(4) << "# m"
	  << std::setw(4) << "n"
	  << std::setw(4) << "cs";
      for (int j=0; j<numW; j++) {
	std::ostringstream sout; sout << "PC" << j;
	out << std::setw(18) << sout.str();
      }
      out << std::endl;

      int n = 0;
      double maxVal = 0.0;
      for (auto u : mean) {
	auto key = u.first;
	for (auto q : key) out << std::setw(4) << q;

	for (int j=0; j<npc; j++) {
	  out << std::setw(18) << sqrt(retF(n, j));
	  maxVal = std::max<double>(maxVal, retF(n, j));
	}
	out << std::endl;
	n++;
      }
      out.close();

      int i = 0;
      for (auto u : mean) {
	auto key = u.first;
	for (int j=0; j<npc; j++) {
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
      for (auto u : mean) {
	std::ostringstream sout;
	sout << u.first;
	out << std::setw(18) << sout.str();
      }
      out << std::endl;

      std::vector<double> norm(npc, 0.0);

      double maxVal = 0.0;
      for (int j=0; j<npc; j++) {
	out << std::setw(8) << j;
	int i = 0;
	for (auto u : mean) {
	  auto key = u.first;
	  out << std::setw(18) << sqrt(retG(i, j));
	  maxVal = std::max<double>(maxVal, retG(i, j));
	}
	out << std::endl;
      }
      out.close();

      for (int j=0; j<npc; j++) {
	int i = 0;
	for (auto u : mean) {
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


  // Return the DFT of the PCs
  //
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
  expMSSA::pcDFT()
  {
    Eigen::VectorXd F, P, fw;
    Eigen::MatrixXd pw;
    {
      double DT = coefDB.times[1] - coefDB.times[0];

      int nfreq = numK/2 + 1;

      fw.resize(nfreq);
      pw.resize(nfreq, npc);
      Eigen::VectorXd in(numK);

      pw.setZero();

      for (int j=0; j<npc; j++) {
	for (int i=0; i<numK; i++) in(i) = PC(i, j);
	TransformFFT fft(DT, in);
	fft.Power(F, P);
	if (j==0) for (int i=0; i<nfreq; i++) fw(i) = F(i);
	for (int i=0; i<nfreq; i++) pw(i, j) = P(i);
      }

      if (powerf) {
	std::ostringstream filename;
	filename << prefix << ".pc_power";
	std::ofstream out(filename.str());
	if (out) {
	  out << "# "
	      << std::setw(13) << "Freq"
	      << std::setw(15) << "Period";
	  for (int j=0; j<npc; j++) {
	    std::ostringstream sout; sout << "PC " << j;
	    out << std::setw(15) << sout.str();
	  }
	  out << "# "
	      << std::setw(13) << "[1]"
	      << std::setw(15) << "[2]";
	  for (int j=0; j<npc; j++) {
	    std::ostringstream sout; sout << '[' << j+3 << ']';
	    out << std::setw(15) << sout.str();
	  }
	  out << std::endl;

	  for (int j=0; j<nfreq; j++) {
	    out << std::setw(15) << std::setprecision(6) << F(j)
		<< std::setw(15) << std::setprecision(6) << 2.0*M_PI/F(j);
	    for (int k=0; k<npc; k++)
	      out << std::setw(15) << std::setprecision(6) << pw(j, k);
	    out << std::endl;
	  }
	  out.close();
	} else {
	  std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	  exit(-1);
	}
      }
    }

    return {fw, pw};
  }

  // Return the DFT of the individual reconstructed data channels
  //
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
  expMSSA::channelDFT()
  {
    Eigen::VectorXd F, P, fw;
    Eigen::MatrixXd pw;

    if (not reconstructed) {

      if (not computed)
	std::cout << "expMSSA::channelDFT: "
		  << "call eigenvalues() or getPC() followed by "
		  << "reconstruct() before channelDFT()"
		  << std::endl;
      else
	std::cout << "expMSSA::channelDFT: "
		  << "call reconstruct() before channelDFT()"
		  << std::endl;

      return {fw, pw};
    }

    {
      double DT = coefDB.times[1] - coefDB.times[0];

      int nfreq = numT/2 + 1;
      int nchan = mean.size();

      fw.resize(nfreq);
      pw.resize(nfreq, nchan);

      Eigen::VectorXd p0(nfreq), in(numT);
      Eigen::MatrixXd pt(nfreq, ncomp);

      int nch = 0;
      for (auto u : mean) {

	// Compute the summed data stream
	//
	for (int i=0; i<numT; i++) {
	  in(i) = 0.0;
	  for (int j=0; j<ncomp; j++) in(i) += RC[u.first](i, j);
	}

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

	// Augment the channel counter
	//
	nch++;

	if (powerf) {

	  // Only need to do the partial DFTs if files are requested
	  //
	  for (int j=0; j<ncomp; j++) {
	    for (int i=0; i<numT; i++) in(i) = RC[u.first](i, j);

	    TransformFFT fft(DT, in);
	    fft.Power(F, P);
	    for (int k=0; k<nfreq; k++) pt(k, j) = P(k);
	  }

	  std::ostringstream filename;
	  filename << prefix << ".power_" << u.first;
	  std::ofstream out(filename.str());
	  if (out) {
	    out << "# " << u.first << std::endl;
	    out << "# " << std::setw(13) << "Freq"
		<< std::setw(15) << "Period"
		<< std::setw(15) << "Summed"
		<< std::setw(15) << "Full";
	    for (int j=0; j<ncomp; j++) {
	      std::ostringstream sout; sout << "PC " << j;
	      out << std::setw(15) << sout.str();
	    }
	    out << "# " << std::setw(13) << "[1]"
		<< std::setw(15) << "[2]"
		<< std::setw(15) << "[3]";
	    for (int j=0; j<ncomp; j++) {
	      std::ostringstream sout; sout << '[' << j+4 << ']';
	      out << std::setw(15) << sout.str();
	    }
	    out << std::endl;

	    for (int j=0; j<nfreq; j++) {
	      out << std::setw(15) << std::setprecision(6) << F(j)
		  << std::setw(15) << std::setprecision(6) << 2.0*M_PI/F(j)
		  << std::setw(15) << std::setprecision(6) << p0(j);
	      for (int k=0; k<ncomp; k++)
		out << std::setw(15) << std::setprecision(6) << pt(j, k);
	      out << std::endl;
	    }
	    out.close();
	  } else {
	    std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	    exit(-1);
	  }
	}
      }
    }

    return {fw, pw};
  }

  // Return the DFT for a single channel for each PC
  //
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
  expMSSA::singleDFT(const Key& key)
  {
    Eigen::VectorXd F, P, fw;
    Eigen::MatrixXd pt;
    {
      double DT = coefDB.times[1] - coefDB.times[0];

      int nfreq = numT/2 + 1;

      fw.resize(nfreq);
      pt.resize(nfreq, ncomp);

      Eigen::VectorXd in(numT);

      auto u = mean.find(key);

      // Make sure that we have the requested key
      //
      if (u == mean.end())
	throw std::runtime_error("expMSSA::singleDFT: requested key not found");

      for (int j=0; j<ncomp; j++) {

	// Pack the data for each PC
	//
	for (int i=0; i<numT; i++) in(i) = RC[u->first](i, j);

	// Perform the DFT
	//
	TransformFFT fft(DT, in);
	fft.Power(F, P);

	// Pack the results for return
	//
	for (int k=0; k<nfreq; k++) {
	  if (j==0) fw(k) = F(k);
	  pt(k, j) = P(k);
	}
      }
    }

    return {fw, pt};
  }

  void expMSSA::wcorrPNG()
  {
#ifdef HAVE_LIBPNGPP
    {
      const int minSize = 400;
      int ndup = 1;
      int nDim = std::min<int>(numW, npc);
      if (nDim < minSize) ndup = minSize/nDim + 1;
      if (ndup > 1) std::cout << "Pixel duplication=" << ndup << std::endl;

      bool use_dist = false;
      if (params["distance"]) use_dist = true;

      for (auto u : RC) {
	Eigen::MatrixXd wc = wCorrKey(u.first);

	png::image< png::rgb_pixel > image(nDim*ndup, nDim*ndup);
	ColorGradient color;

	for (size_t y = 0; y < nDim; y++) {
	  for (size_t x = 0; x < nDim; x++) {
	    double val = wc(y, x);
	    if (use_dist) val = 1.0 - sqrt(val);
	    png::rgb_pixel cval = color(val);
	    for (size_t yy = y*ndup; yy < (y+1)*ndup; yy++) {
	      for (size_t xx = x*ndup; xx < (x+1)*ndup; xx++) {
		if (flip) image[image.get_height()-1-yy][xx] = cval;
		else      image[yy][xx] = cval;

	      }
	    }
	  }
	}

	std::ostringstream sout; sout << prefix + ".wcorr" + "_"
				      << u.first << ".png";

	if (verbose)
	  std::cout << "Attempting to write: " << sout.str() << std::endl;

	image.write(sout.str());
      }

      {
	Eigen::MatrixXd wc = wCorrAll();

	png::image< png::rgb_pixel > image(nDim*ndup, nDim*ndup);
	ColorGradient color;

	for (size_t y = 0; y < nDim; y++) {
	  for (size_t x = 0; x < nDim; x++) {
	    double val = wc(y, x);
	    if (use_dist) val = 1.0 - sqrt(val);

	    png::rgb_pixel cval = color(val);

	    for (size_t yy = y*ndup; yy < (y+1)*ndup; yy++) {
	      for (size_t xx = x*ndup; xx < (x+1)*ndup; xx++) {
		if (flip) image[image.get_height()-1-yy][xx] = cval;
		else      image[yy][xx] = cval;

	      }
	    }
	  }
	}

	std::ostringstream sout; sout << prefix + ".wcorr_allchan.png";

	std::cout << "Attempting to write: " << sout.str() << std::endl;

	image.write(sout.str());
      }

    }
#else
    std::cout << "PNG is not available so I can't make w-correlation images"
	      << std::endl;
#endif
  }

  void expMSSA::kmeans(int clusters, bool toTerm, bool toFile)
  {
    if (clusters==0) {
      std::cout << "expMSSA::kmeans: you need clusters>0" << std::endl;
      return;
    }

    if (toTerm) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- K-means group analysis"              << std::endl;
      std::cout << "---- numT=" << numT << " numW=" << numW   << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }

    std::ofstream out;
    std::string filename;

    if (toFile) {
      filename = prefix + ".kmeans";
      out.open(filename);
      if (not out)
	std::cerr << "Error opening file <" << filename << ">" << std::endl;
    }

    // W-correlation-based distance functor
    //
    KMeans::WcorrDistance dist(numT, numW);

    for (auto u : mean) {
      // Pack point array
      //
      std::vector<KMeans::Ptr> data;
      for (int j=0; j<ncomp; j++) {
	data.push_back(std::make_shared<KMeans::Point>(numT));
	for (int i=0; i<numT; i++) data.back()->x[i] = RC[u.first](i, j);
      }

      // Initialize k-means routine
      //
      KMeans::kMeansClustering kMeans(data);

      // Run 100 iterations
      //
      kMeans.iterate(dist, 100, clusters, 2, false);

      // Retrieve cluster associations
      //
      auto results = kMeans.get_results();

      // Write to file
      //
      if (out) {
	out << std::string(60, '-') << std::endl
	    << " *** n=" << u.first << std::endl
	    << std::string(60, '-') << std::endl;

	for (int j=0; j<results.size(); j++) {
	  out << std::setw(6)  << j
	      << std::setw(12) << std::get<1>(results[j])
	      << std::endl;
	}
      }

      // Write to stdout
      //
      if (toTerm) {
	std::cout << std::string(60, '-') << std::endl
		  << " *** n=" << u.first << std::endl
		  << std::string(60, '-') << std::endl;

	for (int j=0; j<results.size(); j++) {
	  std::cout << std::setw(6)  << j
		    << std::setw(12) << std::get<1>(results[j])
		    << std::endl;
	}
      }
    }

    if (params["allchan"]) {

      // Pack point array
      //
      std::vector<KMeans::Ptr> data;
      int sz = mean.size();
      for (int j=0; j<ncomp; j++) {
	data.push_back(std::make_shared<KMeans::Point>(numT*sz));
	int c = 0;
	for (auto u : mean) {
	  for (int i=0; i<numT; i++) data.back()->x[c++] = RC[u.first](i, j);
	}
      }

      // Initialize k-means routine
      //
      KMeans::kMeansClustering kMeans(data);

      // Run 100 iterations
      //
      KMeans::WcorrDistMulti dist2(numT, numW, sz);
      kMeans.iterate(dist2, 100, clusters, 2, false);

      // Retrieve cluster associations
      //
      auto results = kMeans.get_results();

      // Write to file
      //
      if (out) {
	out << std::string(60, '-') << std::endl
	    << " *** total"         << std::endl
	    << std::string(60, '-') << std::endl;

	for (int j=0; j<results.size(); j++) {
	  out << std::setw(6) << j
	      << std::setw(9) << std::get<1>(results[j])
	      << std::endl;
	}
      }

      // Write to stdout
      //
      if (toTerm) {
	std::cout << std::string(60, '-') << std::endl
		  << " *** total"         << std::endl
		  << std::string(60, '-') << std::endl;

	for (int j=0; j<results.size(); j++) {
	  std::cout << std::setw(6) << j
		    << std::setw(9) << std::get<1>(results[j])
		    << std::endl;
	}
      }

    }

    if (!out) {
      if (toFile)
	std::cout << "Bad output stream for <" << filename << ">" << std::endl;
    }
    out.close();

  }

  std::map<std::string, CoefClasses::CoefsPtr> expMSSA::getReconstructed(bool reconstructmean)
  {
    if (not reconstructed) {
      std::ostringstream sout;
      if (not computed)
	sout << "expMSSA::getReconstructed(): "
	     << "call eigenvalues() or getPC() followed by "
	     << "reconstruct() before getReconstructed()";
      else
	sout << "expMSSA::getReconstructed(): "
	     << "call reconstruct() before getReconstructed()";

      std::runtime_error(sout.str());
    }


    // Copy the original map for return
    //
    auto newdata = data;

    for (int i=0; i<numT; i++) {
      for (auto u : mean) {
	double disp = totVar;
	if (type == TrendType::totPow) disp = totPow;
	if (disp==0.0) disp = var[u.first];

	double acc = 0.0;
	for (int j=0; j<ncomp; j++) acc += RC[u.first](i, j);


  	newdata[u.first][i] = acc*disp;
  if (reconstructmean) {
      	if (useMean) newdata[u.first][i] += u.second;
      }
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
  expMSSA::valid_keys = {
    "verbose",
    "writeCov",
    "Jacobi",
    "BDCSVD",
    "Traj",
    "RedSym",
    "allchan",
    "distance",
    "flip",
    "power",
    "evtol",
    "output",
    "totVar",
    "totPow",
    "noMean"
  };

  void expMSSA::assignParameters(const std::string flags)
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
	throw YamlConfigError("MSSA::expMSSA", "parameter", unmatched, __FILE__, __LINE__);

      // Compute flags
      //
      computed      = false;
      reconstructed = false;

      // Top level parameter flags
      //
      verbose  = bool(params["verbose"   ]);
      flip     = bool(params["flip"      ]);
      powerf   = bool(params["power"     ]);

      if (params["evtol"]  ) evtol    = params["evtol"].as<double>();
      else                   evtol    = 0.01;

      if (params["output"] ) prefix   = params["output"].as<std::string>();
      else                   prefix   = "exp_mssa";

    }
    catch (const YAML::ParserException& e) {
      std::cout << "expMSSA::assignParameters, parsing error=" << e.what()
		<< std::endl;
      throw;
    }
  }


  // Save current MSSA state to an HDF5 file with the given prefix
  void expMSSA::saveState(const std::string& prefix)
  {
    if (not computed) return;	// No point in saving anything

    if (std::filesystem::exists(prefix + "_mssa.h5")) {
      std::ostringstream sout;
      sout << "expMSSA::saveState: the file <" <<  prefix + "_mssa.h5"
	   << "> already exists.\nPlease delete this file or choose a "
	   << "different file name";
      throw std::runtime_error(sout.str());
    }

    try {
      // Create a new hdf5 file
      //
      HighFive::File file(prefix + "_mssa.h5",
			  HighFive::File::ReadWrite |
			  HighFive::File::Create);

      // Write the time dimension
      //
      file.createAttribute<int>("numT", HighFive::DataSpace::From(numT)).write(numT);

      // Write the number of channels
      //
      file.createAttribute<int>("nKeys", HighFive::DataSpace::From(nkeys)).write(nkeys);

      // Write the window size
      //
      file.createAttribute<int>("numW", HighFive::DataSpace::From(numW)).write(numW);

      // Write the number of PCs
      //
      file.createAttribute<int>("numPC", HighFive::DataSpace::From(npc)).write(npc);

      // Save trend state
      //
      int trend = static_cast<std::underlying_type<TrendType>::type>(type);
      file.createAttribute<bool>("trendType", HighFive::DataSpace::From(trend)).write(trend);

      // Save the key list
      //
      std::vector<Key> keylist;
      for (auto k : mean) keylist.push_back(k.first);

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

      // Save mssa_analysis state
      //
      HighFive::Group analysis = file.createGroup("mssa_analysis");

      analysis.createDataSet("Y",  Y );
      analysis.createDataSet("S",  S );
      analysis.createDataSet("U",  U );
      analysis.createDataSet("PC", PC);

      // Save reconstruction
      //
      if (reconstructed) {
	HighFive::Group recon = file.createGroup("reconstruction");

	recon.createAttribute<int>   ("ncomp",  HighFive::DataSpace::From(ncomp) ).write(ncomp);
	recon.createAttribute<double>("totVar", HighFive::DataSpace::From(totVar)).write(totVar);
	recon.createAttribute<double>("totPow", HighFive::DataSpace::From(totVar)).write(totPow);

	for (int n=0; n<keylist.size(); n++) {
	  std::ostringstream scnt;
	  scnt << "RC_" << n;
	  recon.createDataSet(scnt.str(), RC[keylist[n]]);
	}
      }

    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }
  }

  // Restore current MSSA state to an HDF5 file with the given prefix
  void expMSSA::restoreState(const std::string& prefix)
  {
    try {
      // Silence the HDF5 error stack
      //
      HighFive::SilenceHDF5 quiet;

      // Try opening the file as HDF5
      //
      HighFive::File h5file(prefix + "_mssa.h5", HighFive::File::ReadOnly);

      // Read and check parameters
      //
      int nTime, nKeys, numW1, numPC, trend;

      h5file.getAttribute("numT" ).read(nTime);
      h5file.getAttribute("nKeys").read(nKeys);
      h5file.getAttribute("numW" ).read(numW1);
      h5file.getAttribute("numPC").read(numPC);
      h5file.getAttribute("trendType").read(trend);

      // Number of channels
      //
      nkeys = mean.size();

      // Test recovered parameters
      //
      if (nTime != numT) {
	std::ostringstream sout;
	sout << "expMSSA::restoreState: saved state has numT="
	     << nTime << " but expMSSA expects numT=" << numT
	     << ".\nCan't restore mssa state!";
	throw std::runtime_error(sout.str());
      }

      if (nKeys != nkeys) {
	std::ostringstream sout;
	sout << "expMSSA::restoreState: saved state has nkeys="
	     << nKeys << " but expMSSA expects nkeys=" << nkeys
	     << ".\nCan't restore mssa state!";
	throw std::runtime_error(sout.str());
      }

      if (numW != numW1) {
	std::ostringstream sout;
	sout << "expMSSA::restoreState: saved state has numW="
	     << numW1 << " but expMSSA expects numW=" << numW
	     << ".\nCan't restore mssa state!";
	throw std::runtime_error(sout.str());
      }

      if (numPC != npc) {
	std::ostringstream sout;
	sout << "expMSSA::restoreState: saved state has npc="
	     << numPC << " but expMSSA expects npc=" << npc
	     << ".\nCan't restore mssa state!";
	throw std::runtime_error(sout.str());
      }

      int ttype = static_cast<std::underlying_type<TrendType>::type>(type);

      if (trend != ttype) {
	std::ostringstream sout;
	sout << "expMSSA::restoreState: saved state has trend="
	     << getTrendType.at(trend) << " [" << trend
	     << "] but expMSSA expects trend="
	     << getTrendType.at(ttype) << " [" << ttype
	     << "].\nCan't restore mssa state!";
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
	auto it = mean.find(keylist[n]);
	if (it == mean.end()) bad = true;
	else if (it->first.size() != keylist[n].size()) bad = true;
	else {
	  for (int i=0; i<keylist[n].size(); i++) {
	    if (keylist[n][i] != it->first[i]) bad = true;
	  }
	}
      }

      if (bad) {
	std::ostringstream sout;
	sout << "expMSSA::restoreState: keylist mismatch." << std::endl
	     << "Can't restore mssa state! Wanted keylist: ";
	for (auto v : mean) {
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

      auto analysis = h5file.getGroup("mssa_analysis");

      analysis.getDataSet("Y" ).read(Y );
      analysis.getDataSet("S" ).read(S );
      analysis.getDataSet("U" ).read(U );
      analysis.getDataSet("PC").read(PC);

      numK = numT - numW + 1;	// Recompute numK, needed for
				// reconstruction
      computed = true;

      if (h5file.exist("reconstruction")) {
	auto recon = h5file.getGroup("reconstruction");

	recon.getAttribute("ncomp" ).read(ncomp);
	recon.getAttribute("totVar").read(totVar);
	recon.getAttribute("totPow").read(totPow);

	for (int n=0; n<keylist.size(); n++) {
	  std::ostringstream scnt;
	  scnt << "RC_" << n;
	  recon.getDataSet(scnt.str()).read(RC[keylist[n]]);
	}

	reconstructed = true;
      }

    } catch (HighFive::Exception& err) {
      std::cerr << "**** Error opening or reading H5 file ****" << std::endl;
      throw;
    }

  }

  expMSSA::expMSSA(const mssaConfig& config, int nW, int nPC, const std::string flags) : numW(nW), npc(nPC)
  {
    // Parse the YAML string
    //
    assignParameters(flags);

    // Detrending style (totVar and totPow are only useful, so far, for
    // noise computation)
    //
    if (params["totVar"])
      type = TrendType::totVar;
    else if (params["totPow"])
      type = TrendType::totPow;
    else
      type = TrendType::perChannel;


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

      mean[key] = 0.0;
      var [key] = 0.0;

      for (auto y : coefDB.getData(key)) {
	mean[key] += y;
	var [key] += y*y;
	data[key].push_back(y);
      }
    }

    // Normalize and detrend
    //
    totVar = totPow = 0.0;
    useMean = true;

    if (type == TrendType::totPow) {

      // Compute total power from the coefficient database
      //
      for (auto k : coefDB.getKeys()) {
	for (auto v : coefDB.getData(k)) {
	  totPow += v*v;
	}
      }

      // Average total power per time slice
      //
      totPow = std::sqrt(totPow/numT + std::numeric_limits<double>::min());

      for (auto & u : mean) {
	Key k = u.first;
	//--------------
	mean[k] /= numT;
	var[k]   = var[k]/numT - mean[k]*mean[k];
	totVar  += var[k];
	var[k]   = sqrt(fabs(var[k]) + std::numeric_limits<double>::min());
      }

      if (params["noMean"]) useMean = false;

      for (auto & u : mean) {
	Key k = u.first;
	for (auto & v : data[k]) {
	  if (useMean) v -= mean[k];
	  v /= totPow;
	}
      }

    } else if (type == TrendType::totVar) {

      for (auto & u : mean) {
	Key k = u.first;
	//--------------
	mean[k] /= numT;
	var[k]   = std::max<double>(var[k]/numT - mean[k]*mean[k], 0.0);
	totVar  += var[k];
	var[k]   = std::sqrt(fabs(var[k]) + std::numeric_limits<double>::min());
      }

      totVar = std::sqrt(totVar+std::numeric_limits<double>::min());

      for (auto & u : mean) {
	Key k = u.first;
	for (auto & v : data[k]) {
	  v -= mean[k];
	  if (totVar>0.0) v /= totVar;
	}
      }
    } else {
      // the default detrending, by mean and variance, type = TrendType::perChannel
      for (auto & u : mean) {
	Key k = u.first;
	//--------------
	mean[k] /= numT;
	var[k]   = std::max<double>(var[k]/numT - mean[k]*mean[k], 0.0);
	var[k]   = std::sqrt(fabs(var[k])+ std::numeric_limits<double>::min());
	for (auto & v : data[k]) {
	  v -= mean[k];
	  if (var[k]>0.0) v /= var[k];
	}
      }
    }
    computed      = false;
    reconstructed = false;
  }
  // END expMSSA constructor

}
// END namespace MSSA
