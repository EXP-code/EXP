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

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cmath>
#include <map>

#include <Eigen/Dense>

/* For debugging
   #undef eigen_assert
   #define eigen_assert(x)						\
   if (!(x)) { throw (std::runtime_error("Eigen assert error")); }
*/

#include <omp.h>

#include <TransformFFT.H>
#include <ColorGradient.H>
#include <KMeans.H>

#include "CoefDB.H"

#include "RedSVD.H"
#include "config.h"
#include "YamlConfig.H"
#include "libvars.H"
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
  
  
  // Helper ostream manipulator
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
  
  
  // Do the SVD
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
	RedSVD::RedSVD<Eigen::MatrixXd> svd(cov, rank);
	S = svd.singularValues();
	U = svd.matrixU();
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
    
    if (chatty) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- Eigenvalues"                         << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
    
    double tot = 0.0;
    for (int j=0; j<S.size(); j++) tot += S(j);
    if (chatty) std::cout << "Total=" << tot << std::endl;
    
    if (chatty) {
      double cum = 0.0;
      for (int j=0; j<S.size(); j++) {
	std::cout << std::setw( 5) << j
		  << std::setw(18) << S(j)
		  << std::setw(18) << (cum += S(j))/tot
		  << std::endl;
	if (1.0 - cum/tot < evtol) break;
      }
    }
      
    if (dfiles) {
	
      std::string filename = prefix + ".ev";
      std::ofstream out(filename);
      if (out) {
	double cum = 0.0;
	for (int j=0; j<S.size(); j++) {
	  out << std::setw( 5) << j
	      << std::setw(18) << S(j)
	      << std::setw(18) << (cum += S(j))/tot
	      << std::endl;
	}
	out.close();
      } else {
	std::cout << "Could not open <" << filename << ">" << std::endl;
	exit(-1);
      }
    }
	
    npc = std::min<int>(npc, numW*nkeys);
	
    if (chatty) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- Eigenvectors"                        << std::endl;
      std::cout << "----------------------------------------" << std::endl;
      for (int j=0; j<std::min<int>(numW*nkeys, rank); j++) {
	std::cout << std::setw(5) << j;
	for (int i=0; i<numW*nkeys; i++)
	  std::cout << std::setw(15) << std::setprecision(6) << U(i, j);
	std::cout << std::endl;
      }
    }
	
    if (dfiles) {

      // The eigenvectors {E} are composed of nkeys consecutive segments of
      // of length numW
      //
      std::string filename = prefix + ".evec";
      std::ofstream out(filename);
      if (out) {
	for (int j=0; j<std::min<int>(numW*nkeys, rank); j++) {
	  out << std::setw(5) << j;
	  for (int i=0; i<numW*nkeys; i++) out << std::setw(15) << U(i, j);
	  out << std::endl;
	}
	out.close();
      } else {
	std::cout << "Could not open <" << filename << ">" << std::endl;
	exit(-1);
      }
    }
	
    // Compute the PCs by projecting the data
    //
    PC = Y * U;
	
    if (chatty) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- Principal components"                << std::endl;
      std::cout << "----------------------------------------" << std::endl;
      
      for (int i=0; i<numK; i++) {
	std::cout << std::setw(15) << std::setprecision(6) << coefDB.times[i];
	for (int j=0; j<PC.cols(); j++)
	  std::cout << std::setw(15) << std::setprecision(6) << PC(i, j);
	std::cout << std::endl;
      }
    }
	
    if (dfiles) {
      std::string filename = prefix + ".pc";
      std::ofstream out(filename);
      if (out) {
	for (int i=0; i<numK; i++) {
	  out << std::setw(5) << coefDB.times[i];
	  for (int j=0; j<PC.cols(); j++)
	    out << std::setw(15) << PC(i, j);
	  out << std::endl;
	}
	out.close();
      } else {
	std::cout << "Could not open <" << filename << ">" << std::endl;
	exit(-1);
      }
    }
	
    computed = true;
    reconstructed = false;
  }
    
  void expMSSA::reconstruct(const std::vector<int>& evlist)
  {
    // Reconstructed time series
    //
    ncomp = std::min<int>({numW, npc, static_cast<int>(PC.cols())});
    
    // Make a bool vector
    //
    std::vector<bool> I(ncomp, false);
    if (evlist.size()) {
      for (auto v : evlist) if (v<ncomp) I[v] = true;
    } else I = std::vector<bool>(ncomp, true);

    // Deduce the rank
    //
    int rank  = std::min<int>({static_cast<int>(Y.cols()), static_cast<int>(Y.rows()), npc});
    
    
    for (auto u : mean) RC[u.first].resize(numT, ncomp);
    
    // Embedded time series matrix
    //
    Eigen::MatrixXd  Z(numT, numW);
    
    // Split eigenvectors into channel sequences
    //
    std::map<Key, Eigen::MatrixXd> rho;
    
    int n = 0;
    for (auto u : mean) {
      rho[u.first].resize(numW, rank);
      rho[u.first].fill(0.0);
      for (int i=0; i<numW; i++) {
	for (int j=0; j<ncomp; j++)
	  if (I[j]) rho[u.first](i, j) = U(numW*n + i, j);
      }
      n++;
    }
    
    if (chatty) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- Reconstruction"                      << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
    

#pragma omp parallel
    {
      // Parallelize the map iteration by wrapping it in a standard loop.
      // Ungainly for sure.
      //
      int thread_count = omp_get_num_threads();
      int thread_num   = omp_get_thread_num();
      size_t map_size  = mean.size();

      auto u = mean.begin();
      std::advance(u, thread_num);

      for (int q=thread_num; q<map_size; q+=thread_count) {

	for (int w=0; w<ncomp; w++) {
	
	  Z.fill(0.0);
	
	  // Build reverse embedded time series.
	  //
	  for (int j=0; j<numW; j++) {
	    for (int i=0; i<numT; i++)
	      if (i - j >= 0 and i - j < numK) Z(i, j) = PC(i - j, w);
	  }
	  
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
	    for (int j=L; j<U; j++)
	      RC[u->first](i, w) += Z(i, j) * rho[u->first](j, w) * W;
	  }
	}

	if( q+thread_count < map_size ) std::advance(u, thread_count);
      }
    }

    /*
    for (auto u : mean) {

      for (int w=0; w<ncomp; w++) {
	
	Z.fill(0.0);
	
	// Build reverse embedded time series.
	//
	for (int j=0; j<numW; j++) {
	  for (int i=0; i<numT; i++)
	    if (i - j >= 0 and i - j < numK) Z(i, j) = PC(i - j, w);
	}
	  
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
	  
	  RC[u.first](i, w) = 0.0;
	  for (int j=L; j<U; j++)
	    RC[u.first](i, w) += Z(i, j) * rho[u.first](j, w) * W;
	}
      }
    }
    */

    if (chatty) {
      for (int i=0; i<numT; i++) {
	std::cout << std::setw(15) << std::setprecision(6) << coefDB.times[i];
	for (auto u : mean) {
	  double acc = 0.0;
	  for (int j=0; j<ncomp; j++) acc += RC[u.first](i, j);
	  std::cout << std::setw(15) << std::setprecision(6) << acc;
	}
	std::cout << std::endl;
      }
    }
    
    
    if (dfiles) {
      for (auto u : mean) {
	std::ostringstream sout;
	sout << prefix;
	for (auto v : u.first) sout << "_" << v;
	sout << ".elem";
	std::ofstream out(sout.str());
	if (out) {
	  double disp = totVar;
	  if (type == TrendType::totPow) disp = totPow;
	  if (disp==0.0) disp = var[u.first];
	  
	  for (int i=0; i<numT; i++) {
	    out << std::setw(15) << std::setprecision(6) << coefDB.times[i];
	    for (int w=0; w<ncomp; w++)
	      out << std::setw(15) << std::setprecision(6)
		// << RC[u.first](i, w)*disp + u.second;
		  << RC[u.first](i, w);
	    out << std::endl;
	  }
	} else {
	  std::cout << "Could not open <" << sout.str() << ">" << std::endl;
	  exit(-1);
	}
      }
    }

    reconstructed = true;
  }
  
  
  void expMSSA::contributions()
  {
    if (not computed) {
      std::cout << "expMSSA::contributions: "
		<< "call eigenvalues() or getPC() before contributions()"
		<< std::endl;
      return;
    }

#ifdef HAVE_LIBPNGPP
    if (chatty) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- Elemental fractions"                 << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
    
    std::string filename = prefix + ".f_contrib";
    std::ofstream out(filename);
    
    // Used in both computations
    //
    std::map<Key, std::vector<double> > values;
    
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
      
      // For each series, compute ||a^k_j||
      //
      std::map<Key, double> norm;
      
      for (auto u : mean) {
	auto key = u.first;
	values[key].resize(npc, 0.0);
	
	for (int j=0; j<npc; j++) {
	  for (int i=0; i<numT; i++) {
	    values[key][j] += RC[u.first](i, j)*RC[u.first](i, j);
	  }
	}
      }
      
      // This is norm for each series over the entire reconstruction
      // (that is, fixed coefficient channel, summed over all PCs)
      for (auto u : mean) {
	auto key = u.first;
	norm[key] = 0.0;
	for (int j=0; j<npc; j++) {
	  norm[key] += values[key][j];
	}
      }
      
      double maxVal = 0.0;
      for (auto u : mean) {
	auto key = u.first;
	for (auto q : key) out << std::setw(4) << q;
	
	for (int j=0; j<npc; j++) {
	  out << std::setw(18) << sqrt(values[key][j]/norm[key]);
	  maxVal = std::max<double>(maxVal, values[key][j]/norm[key]);
	}
	out << std::endl;
      }
      out.close();
      
      int i = 0;
      for (auto u : mean) {
	auto key = u.first;
	for (int j=0; j<npc; j++) {
	  png::rgb_pixel cval = color(sqrt(values[key][j]/norm[key]/maxVal));
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
      
      // Now, norm is sum over series for each PC
      //
      for (auto u : mean) {
	auto key = u.first;
	for (int j=0; j<npc; j++) {
	  norm[j] += values[key][j];
	}
      }
      
      double maxVal = 0.0;
      for (int j=0; j<npc; j++) {
	out << std::setw(8) << j;
	int i=0;
	for (auto u : mean) {
	  auto key = u.first;
	  out << std::setw(18) << sqrt(values[key][j]/norm[j]);
	  maxVal = std::max<double>(maxVal, values[key][j]/norm[j]);
	}
	out << std::endl;
      }
      out.close();
      
      for (int j=0; j<npc; j++) {
	int i=0;
	for (auto u : mean) {
	  auto key = u.first;
	  png::rgb_pixel cval = color(sqrt(values[key][j]/norm[j]/maxVal));
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
  
  
  Eigen::MatrixXd expMSSA::pcDFT(Eigen::VectorXd& F, Eigen::VectorXd& P)
  {
    
    if (chatty) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- Principal component periodogram"     << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
    
    Eigen::MatrixXd pw;
    {
      double DT = coefDB.times[1] - coefDB.times[0];
      
      int nfreq = numT/2 + 1;

      pw.resize(nfreq, npc);
      Eigen::VectorXd in(numK);
      
      pw.setZero();
      
      for (int j=0; j<npc; j++) {
	for (int i=0; i<numK; i++) in(i) = PC(i, j);
	TransformFFT fft(DT, in);
	fft.Power(F, P);
	for (int i=0; i<nfreq; i++) pw(i, j) = P(i);
      }
      
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
    
    return pw;
  }
  
  Eigen::MatrixXd expMSSA::channelDFT(Eigen::VectorXd& F, Eigen::VectorXd& P)
  {
    if (chatty) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- Coefficient periodogram"             << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
    
    Eigen::MatrixXd pw;
    {
      double DT = coefDB.times[1] - coefDB.times[0];
      
      int nfreq = numT/2 + 1;
      
      pw.resize(numW, nfreq);
      Eigen::VectorXd p0(nfreq), in(numT);
      Eigen::VectorXd pt = Eigen::VectorXd::Zero(nfreq);
      
      for (auto u : mean) {
	
	for (int i=0; i<numT; i++) {
	  in(i) = 0.0;
	  for (int j=0; j<ncomp; j++) in(i) += RC[u.first](i, j);
	}
	TransformFFT fft(DT, in);
	fft.Power(F, p0);
	
	for (int j=0; j<ncomp; j++) {
	  for (int i=0; i<numT; i++) in(i) = RC[u.first](i, j);
	  
	  TransformFFT fft(DT, in);
	  fft.Power(F, P);
	  for (int k=0; k<nfreq; k++) {
	    pt(k)   += P(k);
	    pw(j, k) = P(k);
	  }
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
	  for (int j=0; j<nfreq; j++) {
	    std::ostringstream sout; sout << "PC " << j;
	    out << std::setw(15) << sout.str();
	  }
	  out << "# " << std::setw(13) << "[1]"
	      << std::setw(15) << "[2]"
	      << std::setw(15) << "[3]"
	      << std::setw(15) << "[4]";
	  for (int j=0; j<nfreq; j++) {
	    std::ostringstream sout; sout << '[' << j+5 << ']';
	    out << std::setw(15) << sout.str();
	  }
	  out << std::endl;
	  
	  for (int j=0; j<nfreq; j++) {
	    out << std::setw(15) << std::setprecision(6) << F(j)
		<< std::setw(15) << std::setprecision(6) << 2.0*M_PI/F(j)
		<< std::setw(15) << std::setprecision(6) << pt(j)
		<< std::setw(15) << std::setprecision(6) << p0(j);
	    for (int k=0; k<numW; k++)
	      out << std::setw(15) << std::setprecision(6) << pw(k, j);
	    out << std::endl;
	  }
	  out.close();
	} else {
	  std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	  exit(-1);
	}
      }
    }
    
    return pw;
  }
  
  void expMSSA::wcorrPNG()
  {
#ifdef HAVE_LIBPNGPP
    if (chatty) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- w-correlation"                       << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
    
    {
      const int minSize = 400;
      int ndup = 1;
      int nDim = std::min<int>(numW, npc);
      if (nDim < minSize) ndup = minSize/nDim + 1;
      if (ndup > 1) std::cout << "Pixel duplication=" << ndup << std::endl;
      
      bool use_dist = false;
      if (params["distance"]) use_dist = true;
      
      for (auto u : mean) {
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
      
    if (out) {
      if (chatty)
	std::cout << "Successfully wrote: <" << filename << ">" << std::endl;
    } else {
      if (toFile)
	std::cout << "Bad output stream for <" << filename << ">" << std::endl;
    }
    out.close();
    
  }
  
  
  std::map<std::string, Coefs::CoefsPtr> expMSSA::getReconstructed(bool zero)
  {
    if (dfiles) {

      std::string filename = prefix + ".recon";
      std::ofstream out(filename);
    
      if (out) {
	for (int i=0; i<numT; i++) {
	  out << std::setw(15) << std::setprecision(6) << coefDB.times[i];
	  for (auto u : mean) {
	    double disp = totVar;
	    if (type == TrendType::totPow) disp = totPow;
	    if (disp==0.0) disp = var[u.first];
	    
	    double acc = 0.0;
	    for (int j=0; j<ncomp; j++) acc += RC[u.first](i, j);
	    out << std::setw(15) << std::setprecision(6) << acc*disp + u.second;
	  }
	  out << std::endl;
	}
	out.close();
      } else {
	std::cout << "Could not open <" << filename << ">" << std::endl;
	exit(-1);
      }
      
    
      filename = prefix + ".recon_diff";
      out.open(filename);
      if (out) {
	for (int i=0; i<numT; i++) {
	  out << std::setw(15) << std::setprecision(6) << coefDB.times[i];
	  for (auto k : mean) {
	    double disp = totVar;
	    if (type == TrendType::totPow) disp = totPow;
	    if (disp==0.0) disp = var[k.first];
	    
	    double acc = 0.0;
	    for (int j=0; j<ncomp; j++) acc += RC[k.first](i, j);
	    out << std::setw(15) << std::setprecision(6) << acc*disp + k.second
		<< std::setw(15) << std::setprecision(6) << data[k.first][i]*disp + k.second;
	  }
	  out << std::endl;
	}
	out.close();
      } else {
	std::cout << "Could not open <" << filename << ">" << std::endl;
	exit(-1);
      }
    
      if (params["triples"]) {
      
	const char lab[] = {'c', 's'};
      
	for (auto u : mean) {
	  std::ostringstream fname1, fname2, suffix;
	  auto k = u.first;
	  suffix << "key";
	  for (auto q : k) suffix << "_" << q;
	  fname1 << prefix << ".recon."       << suffix.str();
	  fname2 << prefix << ".recon_accum." << suffix.str();
	  std::ofstream out1 (fname1.str());
	  std::ofstream out2 (fname2.str());
	  if (out1 and out2) {
	    double disp = totVar;
	    if (type == TrendType::totPow) disp = totPow;
	    if (disp==0.0) disp = var[k];
	  
	    out1 << "# " << u.first << std::endl;
	    out2 << "# " << u.first << std::endl;
	    for (int i=0; i<numT; i++) {
	      out1 << std::setw(15) << std::setprecision(6) << coefDB.times[i];
	      out2 << std::setw(15) << std::setprecision(6) << coefDB.times[i];
	      double accum = 0.0;
	      for (int j=0; j<ncomp; j++) {
		out1 << std::setw(15) << std::setprecision(6)
		     << RC[k](i, j)*disp + u.second;
		out2 << std::setw(15) << std::setprecision(6)
		     << (accum += (RC[k](i, j)*disp + u.second));
	      }
	      out1 << std::endl;
	      out2 << std::endl;
	    }
	    out1.close();
	    out2.close();
	  } else {
	    std::cout << "Could not open <" << fname1.str() << ">"
		      << " or <" << fname2.str() << ">" << std::endl;
	    exit(-1);
	  }
	}
      }
    }
    // END dfiles block
    
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
	
	newdata[u.first][i] = acc*disp + u.second;
      }
    }
    
    
    // Use CoefDB to recover names and use newdata to replace data
    //
    coefDB.beginUpdate(zero);
    for (auto v : newdata) {
      if (verbose) std::cout << "Updating for: " << v.first << std::endl;
      coefDB.setData(v.first, v.second);
    }
    
    // Return updated namestr-coefficient map
    //
    return coefDB.endUpdate();
  }
  
  void expMSSA::assignParameters(const std::string flags)
  {
    // Parse the parameters database
    //
    try {
      
      // Load parameters from string
      //
      params = YAML::Load(flags);
      
      // Compute flags
      //
      computed      = false;
      reconstructed = false;
      
      // Top level parameter flags
      //
      chatty   = bool(params["chatty"    ]);
      verbose  = bool(params["verbose"   ]);
      flip     = bool(params["flip"      ]);
      dfiles   = bool(params["writeFiles"]);
      
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
    // if (chatty)
      std::cout << "Eigen is using " << Eigen::nbThreads()
		<< " threads" << std::endl;
    
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
    
    // Print channel info
    //
    if (chatty) {
      std::cout << std::string(60, '-') << std::endl
		<< "Using coefficient keys (l, m, n, c/s)_k or (M, n, c/s)_k"
		<< std::endl
		<< "for spherical (l, m) and cylindrical (M) bases."
		<< std::endl
		<< "c/s denotes cosine/sine and component index is k."
		<< std::endl
		<< std::string(60, '-') << std::endl
		<< std::setw(18) << std::right << "Key"
		<< std::setw(18) << std::right << "Mean values"
		<< std::endl << std::string(60, '-') << std::endl;
      
      // Print key, mean pairs
      //
      for (auto v : mean) std::cout << std::setw(18) << v.first
				    << std::setw(18) << v.second
				    << std::endl;
      std::cout << std::string(60, '-') << std::endl;
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
      
      if (chatty) std::cout << "<total Power>=" << totPow << std::endl;
      
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
      
    } else if (params["totVar"]) {
      
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
    
    if (dfiles) {
      std::string filename(prefix + ".data");
      std::ofstream out(filename);
      
      for (int i=0; i<numT; i++) {
	out << std::setw(6) << coefDB.times[i];
	for (auto u : mean) out << std::setw(18) << data[u.first][i];
	out << std::endl;
      }
      out.close();
    }
  }
  // END expMSSA constructor

}
// END namespace MSSA
  
