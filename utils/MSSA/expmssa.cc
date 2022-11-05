//
// M-SSA for EXP coefficients
//
// Updated to use fixed-rank approximation for the SVD to save time
// and space.  Uses the approximate SVD randomized algorithm from
// Halko, Martinsson, and Tropp by default.  Use BDCSVD or Jacobi
// flags for the standard methods.
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

#include <TransformFFT.H>
#include <ColorGradient.H>
#include <KMeans.H>

#include <CoefDB.H>

#include <RedSVD.H>
#include <config_exp.h>
#include <YamlConfig.H>
#include <libvars.H>

Eigen::MatrixXd wCorr(Eigen::MatrixXd & R)
{
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


Eigen::MatrixXd wCorrAll(std::map<CoefDB::Key, Eigen::MatrixXd>& RC)
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


//! Helper ostream manipulator
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

int
main(int argc, char **argv)
{
  //--------------------------------------------------
  // Parameters
  //--------------------------------------------------
  std::string prefix, grpfile, config, spec;
  int numW, nmin, nmax, npc;
  int stride, skip, clusters;
  std::vector<int> MM;
  double evtol;

  //--------------------------------------------------
  // Command-line parsing
  //--------------------------------------------------

  std::string cmd_line;
  for (int i=0; i<argc; i++) {
    cmd_line += argv[i];
    cmd_line += " ";
  }

  cxxopts::Options options("expmssa","This routine uses M-SSA to analyze basis coefficients\n\nOptions");

  options.add_options()
    ("P,npc",        "maximum number of principle components",
     cxxopts::value<int>(npc)->default_value("99999"))
    ("s,skip",       "number of time steps to skip on reading input",
     cxxopts::value<int>(skip)->default_value("1"))
    ("S,stride",     "number of time steps between coefficient outputs",
     cxxopts::value<int>(stride)->default_value("100"))
    ("t,evtol",      "cut on cumulative variance for eigenvalue sum",
     cxxopts::value<double>(evtol)->default_value("0.01"))
    ("G,group",      "do a full reconstruction for visualization with given group file",
     cxxopts::value<std::string>(grpfile)->default_value("group.list"))
    ("kmeans",       "use k-means analysis for experimental grouping using the w-correlation",
     cxxopts::value<int>(clusters)->default_value("0"))
    ("o,output",     "output file prefix",
     cxxopts::value<std::string>(prefix)->default_value("exp_mssa"))
    ("W,numW",       "window size",
     cxxopts::value<int>(numW))
    ("c,config",      "Input parameter config file",
     cxxopts::value<std::string>(config))
    ("F,spec",       "Component specification config file",
     cxxopts::value<std::string>(spec))
    ("T,template",   "Write template options file with current and all default values")
    ("Jacobi",       "Use the standard accurate but slow Jacobi method for SVD computation")
    ("BDCSVD",       "Use the bidiagonal divide and conquer method rather than Jacobi")
    ("allchan",      "use k-means analysis on appended channel data")
    ("distance",     "use distance rather than correlation in w-corr")
    ("H,histo",      "compute the PC contributions to each coefficient series")
    ("C,coefs",      "save time series of coefficients")
    ("reverse",      "verbose output of inverse embedded matrix")
    ("V,version",    "show version")
    ("v,verbose",    "verbose diagnostics to standard output")
    ("X,noCommand",  "do not save command line")
    ("z,zeropad",    "zero pad trajectory matrix (for testing)")
    ("e,ev",         "exit after computing eigenvalues")
    ("triples",      "print all triples for reconstruction")
    ("f,flip",       "png origin in lower left rather than upper left")
    ("Traj",         "use SVD of trajectory matrix rather than covariance matrix")
    ("totVar",       "use total variance for normalization")
    ("totPow",       "use total power for normalization")
    ("writeCov",     "write the covariance matrix (may be large and time consuing)")
    ("h,help",       "this help message");

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("version")) {
    std::cout << PACKAGE_STRING << std::endl;
    return 0;
  }

  // Write YAML template config file and exit
  //
  if (vm.count("template")) {
    SaveConfig(vm, options, "template.yaml");
    return 0;
  }

  if (vm.count("help")) {
    std::cout << std::string(60, '-') << std::endl;
    std::cout << options.help() << std::endl;
    std::cout << std::string(60, '-') << std::endl << std::endl;
    return 0;
  }

  // Read parameters fron the YAML config file
  //
  if (vm.count("config")) {
    try {
      vm = LoadConfig(options, config);
    } catch (cxxopts::OptionException& e) {
      std::cout << "Option error in configuration file: "
		<< e.what() << std::endl;
      return 0;
    }
  }

  bool quiet = true;

  if (vm.count("verbose")) {
    quiet = false;
  }

  bool flip = false;

  if (vm.count("flip")) {
    flip = true;
  }

  bool zeropad = false;
  if (vm.count("zeropad")) {
    zeropad = true;
  }

  if (vm.count("noCommand")==0) {
    std::string cmdFile = prefix + ".cmd_line";
    std::ofstream cmd(cmdFile.c_str());
    if (!cmd) {
      std::cerr << "exp_disk: error opening <" << cmdFile
		<< "> for writing" << std::endl;
    } else {
      cmd << cmd_line << std::endl;
    }

    cmd.close();
  }

  // Detrending style (totVar and totPow are only useful, so far, for
  // noise computation)
  //
  CoefDB::TrendType type;
  if (vm.count("totVar"))
    type = CoefDB::TrendType::totVar;
  else if (vm.count("totPow"))
    type = CoefDB::TrendType::totPow;
  else
    type = CoefDB::TrendType::perChannel;


  // Eigen OpenMP reporting
  //
  if (not quiet)
    std::cout << "Eigen is using " << Eigen::nbThreads()
	      << " threads" << std::endl;

  // Now open and parse the coefficient files
  //
  CoefDB::CoefContainer coefs(spec);

  std::map<CoefDB::Key, std::vector<double> > data;
  std::map<CoefDB::Key, double> mean, var;

  int numT = coefs.times.size();

  // Generate all the channels
  //
  for (auto key : coefs.getKeys()) {

    mean[key] = 0.0;
    var [key] = 0.0;

    for (auto y : coefs.getData(key)) {
      mean[key] += y;
      var [key] += y*y;
      data[key].push_back(y);
    }
  }

  // Print channel info
  //
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

  // Normalize and detrend
  //
  double totVar = 0.0, totPow = 0.0;
  bool useMean = true;

  if (type == CoefDB::TrendType::totPow) {

    // Compute total power from the coefficient database
    //
    for (auto k : coefs.getKeys()) {
      for (auto v : coefs.getData(k)) {
	totPow += v*v;
      }
    }

    // Average total power per time slice
    //
    totPow = std::sqrt(totPow/numT + std::numeric_limits<double>::min());

    if (not quiet) std::cout << "<total Power>=" << totPow << std::endl;

    for (auto & u : mean) {
      CoefDB::Key k = u.first;
      mean[k] /= numT;
      var[k]   = var[k]/numT - mean[k]*mean[k];
      totVar  += var[k];
      var[k]   = sqrt(fabs(var[k]) + std::numeric_limits<double>::min());
    }

    if (vm.count("noMean")) useMean = false;

    for (auto & u : mean) {
      CoefDB::Key k = u.first;
      for (auto & v : data[k]) {
	if (useMean) v -= mean[k];
	v /= totPow;
      }
    }

  } else if (vm.count("totVar")) {

    for (auto & u : mean) {
      CoefDB::Key k = u.first;
      mean[k] /= numT;
      var[k]   = std::max<double>(var[k]/numT - mean[k]*mean[k], 0.0);
      totVar  += var[k];
      var[k]   = std::sqrt(fabs(var[k]) + std::numeric_limits<double>::min());
    }

    totVar = std::sqrt(totVar+std::numeric_limits<double>::min());

    for (auto & u : mean) {
      CoefDB::Key k = u.first;
      for (auto & v : data[k]) {
	v -= mean[k];
	if (totVar>0.0) v /= totVar;
      }
    }
  } else {
    for (auto & u : mean) {
      CoefDB::Key k = u.first;
      mean[k] /= numT;
      var[k]   = std::max<double>(var[k]/numT - mean[k]*mean[k], 0.0);
      var[k]   = std::sqrt(fabs(var[k])+ std::numeric_limits<double>::min());
      for (auto & v : data[k]) {
	v -= mean[k];
	if (var[k]>0.0) v /= var[k];
      }
    }
  }

  std::string filename(prefix + ".data");
  std::ofstream out(filename);

  for (int i=0; i<numT; i++) {
    out << std::setw(6) << coefs.times[i];
    for (auto u : mean) out << std::setw(18) << data[u.first][i];
    out << std::endl;
  }
  out.close();

  int nkeys = mean.size();

  if (vm.count("numW")==0 or numW<=0) numW = numT/2;

  int numK = numT - numW + 1;
  if (zeropad) numK = numT;

  Eigen::MatrixXd Y(numK, numW*nkeys);
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

  if (vm.count("Traj")>0) {
    rank = std::min<int>({static_cast<int>(Y.cols()), static_cast<int>(Y.rows()), npc});
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
  if (vm.count("Traj") == 0 and vm.count("writeCov")) {
    filename = prefix + ".cov";
    out.open(filename);
    out << cov;
    out.close();
  }

  // Store the solution
  // ------------------
  Eigen::VectorXd S;
  Eigen::MatrixXd U;

  // Use one of the built-in Eigen3 algorithms
  //
  if (vm.count("Jacobi")) {
    // -->Using Jacobi
    Eigen::JacobiSVD<Eigen::MatrixXd>
      svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
      
    S = svd.singularValues();	// Get solution
    U = svd.matrixU();
  } else if (vm.count("BDCSVD")) {
    // -->Using BDC
    Eigen::BDCSVD<Eigen::MatrixXd>
      svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    S = svd.singularValues();	// Get solution
    U = svd.matrixU();
  } else if (vm.count("Traj")) {
    // -->Use Random approximation algorithm from Halko, Martinsson,
    //    and Tropp to compute the SVD of the grand trajectory matrix
    auto YY = Y/Scale;

    RedSVD::RedSVD<Eigen::MatrixXd> svd(YY, rank);

    S = svd.singularValues();	// Get solution
    U = svd.matrixV();
  } else {
    // -->Use Random approximation algorithm from Halko, Martinsson,
    //    and Tropp
    RedSVD::RedSVD<Eigen::MatrixXd> svd(cov, rank);

    S = svd.singularValues();	// Get solution
    U = svd.matrixU();
  }


  std::cout << "shape U = " << U.rows() << " x "
	    << U.cols() << std::endl;

  // Rescale the SVD factorization by the Frobenius norm
  //
  S *= Scale;
  if (vm.count("Traj")) {
    for (int i=0; i<S.size(); i++) S(i) = S(i)*S(i)/numK;
  }

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "---- Eigenvalues"                         << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  double tot = 0.0;
  for (int j=0; j<S.size(); j++) tot += S(j);
  std::cout << "Total=" << tot << std::endl;

  double cum = 0.0;
  for (int j=0; j<S.size(); j++) {
    std::cout << std::setw( 5) << j
	      << std::setw(18) << S(j)
	      << std::setw(18) << (cum += S(j))/tot
	      << std::endl;
    if (1.0 - cum/tot < evtol) break;
  }

  filename = prefix + ".ev";
  out.open(filename);
  if (out) {
    cum = 0.0;
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

  if (vm.count("ev")) return 0;

  npc = std::min<int>(npc, numW*nkeys);

  if (not quiet) {
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

  // The eigenvectors {E} are composed of nkeys consecutive segments of
  // of length numW
  //
  filename = prefix + ".evec";
  out.open(filename);
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

  // Compute the PCs by projecting the data
  //
  Eigen::MatrixXd PC = Y * U;

  if (not quiet) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "---- Principal components"                << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    for (int i=0; i<numK; i++) {
      std::cout << std::setw(15) << std::setprecision(6) << coefs.times[i];
      for (int j=0; j<numW*nkeys; j++)
	std::cout << std::setw(15) << std::setprecision(6) << PC(i, j);
      std::cout << std::endl;
    }
  }

  filename = prefix + ".pc";
  out.open(filename);
  if (out) {
    for (int i=0; i<numK; i++) {
      out << std::setw(5) << coefs.times[i];
      for (int j=0; j<npc; j++)
	out << std::setw(15) << PC(i, j);
      out << std::endl;
    }
    out.close();
  } else {
    std::cout << "Could not open <" << filename << ">" << std::endl;
    exit(-1);
  }

  // Reconstructed time series
  //
  int ncomp = std::min<int>({numW, npc, static_cast<int>(PC.cols())});

  std::map<CoefDB::Key, Eigen::MatrixXd> RC;
  for (auto u : mean) RC[u.first].resize(numT, ncomp);

  // Embedded time series matrix
  //
  Eigen::MatrixXd  Z(numT, numW);

  // Split eigenvectors into channel sequences
  //
  std::map<CoefDB::Key, Eigen::MatrixXd> rho;

  int n = 0;
  for (auto u : mean) {
    rho[u.first].resize(numW, rank);
    rho[u.first].fill(0.0);
    for (int i=0; i<numW; i++) {
      for (int j=0; j<ncomp; j++)
	rho[u.first](i, j) = U(numW*n + i, j);
    }
    n++;
  }

  if (not quiet) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "---- Reconstruction"                      << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  for (auto u : mean) {

    for (int w=0; w<ncomp; w++) {

      Z.fill(0.0);

      // Build reverse embedded time series.
      //
      for (int j=0; j<numW; j++) {
	for (int i=0; i<numT; i++)
	  if (i - j >= 0 and i - j < numK) Z(i, j) = PC(i - j, w);
      }

      if (not quiet and vm.count("reverse")) {
	std::cout << "----------------------------------------" << std::endl;
	std::cout << "---- Z for " << u.first << " w=" << w     << std::endl;
	std::cout << "----------------------------------------" << std::endl;

	std::cout << Z << std::endl << std::endl;
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

  if (not quiet) {
    for (int i=0; i<numT; i++) {
      std::cout << std::setw(15) << std::setprecision(6) << coefs.times[i];
      for (auto u : mean) {
	double acc = 0.0;
	for (int j=0; j<ncomp; j++) acc += RC[u.first](i, j);
	std::cout << std::setw(15) << std::setprecision(6) << acc;
      }
      std::cout << std::endl;
    }
  }

  for (auto u : mean) {
    std::ostringstream sout;
    sout << prefix;
    for (auto v : u.first) sout << "_" << v;
    sout << ".elem";
    std::ofstream out(sout.str());
    if (out) {
      double disp = totVar;
      if (type == CoefDB::TrendType::totPow) disp = totPow;
      if (disp==0.0) disp = var[u.first];
      
      for (int i=0; i<numT; i++) {
	out << std::setw(15) << std::setprecision(6) << coefs.times[i];
	for (int w=0; w<ncomp; w++)
	  out << std::setw(15) << std::setprecision(6)
	      << RC[u.first](i, w)*disp + u.second;
	out << std::endl;
      }
    } else {
      std::cout << "Could not open <" << sout.str() << ">" << std::endl;
      exit(-1);
    }
  }

#ifdef HAVE_LIBPNGPP
  if (vm.count("histo")) {

    if (not quiet) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- Elemental fractions"                 << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }

    filename = prefix + ".f_contrib";
    out.open(filename);

    // Used in both computations
    //
    std::map<CoefDB::Key, std::vector<double> > values;

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
      std::map<CoefDB::Key, double> norm;

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

  }
#endif

  if (not quiet) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "---- Coefficient periodogram"             << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  {
    double DT = coefs.times[1] - coefs.times[0];

    int nfreq = numT/2+1;

    Eigen::MatrixXd pw(numW, nfreq);
    Eigen::VectorXd p0(nfreq), in(numT), F, P;
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
      out.open(filename.str());
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

  if (not quiet) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "---- Principal component periodogram"     << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  {
    double DT = coefs.times[1] - coefs.times[0];

    int nfreq = numK/2 + 1;

    Eigen::MatrixXd pw(nfreq, npc);
    Eigen::VectorXd in(numK), F, P;
    
    pw.setZero();

    for (int j=0; j<npc; j++) {
      for (int i=0; i<numK; i++) in(i) = PC(i, j);
      TransformFFT fft(DT, in);
      fft.Power(F, P);
      for (int i=0; i<nfreq; i++) pw(i, j) = P(i);
    }

    std::ostringstream filename;
    filename << prefix << ".pc_power";
    out.open(filename.str());
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

  if (not quiet) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "---- Coefficient periodogram"             << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  {
    double DT = coefs.times[1] - coefs.times[0];

    int nfreq = numT/2 + 1;

    Eigen::MatrixXd pw(numW, nfreq);
    Eigen::VectorXd p0(nfreq), in(numT), F, P;
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
      out.open(filename.str());
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

#ifdef HAVE_LIBPNGPP
  if (not quiet) {
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
    if (vm.count("distance")) use_dist = true;

    for (auto u : mean) {
      Eigen::MatrixXd wc = wCorr(RC[u.first]);

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
      image.write(sout.str());
    }

    {
      Eigen::MatrixXd wc = wCorrAll(RC);

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
      image.write(sout.str());
    }

  }
#endif

  if (vm.count("kmeans")) {

    if (not quiet) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "---- K-means group analysis"              << std::endl;
      std::cout << "---- numT=" << numT << " numW=" << numW   << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }

    filename = prefix + ".kmeans";
    out.open(filename);
    if (out) {
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
	out << std::string(60, '-') << std::endl
	    << " *** n=" << u.first << std::endl
	    << std::string(60, '-') << std::endl;

	for (int j=0; j<results.size(); j++) {
	  out << std::setw(6)  << j
	      << std::setw(12) << std::get<1>(results[j])
	      << std::endl;
	}

	// Write to stdout
	//
	std::cout << std::string(60, '-') << std::endl
		  << " *** n=" << u.first << std::endl
		  << std::string(60, '-') << std::endl;

	for (int j=0; j<results.size(); j++) {
	  std::cout << std::setw(6)  << j
		    << std::setw(12) << std::get<1>(results[j])
		    << std::endl;
	}
      }

      if (vm.count("allchan")) {

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
	out << std::string(60, '-') << std::endl
	    << " *** total"         << std::endl
	    << std::string(60, '-') << std::endl;

	for (int j=0; j<results.size(); j++) {
	  out << std::setw(6) << j
	      << std::setw(9) << std::get<1>(results[j])
	      << std::endl;
	}

	// Write to stdout
	//
	std::cout << std::string(60, '-') << std::endl
		  << " *** total"         << std::endl
		  << std::string(60, '-') << std::endl;

	for (int j=0; j<results.size(); j++) {
	  std::cout << std::setw(6) << j
		    << std::setw(9) << std::get<1>(results[j])
		    << std::endl;
	}

      }

      out.close();

    } else {
      std::cerr << "Error opening file <" << filename << ">" << std::endl;
    }
  }

  filename = prefix + ".recon";
  out.open(filename);
  if (out) {
    for (int i=0; i<numT; i++) {
      out << std::setw(15) << std::setprecision(6) << coefs.times[i];
      for (auto u : mean) {
	double disp = totVar;
	if (type == CoefDB::TrendType::totPow) disp = totPow;
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
      out << std::setw(15) << std::setprecision(6) << coefs.times[i];
      for (auto k : mean) {
	double disp = totVar;
	if (type == CoefDB::TrendType::totPow) disp = totPow;
	if (disp==0.0) disp = var[k.first];

	double acc = 0.0;
	for (int j=0; j<ncomp; j++) acc += RC[k.first](i, j);
	out << std::setw(15) << std::setprecision(6) << acc*disp + k.second
	    << std::setw(15) << std::setprecision(6) << data[k.first][i]*disp + k.second;
      }
      out << std::endl;
      n++;
    }
    out.close();
  } else {
    std::cout << "Could not open <" << filename << ">" << std::endl;
    exit(-1);
  }

  if (vm.count("triples")) {

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
	if (type == CoefDB::TrendType::totPow) disp = totPow;
	if (disp==0.0) disp = var[k];

	out1 << "# " << u.first << std::endl;
	out2 << "# " << u.first << std::endl;
	for (int i=0; i<numT; i++) {
	  out1 << std::setw(15) << std::setprecision(6) << coefs.times[i];
	  out2 << std::setw(15) << std::setprecision(6) << coefs.times[i];
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
      n++;
    }
  }

  std::vector<std::set<int>> groups;

  if (vm.count("group")) {
    std::ifstream gf(grpfile);
    // Iterator the file, line by line
    //
    if (gf.good()) {
      std::string line;
      while (gf.good()) {
	// Get the line from the file
	//
	std::getline(gf, line);
	// Make a stream
	//
	std::istringstream str_in(line);
	// Make a vector of ints and add to group list
	//
	groups.push_back(std::set<int>((std::istream_iterator<int>(str_in)),
				       std::istream_iterator<int>()));
      }
    }
  } else {
    for (int n=0; n<ncomp; n++) groups.push_back({n});
  }

  int ngroups = groups.size();

  if (ngroups) {
    coefs.write_recon(ncomp, type, useMean, totPow, totVar,
		      mean, var, RC, groups);
  }

  if (vm.count("coefs")) {
    coefs.write_ascii(ncomp, type, useMean, totPow, totVar,
		      mean, var, RC, groups);
  }


  return 0;
}
