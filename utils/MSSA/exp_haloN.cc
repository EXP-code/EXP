//
// M-SSA for EXP [spherical] coefficients
//
// This version can combine multiple simulations.
//
// Updated to use fixed-rank approximation for the SVD to save time
// and space.  Uses the approximate SVD randomized algorithm from
// Halko, Martinsson, and Tropp by default.  Use BDCSVD or Jacobo
// flags for the standard methods.
//

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <regex>
#include <cmath>
#include <set>
#include <map>

#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>

#include <TransformFFT.H>
#include <ColorGradient.H>

#include "Coefs.H"
#include "RedSVD.H"

#include <config_exp.h>
#include <YamlConfig.H>


const double PI2() { return std::atan(1)*8; }

typedef std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned> Key;

// For reading multiple coefficient files
//
struct CoefData {
  std::vector<double> times;	// Coefficient times
  std::map<Key, std::vector<double>> coefs; // Coefficients
  std::string file;			    // Name of the outcoef file
  int index;			// This is redundant, but makes for
				// easy iteration
};

// Save input data as shared-pointer stanza to remove the copying burden
//
using CoefDataPtr = std::shared_ptr<CoefData>;

// Coefficient data store
//
std::vector<CoefDataPtr> Coefs;

// Helper for key output
//
template<typename T, typename U, typename V, typename W, typename X>
std::ostream& operator<< (std::ostream& out, const std::tuple<T, U, V, W, X>& t)
{
  const char lab[] = {'c', 's'};
  std::ostringstream sout;
  sout << '(' << std::get<0>(t)
       << ',' << std::get<1>(t)
       << ',' << std::get<2>(t)
       << ',' << lab[std::get<3>(t)]
       << ',' << std::get<4>(t)
       << ')';
  out << sout.str();
  return out;
}


typedef std::tuple<unsigned, unsigned> LMkey;

//! Read a single set of coefficients from a file, add to global data
//! base
std::map<std::string, int>
spherical_read(const std::string& file, unsigned stride,
	       double tmin=-std::numeric_limits<double>::max(),
	       double tmax= std::numeric_limits<double>::max())
{
  // Make data stanza
  //
  auto cp = std::make_shared<CoefData>();
  cp->index = Coefs.size();
  cp->file  = file;

  std::map<std::string, int> ret;
  std::ifstream in(file);

  std::map<double, SphCoefsPtr> data;
  unsigned counter = 0;

  while (in.good()) {
    SphCoefsPtr c = std::make_shared<SphCoefs>();
    if (not c->read(in)) break;
    if (c->header.tnow >= tmin and c->header.tnow <= tmax) {
      if (counter++ % stride == 0) data[c->header.tnow] = c;
    }
  }

  int lmax   = data.begin()->second->header.Lmax;
  int nmax   = data.begin()->second->header.nmax;
  int ntimes = data.size();

  for (auto v : data) {
    cp->times.push_back(v.second->header.tnow);
  }

  for (int l=0; l<=lmax; l++) {
    for (int m=0; m<=l; m++) {
      for (int n=0; n<nmax; n++) {
	LMkey lmk  = {l, m};
	Key   key0 = {l, m, n, 0, cp->index};
	Key   key1 = {l, m, n, 1, cp->index};
	cp->coefs[key0].resize(ntimes);
	if (m) cp->coefs[key1].resize(ntimes);
	for (int i=0; i<cp->times.size(); i++) {
	  cp->coefs[key0][i] = data[cp->times[i]]->cos_c[lmk][n];
	  if (m) cp->coefs[key1][i] = data[cp->times[i]]->sin_c[lmk][n];
	}
      }
    }
  }

  Coefs.push_back(cp);		// Add to global data

  ret["lmax"]   = lmax;
  ret["nmax"]   = nmax;
  ret["ntimes"] = ntimes;

  return ret;
}


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



int
main(int argc, char **argv)
{
  //--------------------------------------------------
  // Parameters
  //--------------------------------------------------
  std::vector<std::string> dfiles;
  std::string prefix, keyfile, grpfile, config;
  int numT, numW, L, M, nmin, nmax, npc, stride, skip, lmin, lmax;
  double evtol, tmin, tmax;

  //--------------------------------------------------
  // Command-line parsing
  //--------------------------------------------------

  std::string cmd_line;
  for (int i=0; i<argc; i++) {
    cmd_line += argv[i];
    cmd_line += " ";
  }

  cxxopts::Options options("exp_halo","This routine uses M-SSA to analyze basis coefficient\n\nOptions");

  options.add_options()
    ("o,output",     "output file prefix",
     cxxopts::value<std::string>(prefix)->default_value("exp_mssa"))
    ("d,files",      "EXP coefficient data file(s)",
     cxxopts::value<std::vector<std::string>>(dfiles))
    ("k,keyfile",    "list of keys to process",
     cxxopts::value<std::string>(keyfile))
    ("W,numW",       "window size",
     cxxopts::value<int>(numW)->default_value("10"))
    ("lmin",         "spherical harmonic lower limit",
    cxxopts::value<int>(lmin)->default_value("-1"))
    ("lmax",         "spherical harmonic upper limit",
    cxxopts::value<int>(lmax)->default_value("-1"))
    ("L,harmonicL",  "spherical harmonic l",
     cxxopts::value<int>(L)->default_value("0"))
    ("M,harmonicM",  "spherical harmonic m",
     cxxopts::value<int>(M)->default_value("0"))
    ("n,nmin",       "minimum radial order",
     cxxopts::value<int>(nmin)->default_value("0"))
    ("N,nmax",       "maximum radial order (-1 means use coefficient maximum)",
     cxxopts::value<int>(nmax)->default_value("-1"))
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
    ("tmin",         "minimum simulation snapshot time for processing",
     cxxopts::value<double>(tmin)->default_value("-1.0e42"))
    ("tmax",         "maximum simulation snapshot time for processing",
     cxxopts::value<double>(tmax)->default_value("1.0e42"))
    ("c,config",      "Input parameter config file",
     cxxopts::value<std::string>(config))
    ("T,template",   "Write template options file with current and all default values")
    ("Jacobi",       "Use the standard accurate but slow Jacobi method for SVD cmputation")
    ("BDCSVD",       "Use the bidiagonal divide and conquer method rather than Jacobi")
    ("Traj",         "Use RedSVD to decompose the trajectory matrix")
    ("H,histo",      "compute the PC contributions to each coefficient series")
    ("C,coefs",      "save time series of coefficients")
    ("reverse",      "verbose output of inverse embedded matrix")
    ("v,version",    "show version")
    ("X,noCommand",  "do not save command line")
    ("z,zeropad",    "zero pad trajectory matrix (for testing)")
    ("e,ev",         "exit after computing eigenvalues")
    ("triples",      "print all triples for reconstruction")
    ("channels",     "print principle components by channels")
    ("q,quiet",      "only print ev values to stdout")
    ("p,pairs",      "group eigenvalues by pairs")
    ("f,flip",       "png origin in lower left rather than upper left")
    ("Lfull",        "use all m value for the specified l")
    ("totVar",       "use total variance for normalization")
    ("compare",      "print reconstruction comparison for each simulation")
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

  // Print help info
  //
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

  bool quiet = false;

  if (vm.count("quiet")) {
    quiet = true;
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
      std::cerr << "test_exp: error opening <" << cmdFile
		<< "> for writing" << std::endl;
    } else {
      cmd << cmd_line << std::endl;
    }

    cmd.close();
  }

  // Sanity check on data files: make we have at least one and sure
  // each one exists
  //
  if (vm.count("files")==0) {
    std::cerr << argv[0]
	      << ": **error** you must specify at least one data file using --files"
	      << std::endl;
    exit(-1);
  }

  // Check for existence of files
  //
  bool okay = true;
  for (auto file : dfiles) {
    if (not std::filesystem::exists(file)) {
      std::cerr << argv[0]
		<< ": **error** file <" << file << "> does not exist?"
		<< std::endl;
      okay = false;
    }
  }
  if (not okay) {
    exit(-1);
  }

  // Attempt to read data from all data files
  //
  int NMAX, LMAX, NTIM, count = 0;

  for (auto file : dfiles) {

    auto param = spherical_read(file, skip, tmin, tmax);

    if (count==0) {
      NMAX = param["nmax"];
      LMAX = param["lmax"];
      NTIM = param["ntimes"];
    } else {
      if (NMAX != param["nmax"]) {
	std::cerr << "Coefficient file mismatch for file <" << file << "> : "
		  << "found Nmax=" << param["nmax"] << " but expected " << NMAX
		  << std::endl;
	exit(-2);
      }
      if (LMAX != param["lmax"]) {
	std::cerr << "Coefficient file mismatch for file <" << file << "> : "
		  << "found Lmax=" << param["lmax"] << " but expected " << LMAX
		  << std::endl;
	exit(-2);
      }
      if (NTIM != param["ntimes"]) {
	std::cerr << "Coefficient file mismatch for file <" << file << "> : "
		  << "found Ntim=" << param["ntimes"] << " but expected " << NTIM
		  << std::endl;
	exit(-2);
      }
    }

    numT = Coefs.back()->times.size();

    assert(NTIM == numT);		// Sanity check
  }

  std::cout << "Number of times: " << numT << std::endl;

  std::vector<double> times = Coefs.front()->times;
  std::cout << "T in [" << times.front() << ", " << times.back() << "]"
	    << std::endl;


				// Make time series
  std::map<Key, std::vector<double> > data;
  std::map<Key, double> mean, var;
  int ndim = 0;			// Data dimension

  std::set<std::pair<unsigned, unsigned>> LMset;

  if (vm.count("keyfile")) {

    std::ifstream in(keyfile);

    while (in.good()) {
      int cs, l, m, n;
      std::string line;
      std::getline(in, line);
				// Forget comments
      if (line[0] == '#') continue;
				// Truncate comment
      if (line.find('#') != std::string::npos)
	line = line.substr(line.find('#'));

      std::istringstream sin(line);
      sin >> cs;
      sin >> l;
      sin >> m;
      sin >> n;

      // Compute maximum n value from the key list
      //
      nmax = std::max<int>(nmax, n);

      // Make keys for every file
      //
      for (auto cf : Coefs) {
	Key key = {l, m, n, cs, cf->index};

	if (cf->coefs.find(key) == cf->coefs.end()) {
	  std::cout << "Key " << key << " not found in coefs, keys are" << std::endl;
	  for (auto v : cf->coefs) std::cout << v.first << std::endl;
	  exit(-1);
	}

	LMset.insert({l, m});

	mean[key] = 0.0;
	var [key] = 0.0;

	for (int i=0; i<numT; i++) {
	  double y = cf->coefs[key][i];
	  mean[key] += y;
	  var [key] += y*y;
	  data[key].push_back(y);
	}
      }
    }

    ndim = mean.size();

  }
  // No specified key file
  //
  else {

    if (L>LMAX) {
      std::cout << "L=" << L << " is larger than LMAX=" << LMAX
		<< " . . . quitting" << std::endl;
      exit(1);
    }

    if (M>L) {
      std::cout << "M=" << M << " is larger than L=" << L
		<< " . . . quitting" << std::endl;
      exit(1);
    }

    if (nmax<0) nmax = NMAX;
    nmax = std::min<double>(NMAX, nmax);

    if (vm.count("group")) {
      nmin = 0;
    } else {
      ndim = nmax - nmin + 1;
    }

    if (vm.count("group")) nmin = 0;

    ndim = 0;

    if (lmin>=0 and lmax>=0) {

      for (int L=lmin; L<=std::min<int>(lmax, LMAX); L++) {

	int Mmin = M, Mmax = M;
	if (vm.count("Lfull")) {
	  Mmin = 0;
	  Mmax = L;
	}

	for (int mm=Mmin; mm<=Mmax; mm++) {

	  LMset.insert({L, mm});

	  for (int n=nmin; n<nmax; n++) {

	    // File loop
	    //
	    for (auto cf : Coefs) {

	      int cmax = 2;		// Default: cosines and sines
	      if (mm==0) cmax = 1;	// Cosines only of m=0
	      for (int c=0; c<cmax; c++) {
		Key key = {L, mm, n, c, cf->index};
		if (cf->coefs.find(key) == cf->coefs.end()) {
		  std::cout << "Key " << key << " not found in coefs, keys are" << std::endl;
		  for (auto v : cf->coefs) std::cout << v.first << std::endl;
		  exit(-1);
		}

		mean[key] = 0.0;
		var [key] = 0.0;

		for (int i=0; i<numT; i++) {
		  double y = cf->coefs[key][i];
		  mean[key] += y;
		  var [key] += y*y;
		  data[key].push_back(y);
		}

		ndim++;
	      }
	      // Cosine/sine loop
	    }
	    // File loop
	  }
	  // N loop
	}
	// M loop
      }

    } else {

      int Mmin = M, Mmax = M;
      if (vm.count("Lfull")) {
	Mmin = 0;
	Mmax = L;
      }

      for (int mm=Mmin; mm<=Mmax; mm++) {

	LMset.insert({L, mm});

	for (int n=nmin; n<nmax; n++) {

	  // File loop
	  //
	  for (auto cf : Coefs) {

	    int cmax = 2;		// Default: cosines and sines
	    if (mm==0) cmax = 1;	// Cosines only if m=0
	    for (int c=0; c<cmax; c++) {
	      Key key = {L, mm, n, c, cf->index};
	      if (cf->coefs.find(key) == cf->coefs.end()) {
		std::cout << "Key " << key << " not found in coefs, keys are" << std::endl;
		for (auto v : cf->coefs) std::cout << v.first << std::endl;
		exit(-1);
	      }

	      mean[key] = 0.0;
	      var [key] = 0.0;

	      for (int i=0; i<numT; i++) {
		double y = cf->coefs[key][i];
		mean[key] += y;
		var [key] += y*y;
		data[key].push_back(y);
	      }
	      
	      ndim++;
	    }
	    // Cosine/sine loop
	  }
	  // File loop
	}
	// N loop
      }
      // M loop
    }
  }
  // Subspace

  if (not quiet) {
    std::cout << "Using coefficient keys (L, M, n, cs, i):";
    for (auto v : mean) std::cout << " " << v.first;
    std::cout << std::endl;
  }

  // Normalize
  //
  double totVar = 0.0;

  if (vm.count("totVar")) {
    for (auto & u : mean) {
      Key k = u.first;
      mean[k] /= numT;
      var[k]   = var[k]/numT - mean[k]*mean[k];
      totVar  += var[k];
    }

    for (auto & u : mean) {
      Key k = u.first;
      for (auto & v : data[k]) {
	v -= mean[k];
	v /= sqrt(totVar);
      }
    }
  } else {
    for (auto & u : mean) {
      Key k = u.first;
      mean[k] /= numT;
      var[k]   = var[k]/numT - mean[k]*mean[k];
      for (auto & v : data[k]) {
	v -= mean[k];
	v /= sqrt(var[k]);
      }
    }
  }

  std::string filename(prefix + ".data");
  std::ofstream out(filename);

  for (int i=0; i<numT; i++) {
    out << std::setw(6) << times[i];
    for (auto u : mean) out << std::setw(18) << data[u.first][i];
    out << std::endl;
  }
  out.close();


  int numK = numT - numW + 1;
  if (zeropad) numK = numT;

  Eigen::MatrixXd Y(numK, numW*ndim);
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
  int rank = npc;

  if (vm.count("Traj")>0) {
    rank = std::min<int>({static_cast<int>(Y.cols()), static_cast<int>(Y.rows()), npc});
    Scale = Y.norm();
    if (Scale<=0.0) {
      std::cout << "Frobenius norm of trajectory matrix is <= 0!" << std::endl;
      exit(-1);
    }
  } else {
    cov   = Y.transpose() * Y/numK;
    Scale = cov.norm();

    filename = prefix + ".cov";
    out.open(filename);
    out << cov;
    out.close();

    if (Scale<=0.0) {
      std::cout << "Frobenius norm of covariance matrix is <= 0!" << std::endl;
      exit(-1);
    } else {
      cov /= Scale;
    }
    
  }

  std::cout << "Eigen is using " << Eigen::nbThreads()
	    << " threads" << std::endl;

  // Store the solution
  // ------------------
  Eigen::VectorXd S;
  Eigen::MatrixXd U;

  // Use one of the built-in Eigen3 algorithms
  //
  if (vm.count("Jacobi")) {
    // -->Using Jacobi
    std::cout << "Using JacobiSVD" << std::endl;

    Eigen::JacobiSVD<Eigen::MatrixXd>
      svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    S = svd.singularValues();	// Get solution
    U = svd.matrixU();
  } else if (vm.count("BDCSVD")) {
    // -->Using BDC
    std::cout << "Using BDCSVD" << std::endl;

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
    std::cout << "Using RedSVD" << std::endl;

    RedSVD::RedSVD<Eigen::MatrixXd> svd(cov, rank);

    S = svd.singularValues();	// Get solution
    U = svd.matrixU();
  }

  

  std::cout << "shape U = " << U.rows() << " x "
	    << U.cols() << std::endl;

  // Rescale the SVD factorization by the Frobenius norm of the input
  // matrix
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

  npc = std::min<int>(npc, numW*ndim);

  if (not quiet) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "---- Eigenvectors"                        << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    for (int j=0; j<std::min<int>(numW*ndim, rank); j++) {
      std::cout << std::setw(5) << j;
      for (int i=0; i<numW*ndim; i++)
	std::cout << std::setw(15) << std::setprecision(6) << U(i, j);
      std::cout << std::endl;
    }
  }

  int ncomp = std::min<int>({numW, npc});
  
  // The eigenvectors {E} are composed of ndim consecutive segments of
  // of length numW

  filename = prefix + ".evec";
  out.open(filename);
  if (out) {
    for (int j=0; j<std::min<int>(numW*ndim, rank); j++) {
      out << std::setw(5) << j;
      for (int i=0; i<numW*ndim; i++) out << std::setw(15) << U(i, j);
      out << std::endl;
    }
    out.close();
  } else {
    std::cout << "Could not open <" << filename << ">" << std::endl;
    exit(-1);
  }

  Eigen::MatrixXd PC = Y * U;

  if (vm.count("channels")) {
    std::vector<Eigen::MatrixXd> PCN(ndim);
    for (int n=0; n<ndim; n++) {
      Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(numW*ndim, numW*ndim);
      for (int i=0; i<numW; i++) D1(n*numW+i, n*numW+i) = 1.0;
      PCN[n] = Y * D1 * U;
      std::cout << "PCN[" << n << "]: "
		<< PCN[n].rows() << " x " << PCN[n].cols() << std::endl;
    }

    for (int k=0; k<npc; k++) {
      std::ostringstream fname;
      fname << prefix << ".pc." << k;
      out.open(fname.str());
      if (out) {
	auto it = mean.begin();	// Header
	for (int n=0; n<ndim; n++, it++)
	  out << "# " << std::setw(5) << n+1 << it->first << std::endl;
				// Data
	for (int i=0; i<numK; i++) {
	  out << std::setw(5) << times[i];
	  for (int n=0; n<ndim; n++) {
	    out << std::setw(15) << PCN[n](i, k);
	  }
	  out << std::endl;
	}
	out.close();
      } else {
	std::cout << "Could not open <" << fname.str() << ">" << std::endl;
	exit(-1);
      }
    }
  }

  if (not quiet) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "---- Principal components"                << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    for (int i=0; i<numK; i++) {
      std::cout << std::setw(15) << std::setprecision(6) << times[i];
      for (int j=0; j<ncomp; j++)
	std::cout << std::setw(15) << std::setprecision(6) << PC(i, j);
      std::cout << std::endl;
    }
  }

  filename = prefix + ".pc";
  out.open(filename);
  if (out) {
    for (int i=0; i<numK; i++) {
      out << std::setw(5) << times[i];
      for (int j=0; j<ncomp; j++)
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
  std::map<Key, Eigen::MatrixXd> RC;

  for (auto u : mean) RC[u.first].resize(numT, ncomp);

  // Embedded time series matrix
  //
  Eigen::MatrixXd  Z(numT, numW);

  // Split eigenvectors into channel sequences
  //
  std::map<Key, Eigen::MatrixXd> rho;

  int n = 0;
  for (auto u : mean) {
    rho[u.first].resize(numW, npc);
    rho[u.first].fill(0.0);
    for (int i=0; i<numW; i++) {
      for (int j=0; j<npc; j++)
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
      std::cout << std::setw(15) << std::setprecision(6) << times[i];
      for (auto u : mean) {
	double acc = 0.0;
	for (int j=0; j<ncomp; j++) acc += RC[u.first](i, j);
	std::cout << std::setw(15) << std::setprecision(6) << acc;
      }
      std::cout << std::endl;
    }
  }

  filename = prefix + ".elem";
  out.open(filename);
  if (out) {
    for (int i=0; i<numT; i++) {
      out << std::setw(15) << std::setprecision(6) << times[i];
      for (auto u : mean) {
	double acc = 0.0;
	for (int j=0; j<ncomp; j++)
	  out << std::setw(15) << std::setprecision(6) << RC[u.first](i, j);
      }
      out << std::endl;
    }
    out.close();
  } else {
    std::cout << "Could not open <" << filename << ">" << std::endl;
    exit(-1);
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
    std::map<Key, std::vector<double> > values;

    // Set up for PNG images
    //
    const int minSize = 600;
    int ndupX = 1, ndupY = 1;
    if (ncomp < minSize)       ndupX = minSize/ncomp + 1;
    if (mean.size() < minSize) ndupY = minSize/mean.size() + 1;

    std::cout << "Matrix size: " << ncomp << " X " << mean.size() << std::endl;

    int ysiz = mean.size()*ndupY;
    int xsiz = ncomp*ndupX;

    png::image< png::rgb_pixel > image1(xsiz, ysiz);
    png::image< png::rgb_pixel > image2(xsiz, ysiz);

    ColorGradient color;
    color.createFiveColorHeatMapGradient();

    if (out) {

      // Column header
      //
      out << std::setw(4) << "# m"
	  << std::setw(4) << "n"
	  << std::setw(4) << "cs";
      for (int j=0; j<ncomp; j++) {
	std::ostringstream sout; sout << "PC" << j;
	out << std::setw(18) << sout.str();
      }
      out << std::endl;

      // For each series, compute ||a^k_j||
      //
      std::map<Key, double> norm;

      for (auto u : mean) {
	auto key = u.first;
	values[key].resize(ncomp, 0.0);

	for (int i=0; i<numT; i++) {
	  for (int j=0; j<ncomp; j++) {
	    values[key][j] += RC[key](i, j)*RC[key](i, j);
	  }
	}
      }

      // This is norm for each series over the entire reconstruction
      // (that is, fixed coefficient channel, summed over all PCs)
      for (auto u : mean) {
	auto key = u.first;
	norm[key] = 0.0;
	for (int j=0; j<ncomp; j++) {
	  norm[key] += values[key][j];
	}
      }

      double maxVal = 0.0;
      for (auto u : mean) {
	auto key = u.first;
	out << std::setw(4) << std::get<0>(key)
	    << std::setw(4) << std::get<1>(key)
	    << std::setw(4) << std::get<2>(key);

	for (int j=0; j<ncomp; j++) {
	  out << std::setw(18) << sqrt(values[key][j]/norm[key]);
	  maxVal = std::max<double>(maxVal, values[key][j]/norm[key]);
	}
	out << std::endl;
      }
      out.close();

      int i = 0;
      for (auto u : values) {
	auto key = u.first;

	// Sanity check
	{
	  if (values.find(key) == values.end()) {
	    std::cout << "Could not find key=" << key << " in values"
		      << std::endl;
	  }
	  if (norm.find(key) == norm.end()) {
	    std::cout << "Could not find key=" << key << " in norm"
		      << std::endl;
	  }
	}

	for (int j=0; j<ncomp; j++) {
	  // Sanity check
	  {
	    if (j >= values[key].size()) {
	      std::cout << "Could not find j=" << j << " in values[key] "
			<< " which has size=" << values[key].size()
			<< std::endl;
	    }

	  }

	  png::rgb_pixel cval = color(sqrt(values[key][j]/norm[key]/maxVal));
	  for (size_t yy = i*ndupY; yy < (i+1)*ndupY; yy++) {
	    for (size_t xx = j*ndupX; xx < (j+1)*ndupX; xx++) {
	      if (xx<xsiz and yy<ysiz) {
		try {
		  image1.set_pixel(xx, yy, cval);
		}
		catch (png::error& error) {
		  std::cerr << "PNG error: " << error.what() << std::endl
			    << "PNG error: xx, yy=" << xx << ", " << yy
			    << "[" << xsiz << ", " << ysiz << "]" << std::endl;
		}
	      } else {
		std::cout << "Image bounds: xx, yy=" << xx << ", " << yy
			  << "[" << xsiz << ", " << ysiz << "]" << std::endl;
	      }
	    }
	  }
	}
	i++;
      }

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
	for (int j=0; j<ncomp; j++) {
	  norm[j] += values[key][j];
	}
      }

      double maxVal = 0.0;
      for (int j=0; j<ncomp; j++) {
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

      for (int j=0; j<ncomp; j++) {
	int i=0;
	for (auto u : mean) {
	  auto key = u.first;
	  png::rgb_pixel cval = color(sqrt(values[key][j]/norm[j]/maxVal));
	  for (size_t yy = i*ndupY; yy < (i+1)*ndupY; yy++) {
	    for (size_t xx = j*ndupX; xx < (j+1)*ndupX; xx++) {
	      try {
		image2.set_pixel(xx, yy, cval);

	      } catch (png::error& error) {
		std::cerr << "PNG error: " << error.what() << std::endl
			  << "PNG error: xx, yy=" << xx << ", " << yy
			  << "[" << xsiz << ", " << ysiz << "]" << std::endl;
	      }
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
    std::cout << "---- Periodogram"                         << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  {
    int nfreq = numT/2;
    if (numT % 2 == 0) nfreq++;

    Eigen::MatrixXd pw(numW, nfreq);
    Eigen::VectorXd p0(nfreq);
    Eigen::VectorXd in(numT), F, P;

    for (auto u : mean) {

      for (int i=0; i<numT; i++) {
	in(i) = 0.0;
	for (int j=0; j<ncomp; j++) in(i) += RC[u.first](i, j);
	TransformFFT fft(times[1] - times[0], in);
	fft.Power(F, p0);
      }

      for (int j=0; j<ncomp; j++) {
	for (int i=0; i<numT; i++) in(i) = RC[u.first](i, j);

	TransformFFT fft(times[1] - times[0], in);
	fft.Power(F, P);
	for (int k=0; k<nfreq; k++) pw(j, k) = P(k);
      }

      std::ostringstream filename;
      filename << prefix << ".power_" << u.first;
      out.open(filename.str());
      if (out) {
	out << "# " << u.first << std::endl;
	out << "# " << std::setw(13) << "Freq"
	    << std::setw(15) << "Period"
	    << std::setw(15) << "Full";
	for (int j=0; j<nfreq; j++) {
	  std::ostringstream sout; sout << "PC " << j;
	  out << std::setw(15) << sout.str();
	}
	out << "# " << std::setw(13) << "[1]"
	    << std::setw(15) << "[2]"
	    << std::setw(15) << "[3]";
	for (int j=0; j<nfreq; j++) {
	  std::ostringstream sout; sout << '[' << j+4 << ']';
	  out << std::setw(15) << sout.str();
	}
	out << std::endl;

	for (int j=0; j<nfreq; j++) {
	  out << std::setw(15) << std::setprecision(6) << F(j)
	      << std::setw(15) << std::setprecision(6) << 2.0*M_PI/F(j)
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
    int ndup = 1, ncnt = 0;
    int nDim = std::min<int>(numW, npc);
    if (numW < minSize) ndup = minSize/nDim + 1;
    if (ndup > 1) std::cout << "Pixel duplication=" << ndup << std::endl;

    for (auto u : mean) {
      Eigen::MatrixXd wc = wCorr(RC[u.first]);

      png::image< png::rgb_pixel > image(nDim*ndup, nDim*ndup);
      ColorGradient color;

      for (size_t y = 0; y < nDim; y++) {
	for (size_t x = 0; x < nDim; x++) {
	  png::rgb_pixel cval = color(wc(y, x));
	  for (size_t yy = y*ndup; yy < (y+1)*ndup; yy++) {
	    for (size_t xx = x*ndup; xx < (x+1)*ndup; xx++) {
	      try {
		if (flip) image.set_pixel(xx, image.get_height()-1-yy, cval);
		else      image.set_pixel(xx, yy, cval);
	      }
	      catch (png::error& error) {
		std::cerr << "PNG error: " << error.what() << std::endl;
	      }
	    }
	  }
	}
      }

      std::ostringstream sout; sout << prefix + "_" << ncnt++ << ".png";
      image.write(sout.str());
    }
  }
#endif


  filename = prefix + ".recon";
  out.open(filename);
  if (out) {
    const int wid = 18;		// Field width

    // Header #1
    //
    out << "#" << std::setw(wid-1) << std::right << "Time";
    for (auto k : mean) {
      std::ostringstream sout; sout << k.first << "(r)";
      out << std::setw(wid) << std::right << sout.str();
    }
    out << std::endl;

    // Header #2
    //
    out << "#" << std::setw(wid-8) << std::right << "[1] ";
    int hcnt = 2;
    for (auto k : mean) {
      std::ostringstream sout; sout << "[" << hcnt++ << "] ";
      out << std::setw(wid) << std::right << sout.str();
    }
    out << std::endl;

    for (int i=0; i<numT; i++) {
      out << std::setw(wid) << std::setprecision(6) << times[i];
      for (auto u : mean) {
	auto k = u.first;
	double acc = 0.0;
	for (int j=0; j<ncomp; j++) acc += RC[k](i, j);

	// Retrend
	//
	if (vm.count("totVar"))
	  acc = acc*sqrt(totVar) + mean[k];
	else
	  acc = acc*sqrt(var[k]) + mean[k];

	out << std::setw(wid) << std::setprecision(6) << acc;
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
    const int wid = 18;		// Field width

    // Header #1
    //
    out << "#" << std::setw(wid-1) << std::right << "Time";
    for (auto k : mean) {
      std::ostringstream sout; sout << k.first << "(r)";
      out << std::setw(wid) << std::right << sout.str();
      sout.str(""); sout << k.first << "(d)";
      out << std::setw(wid) << std::right << sout.str();
    }
    out << std::endl;

    // Header #2
    //
    out << "#" << std::setw(14) << std::right << "[1] ";
    int hcnt = 2;
    for (auto k : mean) {
      std::ostringstream sout; sout << "[" << hcnt++ << "] ";
      out << std::setw(wid) << std::right << sout.str();
      sout.str(""); sout << "[" << hcnt++ << "] ";
      out << std::setw(wid) << std::right << sout.str();
    }
    out << std::endl;

    // Data output
    //
    for (int i=0; i<numT; i++) {
      out << std::setw(wid) << std::setprecision(6) << times[i];
      for (auto u : mean) {
	auto k = u.first;
	double acc = 0.0, dat = data[k][i];

	for (int j=0; j<ncomp; j++) acc += RC[k](i, j);

	// Unnormalize
	//
	if (vm.count("totVar")) {
	  acc = acc*sqrt(totVar) + mean[k];
	  dat = dat*sqrt(totVar) + mean[k];
	} else {
	  acc = acc*sqrt(var[k]) + mean[k];
	  dat = dat*sqrt(var[k]) + mean[k];
	}

	out << std::setw(wid) << std::setprecision(6) << acc
	    << std::setw(wid) << std::setprecision(6) << dat;
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

    int n = 0;
    for (auto u : mean) {
      std::ostringstream fname1, fname2, suffix;
      auto k = u.first;
      suffix << std::get<0>(k) << "_" << std::get<1>(k) << "_" << std::get<2>(k) << "_" << lab[std::get<3>(k)];
      fname1 << prefix << ".recon."       << suffix.str();
      fname2 << prefix << ".recon_accum." << suffix.str();
      std::ofstream out1 (fname1.str());
      std::ofstream out2 (fname2.str());
      if (out1 and out2) {
	out1 << "# " << u.first << std::endl;
	out2 << "# " << u.first << std::endl;
	for (int i=0; i<numT; i++) {
	  out1 << std::setw(15) << std::setprecision(6) << times[i];
	  out2 << std::setw(15) << std::setprecision(6) << times[i];
	  double accum = 0.0;
	  for (int j=0; j<ncomp; j++) {
	    out1 << std::setw(15) << std::setprecision(6) << RC[k](i, j);
	    out2 << std::setw(15) << std::setprecision(6) << (accum += RC[k](i, j));
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

  if (vm.count("group")) {

				// Triple grouping
    std::vector<std::set<int>> groups;
				// Line buffer
    std::string line;
    int elem;
				// Open group file
    std::ifstream gf(grpfile);
				// Get first line
    std::getline(gf, line);

    while (gf.good()) {
				// Remove leading, trailing and
				// repeated white space from line
      line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");

      if (line[0] != '#') {	// Skip if commented

				// Remove trailing comments
	auto loc = line.find("#");
	if (loc != std::string::npos) line = line.substr(0, loc);

	std::set<int> group;
	std::istringstream sin(line);
	while (sin.good()) {	// Enter elements into the set
	  try {
	    sin >> elem;
	    group.insert(elem);
	  }
	  catch (std::exception& e) {
	    break;
	  }
	}
	groups.push_back(group); // Save the group
      }
      std::getline(gf, line);  // Get line containing next group . . .
    }

    if (groups.size()) {

      const unsigned int magic_word = 0x5ecede4;

      for (int N=0; N<Coefs.size(); N++) {

	std::ostringstream fname;
	fname << prefix << "." << N << ".recon_bin";
	std::ofstream out (fname.str());

	if (out) {

	  unsigned LMsize = LMset.size();
	  int      Gsize  = groups.size();

	  out.write((const char *)&magic_word, sizeof(unsigned int));
	  out.write((const char *)&LMsize,     sizeof(unsigned int));
	  out.write((const char *)&numT,       sizeof(int));
	  out.write((const char *)&nmax,       sizeof(int));
	  out.write((const char *)&Gsize,      sizeof(int));
	  out.write((const char *)&times[0],   sizeof(double)*numT);

	  for (auto lm : LMset) {
	    int LL = lm.first;
	    int MM = lm.second;
	    out.write((const char *)&LL, sizeof(int));
	    out.write((const char *)&MM, sizeof(int));
	  }

	  int igrp = 0;
	  for (auto v : groups) {
	    
	    for (unsigned i=0; i<numT; i++) {

	      for (auto lm : LMset) {
		int LL = lm.first;
		int MM = lm.second;

		for (int n=0; n<nmax; n++) {

		  Key c(LL, MM, n, 0, N), s(LL, MM, n, 1, N);
		  auto rc = RC.find(c);
		  auto rs = RC.find(s);
		  double valC = 0.0, valS = 0.0;
		  if (rc != RC.end()) {
		    for (auto u : v) {
		      if (u >-1 and u < ncomp) valC += RC[c](i, u);
		    }
		  }
		  if (rs != RC.end() and MM>0) {
		    for (auto u : v) {
		      if (u >-1 and u < ncomp) valS += RC[s](i, u);
		    }
		  }
		  
		  // Retrend
		  //
		  if (vm.count("totVar")) {
		    valC = valC*sqrt(totVar) + mean[c];
		    valS = valS*sqrt(totVar) + mean[s];
		  } else {
		    valC = valC*sqrt(var[c]) + mean[c];
		    valS = valS*sqrt(var[s]) + mean[s];
		  }
		  
		  out.write((const char *)&valC, sizeof(double));
		  out.write((const char *)&valS, sizeof(double));
		}
		// radial order loop
	      }
	      // L, M loop
	    }
	    // T loop
	    igrp++;
	  }
	  // Group loop
	}
	// Output file okay
	else {
	  std::cout << "Could not open <" << fname.str() << ">" << std::endl;
	}
      }

      // File loop
      //
      for (int N=0; N<Coefs.size(); N++) {

	std::ostringstream fname;
	fname << prefix << "." << N << ".recon_cmpl";
	std::ofstream out (fname.str());

	if (out) {

	  unsigned LMsize = LMset.size();
	  int      Gsize  = 1;

	  out.write((const char *)&magic_word, sizeof(unsigned int));
	  out.write((const char *)&LMsize,     sizeof(unsigned int));
	  out.write((const char *)&numT,       sizeof(int));
	  out.write((const char *)&nmax,       sizeof(int));
	  out.write((const char *)&Gsize,      sizeof(int));
	  out.write((const char *)&times[0],   sizeof(double)*numT);
	  
	  for (auto lm : LMset) {
	    int LL = lm.first;
	    int MM = lm.second;
	    out.write((const char *)&LL, sizeof(int));
	    out.write((const char *)&MM, sizeof(int));
	  }

	  std::set<int> comple;
	  for (int n=0; n<ncomp; n++) {
	    bool found = false;
	    for (auto s : groups) {
	      if (s.find(n) != s.end()) found = true;
	    }
	    if (!found) comple.insert(n);
	  }
	  
	  for (unsigned i=0; i<numT; i++) {

	    for (auto lm : LMset) {
	      int LL = lm.first;
	      int MM = lm.second;
	      
	      for (int n=0; n<nmax; n++) {
		Key c(LL, MM, n, 0, N), s(LL, MM, n, 1, N);
		auto rc = RC.find(c);
		auto rs = RC.find(s);
		double valC = 0.0, valS = 0.0;
		if (rc != RC.end()) {
		  for (auto u : comple) {
		    if (u >-1 and u < ncomp) valC += RC[c](i, u);
		  }
		}
		if (rs != RC.end() and MM>0) {
		  for (auto u : comple) {
		    if (u >-1 and u < ncomp) valS += RC[s](i, u);
		  }
		}
		
		// Retrend
		//
		if (vm.count("totVar")) {
		  valC = valC*sqrt(totVar) + mean[c];
		  valS = valS*sqrt(totVar) + mean[s];
		} else {
		  valC = valC*sqrt(var[c]) + mean[c];
		  valS = valS*sqrt(var[s]) + mean[s];
		}
		
		out.write((const char *)&valC, sizeof(double));
		out.write((const char *)&valS, sizeof(double));
	      }
	      // radial order loop
	    }
	    // L, M loop
	  }
	  // T loop
	}
	// Open file
	else {
	  std::cout << "Could not open <" << fname.str() << ">" << std::endl;
	}
      }
      // File loop
    }
    // Have groups
  }

  if (vm.count("coefs")) {

    for (auto lm : LMset) {
      int LL = lm.first;
      int MM = lm.second;

      // File loop
      //
      for (int N=0; N<Coefs.size(); N++) {

	std::ostringstream filename;
	filename << prefix << "_tot_" << LL << "_" << MM << "_" << N << ".coefs";
	out.open(filename.str());

	int q = 1;
	if (MM) q = 2;

	if (out) {
	  out << "# L=" << LL << " M=" << M << " coefs=" << dfiles[N]
	      << std::endl;

	  for (int i=0; i<numT; i+=stride) {
	    out << std::setw(15) << std::setprecision(6) << times[i] << std::endl;
	    out << std::setw(15) << NMAX << std::endl;

	    for (int n=0; n<NMAX; n++) {
	      Key c(LL, MM, n, 0, N);
	      auto rc = RC.find(c);
	      if (rc == RC.end())
		out << std::setw(15) << std::setprecision(6) << 0.0;
	      else {
		double acc = 0.0;
		for (int j=0; j<ncomp; j++) acc += RC[c](i, j);
		out << std::setw(15) << std::setprecision(6) << acc;
	      }
	    }
	    out << std::endl;

	    if (MM) {
	      for (int n=0; n<NMAX; n++) {
		Key s(LL, MM, n, 1, N);
		auto rs = RC.find(s);
		if (rs == RC.end())
		  out << std::setw(15) << std::setprecision(6) << 0.0;
		else {
		  double acc = 0.0;
		  for (int j=0; j<ncomp; j++) acc += RC[s](i, j);
		  out << std::setw(15) << std::setprecision(6) << acc;
		}
	      }
	      out << std::endl;
	    }
	  }
	  out.close();
	} else {
	  std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	  exit(-1);
	}
      }
      // File loop
    }
    // LM set loop

    int pr = vm.count("pairs") ? 2 : 1;

    for (int j=0; j<ncomp; j+=pr) {

      for (auto lm : LMset) {
	int LL = lm.first;
	int MM = lm.second;

	// File loop
	//
	for (int N=0; N<Coefs.size(); N++) {

	  std::ostringstream filename;
	  filename << prefix << "_" << LL << "_" << MM << "_" << j
		   << "_" << N << ".coefs";
	  out.open(filename.str());

	  if (out) {
	    for (int i=0; i<numT; i+=stride) {
	      out << std::setw(15) << std::setprecision(6) << times[i] << std::endl;
	      out << std::setw(15) << NMAX << std::endl;
	      
	      for (int n=0; n<NMAX; n++) {
		Key c(LL, MM, n, 0, N);
		auto rc = RC.find(c);
		if (rc == RC.end())
		  out << std::setw(15) << std::setprecision(6) << 0.0;
		else {
		  double val = RC[c](i, j);
		  if (pr>1) val += RC[c](i, j+1);
		  out << std::setw(15) << std::setprecision(6) << val;
		}
	      }
	      out << std::endl;

	      if (MM) {
		for (int n=0; n<NMAX; n++) {
		  Key s(LL, MM, n, 1, N);
		  auto rs = RC.find(s);
		  if (rs == RC.end())
		    out << std::setw(15) << std::setprecision(6) << 0.0;
		  else {
		    double val = RC[s](i, j);
		    if (pr>1) val += RC[s](i, j+1);
		    out << std::setw(15) << std::setprecision(6) << val;
		  }
		}
		out << std::endl;
	      }
	    }
	    out.close();
	  } else {
	    std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	    exit(-1);
	  } // file open failure
	}
	// File loop
	
      }
      // LM pairs

    }
    // ncomp loop

  }
  // coefs

  return 0;
}
