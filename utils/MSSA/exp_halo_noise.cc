//
// M-SSA for EXP spherical coefficients
// Noise filter version
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>


#include <Eigen/Dense>

#include "Coefs.H"
#include "RedSVD.H"
#include "config.h"
#include "cxxopts.H"


typedef std::tuple<unsigned, unsigned, unsigned, unsigned> Key;
std::map<Key, std::vector<double>> coefs;
std::map<double, SphCoefsPtr> rawc;
std::vector<double> times;

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::tuple<T, T>& t)
{
  std::ostringstream sout;
  sout << '{' << std::get<0>(t)
       << ',' << std::get<1>(t)
       << '}';
  out << sout.str();
  return out;
}

template<typename T, typename U, typename V, typename W>
std::ostream& operator<< (std::ostream& out, const std::tuple<T, U, V, W>& t)
{
  const char lab[] = {'c', 's'};
  std::ostringstream sout;
  sout << '(' << std::get<0>(t)
       << ',' << std::get<1>(t)
       << ',' << std::get<2>(t)
       << ',' << lab[std::get<3>(t)]
       << ')';
  out << sout.str();
  return out;
}


typedef std::tuple<unsigned, unsigned> LMkey;

//! Return structure for EXP coefficient table read
struct SphHeader
{
  std::string id;
  double      scale;
  int         lmax;
  int         nmax;
  int         ntimes;
};

//! Read EXP spherical coefficient table
SphHeader
spherical_read(const std::string& file,
	       const double tmin,
	       const double tmax,
	       unsigned stride=1)
{
  SphHeader ret;

  std::ifstream in(file);

  unsigned counter = 0;

  while (in.good()) {
    if (in.eof()) break;	// Check for end of file
    SphCoefsPtr c = std::make_shared<SphCoefs>();
    if (not c->read(in)) break;
    if (counter==0) {
      ret.id    = c->header.id;
      ret.scale = c->header.scale;
      ret.lmax  = c->header.Lmax;
      ret.nmax  = c->header.nmax;
    }
    if (c->header.tnow > tmax) break;
    if (c->header.tnow > tmin) {
      if (counter++ % stride == 0) rawc[c->header.tnow] = c;
    }
  }

  int ntimes = 0;

  if (counter) {

    int lmax = rawc.begin()->second->header.Lmax;
    int nmax = rawc.begin()->second->header.nmax;

    ntimes = rawc.size();

    for (auto v : rawc) {
      times.push_back(v.second->header.tnow);
    }

    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	for (int n=0; n<nmax; n++) {
	  LMkey lmk  = {l, m};
	  Key   key0 = {l, m, n, 0}, key1 = {l, m, n, 1};
	  coefs[key0].resize(ntimes);
	  if (m) coefs[key1].resize(ntimes);
	  for (int i=0; i<times.size(); i++) {
	    coefs[key0][i] = rawc[times[i]]->cos_c[lmk][n];
	    if (m) coefs[key1][i] = rawc[times[i]]->sin_c[lmk][n];
	  }
	}
      }
    }
  }

  ret.ntimes = ntimes;

  return ret;
}


int
main(int argc, char **argv)
{
  //--------------------------------------------------
  // Parameters
  //--------------------------------------------------
  std::string prefix, datafile, keyfile, grpfile;
  int numT, numW, ndim, nmin, nmax, npc;
  std::vector<int> Lvec;
  double evtol, tmin = 0.0, tmax = std::numeric_limits<double>::max();


  //--------------------------------------------------
  // Command-line parsing
  //--------------------------------------------------

  std::string cmd_line;
  for (int i=0; i<argc; i++) {
    cmd_line += argv[i];
    cmd_line += " ";
  }

  cxxopts::Options options("exp_halo_noise","This routine uses M-SSA to analyze basis coefficients\nNB: *.recon and *.recon_diff are in EXP coefficient output format\nOptions");

  options.add_options()
    ("o,output",     "output file prefix",
     cxxopts::value<std::string>(prefix)->default_value("halo_noise"))
    ("d,datafile",   "EXP coefficient data file",
     cxxopts::value<std::string>(datafile)->default_value("outcoef.run0.dark"))
    ("W,numW",       "window size",
     cxxopts::value<int>(numW)->default_value("10"))
    ("L,Lvec",       "vector list of harmonic orders",
     cxxopts::value<std::vector<int>>(Lvec))
    ("n,nmin",       "minimum radial order",
     cxxopts::value<int>(nmin)->default_value("0"))
    ("N,nmax",       "maximum radial order",
     cxxopts::value<int>(nmax)->default_value("99999"))
    ("P,npc",        "maximum number of eigenvectors",
     cxxopts::value<int>(npc)->default_value("99999"))
    ("e,evtol",      "cut on cumulative variance for eigenvalue sum",
     cxxopts::value<double>(evtol)->default_value("0.01"))
    ("t,Tmin",       "Minimum time in series",
     cxxopts::value<double>(tmax))
    ("T,Tmax",       "Maximum time in series",
     cxxopts::value<double>(tmax))
    ("Jacobi",     "Use the standard accurate but slow Jacobi method for SVD computation")
    ("BDCSVD",     "Use the bidiagonal divide and conquer method rather than Jacobi")
    ("Traj",       "Use RedSVD to decompose the trajectory matrix")
    ("v,version",    "show version")
    ("X,noCommand",  "do not save command line")
    ("E,ev",         "exit after computing eigenvalues")
    ("z,zero",       "zero unfiltered coefficients")
    ("channels",     "print principle components by channels")
    ("totVar",       "use total variance for normalization")
    ("totPow",       "use total power for normalization")
    ("noMean",       "no mean detrending for totPow")
    ("select",       "truncate EV series using relative power threshold")
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

  if (vm.count("help")) {
    std::cout << std::string(60, '-') << std::endl;
    std::cout << options.help() << std::endl;
    std::cout << std::string(60, '-') << std::endl << std::endl;
    return 0;
  }
  
  if (vm.count("noCommand")==0) {
    std::string cmdFile = prefix + ".cmd_line";
    std::ofstream cmd(cmdFile.c_str());
    if (!cmd) {
      std::cerr << argv[0] << ": error opening <" << cmdFile
		<< "> for writing" << std::endl;
    } else {
      cmd << cmd_line << std::endl;
    }
    
    cmd.close();
  }

  std::cout << "Eigen is using " << Eigen::nbThreads()
	    << " threads" << std::endl;

  // Attempt to read data
  //
  SphHeader param = spherical_read(datafile, tmin, tmax);

  // Check that coefficient data series exist
  //
  if (rawc.size()==0) {
    std::cout << "No valid data in <" << datafile << ">" << std::endl;
    exit(-1);
  }

  int LMAX = param.lmax;
  int NMAX = param.nmax;
  int NTIM = param.ntimes;

  if (nmax==0) nmax = NMAX;
  else         nmax = std::min<int>(nmax, NMAX);

  if (Lvec.size()==0) {
    for (int l=0; l<=LMAX; l++) Lvec.push_back(l);
  }

  numT = times.size();

  assert(NTIM == numT);		// Sanity check

  std::cout << "Number of times: " << numT << std::endl;


  // Coefficient copies
  //
  int NRAW = rawc.size();
  
  assert(NTIM == NRAW);

  std::vector<SphCoefs> cur(NRAW), dff(NRAW);

  {
    int n = 0;
    for (auto it : rawc) {

      // Clone coefficient stanzas and zero coefficients
      //
      cur[n].clone(it.second);
      if (vm.count("zero")) cur[n].zero();

      dff[n].clone(it.second);
      if (vm.count("zero")) dff[n].zero();
      
      n++;
    }
  }

  // Compute total power from coefficient database
  //
  double totPow    = 0.0;
  bool   useTotPow = false;
  bool   useMean   = true;
    
  // Average total power per time slice for time-series normalization
  //
  if (vm.count("totPow")) {
    useTotPow = true;
    if (vm.count("noMean")) useMean = false;
    for (auto & u : coefs) {
      for (auto & v : u.second) totPow += v*v;
    }

    totPow = std::sqrt(totPow/numT + std::numeric_limits<double>::min());
  }

  // L value loop
  //
  for (auto L : Lvec) {
    
    std::map<Key, std::vector<double> > data;
    std::map<Key, double> mean, var;

    int Mmin = 0, Mmax = L;

    for (int mm=Mmin; mm<=Mmax; mm++) {

      for (int n=nmin; n<nmax; n++) {

	int cmax = 2;		// Default: cosines and sines
	if (mm==0) cmax = 1;	// Cosines only of m=0
	for (int c=0; c<cmax; c++) {

	  Key key({L, mm, n, c});

	  if (coefs.find(key) == coefs.end()) {
	    std::cout << "Key " << key << " not found in coefs, keys are" << std::endl;
	    for (auto v : coefs) std::cout << v.first << std::endl;
	    exit(-1);
	  }

	  mean[key] = 0.0;
	  var [key] = 0.0;

	  for (int i=0; i<numT; i++) {
	    double y = coefs[key][i];
	    mean[key] += y;
	    var [key] += y*y;
	    data[key].push_back(y);
	  }
	  // END: time loop

	}
	// END: Cosine/sine loop

	ndim++;
      }
      // END: N loop
    }
    // END: M loop

    if (false) {
      std::cout << "Using coefficient keys (l, m, n, cs):";
      for (auto v : mean) std::cout << " " << v.first;
      std::cout << std::endl;
    }

    // Detrend all channels using specified method
    //
    double totVar = 0.0;
    bool useTotVar = false;
    
    if (vm.count("totPow")) {

      for (auto & u : mean) {
	Key k = u.first;
	mean[k] /= numT;
      }
      
      for (auto & u : mean) {
	Key k = u.first;
	for (auto & v : data[k]) {
	  if (useMean) v -= mean[k];
	  v /= totPow;
	}
      }

    } else if (vm.count("totVar")) {
      useTotVar = true;
      
      for (auto & u : mean) {
	Key k = u.first;
	mean[k] /= numT;
	var[k]   = std::max<double>(var[k]/numT - mean[k]*mean[k], 0.0);
	totVar  += var[k];
	var[k]   = std::sqrt(fabs(var[k]) + std::numeric_limits<double>::min());
      }
      
      totVar = std::sqrt(totVar + std::numeric_limits<double>::min());

      for (auto & u : mean) {
	Key k = u.first;
	for (auto & v : data[k]) {
	  v -= mean[k];
	  v /= totVar;
	}
      }
    } else {
      for (auto & u : mean) {
	Key k = u.first;
	mean[k] /= numT;
	var[k]   = std::max<double>(var[k]/numT - mean[k]*mean[k], 0.0);
	var[k]   = std::sqrt(fabs(var[k]) + std::numeric_limits<double>::min());
	for (auto & v : data[k]) {
	  v -= mean[k];
	  v /= var[k];
	}
      }
    }

    int nkeys = mean.size();
    int numK  = numT - numW + 1;

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

    // Store the solution
    // ------------------
    Eigen::VectorXd S;
    Eigen::MatrixXd U;

    // Compute the solution
    // --------------------

    // Use one of the built-in Eigen3 algorithms
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

      int rank = std::min<int>({static_cast<int>(cov.cols()), npc});

      RedSVD::RedSVD<Eigen::MatrixXd> svd(cov, rank);
      
      S = svd.singularValues();	// Get solution
      U = svd.matrixU();
    }
    
    std::cout << "shape U = " << U.rows() << " x "
	      << U.cols() << std::endl;
    
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "---- Eigenvalues L=" << L                 << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    
    // Rescale the SVD factorization by the Frobenius norm
    //
    S *= Scale;
    if (vm.count("Traj")) {
      for (int i=0; i<S.size(); i++) S(i) = S(i)*S(i)/numK;
    }
    
    double tot = 0.0;		// Compute total variance
    for (int j=0; j<S.size(); j++) tot += S(j);
    std::cout << "Frob norm=" << Scale << std::endl;
    std::cout << "Total var=" << tot   << std::endl;

    if (std::isnan(tot)) {
      std::cout << "Unexpected ERROR: bad total variance" << std::endl;
      exit(-1);
    }

    double cum = 0.0;
    int maxI = S.size();
    for (int j=0; j<S.size(); j++) {
      if (j<npc)
	std::cout << std::setw( 5) << j
		  << std::setw(18) << S(j)
		  << std::setw(18) << (cum += S(j))/tot
		  << std::endl;

      // Set maxI?
      if (maxI==S.size()) {
	if (vm.count("totPow") and vm.count("select")) {
	  if (S(j) < evtol*totPow) maxI = j;
	} else {
	  if (1.0 - cum/tot < evtol) maxI = j;
	}
      }
    }

    // Report the EV truncation value
    //
    std::cout << std::endl
	      << "EV[" << L << "] cut at n=" << maxI
	      << std::endl << std::endl;

    // Cut the entire subspace; continue to next L value
    //
    if (maxI==0) continue;

    std::ostringstream filename;

    // Save the eigenvalues to a file
    //
    filename << prefix << "_" << L << ".ev";
    std::ofstream out(filename.str());
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
      std::cout << "Could not open <" << filename.str() << ">" << std::endl;
      exit(-1);
    }

    if (vm.count("ev")) return 0;

    // Compute the PCs by projecting the data
    //
    Eigen::MatrixXd PC = Y * U;

    // Reconstructed time series
    //
    int ncomp = std::min<int>({numW, maxI, npc, static_cast<int>(PC.cols())});

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
      rho[u.first].resize(numW, ncomp);
      rho[u.first].fill(0.0);
      for (int i=0; i<numW; i++) {
	for (int j=0; j<ncomp; j++)
	  rho[u.first](i, j) = U(numW*n + i, j);
      }
      n++;
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
    
    std::map<double, SphCoefsPtr>::iterator it = rawc.begin();
      
    double disp;
    if (useTotVar) disp = totVar;
    if (useTotPow) disp = totPow;
      
    for (int i=0; i<numT; i++) {
	
      auto cf = it->second;
	
      // Write reconstruction to coefficient stanzas
      //
      for (auto u : mean) {
	if (not useTotVar and not useTotPow) disp = var[u.first];
	  
	double acc = 0.0;
	for (int j=0; j<ncomp; j++) acc += RC[u.first](i, j);
	
	LMkey lm = {std::get<0>(u.first), std::get<1>(u.first)};
	int    n = std::get<2>(u.first);
	
	if (std::get<3>(u.first) == 0) {
	  double val = acc*disp;	// Restore variance
	  if (useMean) val += u.second; // Restore mean
	  cur[i].cos_c[lm][n] = val;	// Reconstructed
	  dff[i].cos_c[lm][n] = val - cf->cos_c[lm][n]; // Difference from orig
	}
	else {
	  double val = acc*disp;	    // Restore variance
	  if (useMean) val += u.second;	    // Restore mean
	  cur[i].sin_c[lm][n] = val;	    // Reconstructed
	  dff[i].sin_c[lm][n] = val - cf->sin_c[lm][n]; // Difference from orig
	}
      }
      // END: key loop
	
      it++;			// Next stanza
    }
    // END: time loop

  }
  // END: L loop

  std::ostringstream filename1, filename2;
  filename1 << prefix << ".recon";
  filename2 << prefix << ".recon_diff";

  std::ofstream out1(filename1.str());
  std::ofstream out2(filename2.str());
    
  if (out1) {
    for (auto & c : cur) c.write(out1);
  } else {
    std::cout << "Could not open <" << filename1.str() << ">" << std::endl;
  }

  if (out2) {
    for (auto & c : dff) c.write(out2);
  } else {
    std::cout << "Could not open <" << filename2.str() << ">" << std::endl;
  }


  return 0;
}
