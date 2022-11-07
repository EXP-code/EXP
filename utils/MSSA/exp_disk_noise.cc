//
// M-SSA for EXP [cylindrical] coefficients
//
// Updated to use fixed-rank approximation for the SVD to save time
// and space.  Uses the approximate SVD randomized algorithm from
// Halko, Martinsson, and Tropp by default.  Use BDCSVD or Jacobi
// flags for the standard methods.
//
// Filter version
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

#include "Coefs.H"
#include "RedSVD.H"
#include <config_exp.h>
#include <cxxopts.H>

typedef std::tuple<unsigned, unsigned, unsigned> Key;
std::map<Key, std::vector<double>> coefs;
std::map<double, CylCoefsPtr> rawc;
std::vector<double> times;

template<typename T, typename U, typename V>
std::ostream& operator<< (std::ostream& out, const std::tuple<T, U, V>& t)
{
  const char lab[] = {'c', 's'};
  std::ostringstream sout;
  sout << '(' << std::get<0>(t)
       << ',' << std::get<1>(t)
       << ',' << lab[std::get<2>(t)]
       << ')';
  out << sout.str();
  return out;
}


std::map<std::string, int>
cylinder_read(const std::string& file, unsigned stride=1)
{
  std::map<std::string, int> ret;
  std::ifstream in(file);

  unsigned counter = 0;

  while (in.good()) {
    if (in.eof()) break;	// Check for end of file
    CylCoefsPtr c = std::make_shared<CylCoefs>();
    if (not c->read(in)) break;
    if (counter++ % stride == 0) rawc[c->time] = c;
  }

  int ntimes = 0;

  if (counter) {

    int mmax = rawc.begin()->second->mmax;
    int nmax = rawc.begin()->second->nmax;

    ntimes = rawc.size();

    for (auto v : rawc) {
      times.push_back(v.second->time);
    }

    for (int m=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	Key key0 = {m, n, 0}, key1 = {m, n, 1};
	coefs[key0].resize(ntimes);
	if (m) coefs[key1].resize(ntimes);
	for (int i=0; i<times.size(); i++) {
	  coefs[key0][i] = rawc[times[i]]->cos_c[m][n];
	  if (m) coefs[key1][i] = rawc[times[i]]->sin_c[m][n];
	}
      }
    }

    ret["mmax"]   = mmax;
    ret["nmax"]   = nmax;
    ret["ntimes"] = ntimes;
  }

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
  std::vector<int> Mvec;
  double tmin = -std::numeric_limits<double>::max();
  double tmax =  std::numeric_limits<double>::max();
  double evtol;


  //--------------------------------------------------
  // Command-line parsing
  //--------------------------------------------------

  std::string cmd_line;
  for (int i=0; i<argc; i++) {
    cmd_line += argv[i];
    cmd_line += " ";
  }

  cxxopts::Options options("exp_disk_noise","This routine uses M-SSA to analyze basis coefficients\nNB: *.recon and *.recon_diff are in EXP coefficient output format\nOptions");

  options.add_options()
    ("o,output",     "output file prefix",
     cxxopts::value<std::string>(prefix)->default_value("disk_noise"))
    ("d,datafile",   "input data file name",
     cxxopts::value<std::string>(datafile)->default_value("eof_coef_ascii.dat"))
    ("W,numW",       "window size",
     cxxopts::value<int>(numW)->default_value("10"))
    ("M,Mvec",       "vector list of harmonic orders",
     cxxopts::value<std::vector<int>>(Mvec))
    ("n,nmin",       "minimum radial order",
     cxxopts::value<int>(nmin)->default_value("0"))
    ("N,nmax",       "maximum radial order",
     cxxopts::value<int>(nmax)->default_value("99999"))
    ("P,npc",        "maximum number of eigenvectors",
     cxxopts::value<int>(npc)->default_value("99999"))
    ("e,evtol",      "cut on cumulative variance for eigenvalue sum",
     cxxopts::value<double>(evtol)->default_value("0.01"))
    ("t,Tmin",       "Minimum time in series",
     cxxopts::value<double>(tmin))
    ("T,Tmax",       "Maximum time in series",
     cxxopts::value<double>(tmax))
    ("Jacobi",       "Use the standard accurate but slow Jacobi method for SVD computation")
    ("BDCSVD",       "Use the bidiagonal divide and conquer method rather than Jacobi")
    ("Traj",         "Use RedSVD to decompose the trajectory matrix")
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
  std::map<std::string, int> param = cylinder_read(datafile);

  // Check that coefficient data series exist
  //
  if (rawc.size()==0) {
    std::cout << "No valid data in <" << datafile << ">" << std::endl;
    exit(-1);
  }

  int NMAX = param["nmax"];
  int MMAX = param["mmax"];
  int NTIM = param["ntimes"];

  if (nmax==0) nmax = NMAX;
  else         nmax = std::min<int>(nmax, NMAX);

  if (Mvec.size()==0) {
    for (int l=0; l<=MMAX; l++) Mvec.push_back(l);
  }

  numT = times.size();

  assert(NTIM == numT);		// Sanity check

  std::cout << "Number of times: " << numT << std::endl;


  // Coefficient copies
  //
  int NRAW = rawc.size();
  
  assert(NTIM == NRAW);

  std::vector<CylCoefs> cur(NRAW), dff(NRAW);

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

  // Cache M=0 value
  //
  double totS0 = 0.0;

  // M value loop
  //
  for (auto M : Mvec) {
    
    std::map<Key, std::vector<double> > data;
    std::map<Key, double> mean, var;

    int ndim = nmax - nmin;
    if (M>0) ndim *= 2;

    for (int n=0; n<ndim; n++) {
      Key key;
      int N = nmin + n;
      if (M) key = {M, N/2, N%2}; // Cosines (0) and Sines (1)
      else   key = {M, N,     0}; // Cosines (0)

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
    }

    if (false) {
      std::cout << "Using coefficient keys (M, n, cs):";
      for (auto v : mean) std::cout << " " << v.first;
      std::cout << std::endl;
    }

    if (false) {
      std::cout << std::string(60, '-') << std::endl
		<< "Mean values" << std::endl;
      for (auto v : mean) std::cout << std::setw(18) << v.first
				    << std::setw(18) << v.second
				    << std::endl;
      std::cout << std::string(60, '-') << std::endl;
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
	  if (totPow>0.0) v /= totPow;
	}
      }
      
    } else if (vm.count("totVar")) {
      useTotVar = true;
      
      for (auto & u : mean) {
	Key k = u.first;
	mean[k] /= numT;
	var[k]   = std::max<double>(var[k]/numT - mean[k]*mean[k], 0.0);
	totVar  += var[k];
	var[k]   = std::sqrt(var[k]);
      }
      
      totVar = std::sqrt(totVar);

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
	mean[k] /= numT;
	var[k]   = std::max<double>(var[k]/numT - mean[k]*mean[k], 0.0);
	var[k]   = std::sqrt(var[k]);
	for (auto & v : data[k]) {
	  v -= mean[k];
	  if (var[k]>0.0) v /= var[k];
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
    std::cout << "---- Eigenvalues M=" << M                 << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    
    // Compute total variance
    double tot = 0.0;
    for (int j=0; j<S.size(); j++) tot += S(j);
    std::cout << "Total=" << tot << std::endl;

    if (M==0) totS0 = tot;
    
    double cum = 0.0;
    int maxI = S.size();
    for (int j=0; j<S.size(); j++) {
      std::cout << std::setw( 5) << j
		<< std::setw(18) << S(j)
		<< std::setw(18) << (cum += S(j))/tot
		<< std::endl;

      // Set maxI?
      if (maxI==S.size()) {
	if (vm.count("totPow") and vm.count("select")) {
	  if (S(j) < evtol*totS0) maxI = j;
	} else {
	  if (1.0 - cum/tot < evtol) maxI = j;
	}
      }
    }

    // Report the EV truncation value
    //
    std::cout << std::endl
	      << "EV[" << M << "] cut at n=" << maxI
	      << std::endl << std::endl;

    // Cut the entire subspace; continue to next L value
    //
    if (maxI==0) continue;


    // Save the eigenvalues to a file
    //
    std::ostringstream filename;
    filename << prefix << "." << M << ".ev";
    std::ofstream out(filename.str());
    if (out) {
      cum = 0.0;
      for (int j=0; j<S.size(); j++) {
	out << std::setw( 5) << M
	    << std::setw( 5) << j
	    << std::setw(18) << S(j)
	    << std::setw(18) << (cum += S(j))/tot
	    << std::endl;
      }
      out.close();
    } else {
      std::cout << "Could not open <" << filename.str() << ">" << std::endl;
      exit(-1);
    }
    
    if (vm.count("ev")==0) {

      // Principal components
      //
      Eigen::MatrixXd PC = Y * U;

      int ncomp = std::min<int>({numW, npc, maxI, static_cast<int>(PC.cols())});

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
	    for (int j=L; j<U; j++)
	      RC[u.first](i, w) += Z(i, j) * rho[u.first](j, w) * W;
	  }
	}
      }

      std::map<double, CylCoefsPtr>::iterator it = rawc.begin();
      
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

	  int m = std::get<0>(u.first), n = std::get<1>(u.first);

	  if (std::get<2>(u.first) == 0) {
	    double val = acc*disp;	  // Restore variance
	    if (useMean) val += u.second; // Restore mean
	    cur[i].cos_c[m][n] = val;	  // Reconstructed
	    dff[i].cos_c[m][n] = val - cf->cos_c[m][n]; // Difference from orig
	  }
	  else {
	    double val = acc*disp;        // Restore variance
	    if (useMean) val += u.second; // Restore mean
	    cur[i].sin_c[m][n] = val;	  // Reconstructed
	    dff[i].sin_c[m][n] = val - cf->sin_c[m][n]; // Difference from orig
	  }
	}
	// END: key loop
	
	it++;
      }
      // END: time loop


      // This is for debugging
      //
      if (true) {
	std::ostringstream filename;
	filename << prefix << "." << M << ".recon_ascii";
	std::ofstream out(filename.str());
	if (out) {
	  const int wid = 16;
	  
	  out << "# " << std::right << std::setw(wid-1) << "Time ";
	  for (auto u : mean) {
	    std::ostringstream sout; sout << u.first << "(r)";
	    out << std::right << std::setw(wid-1) << sout.str() << " ";
	    sout.str(""); sout << u.first << "(d)";
	    out << std::right << std::setw(wid-1) << sout.str() << " ";
	  }
	  
	  out << std::endl;
      
	  out << "# " << std::right << std::setw(wid-1) << "[1] ";
	  int i=2;
	  for (auto u : mean) {
	    std::ostringstream sout; sout << "[" << i++ << "]";
	    out << std::right << std::setw(wid-1) << sout.str() << " ";
	    sout.str(""); sout << "[" << i++ << "]";
	    out << std::right << std::setw(wid-1) << sout.str() << " ";
	  }
	  out << std::endl;
	
	  double disp = 0.0;
	  if (useTotVar) disp = totVar;
	  if (useTotPow) disp = totPow;
	  
	  for (int i=0; i<numT; i++) {
	    out << std::right << std::setw(wid) << std::setprecision(6) << times[i];
	    for (auto u : mean) {
	      if (not useTotVar and not useTotPow) disp = var[u.first];
	      double off = 0.0;
	      if (useTotVar and useMean) off = u.second;
	      
	      double acc = 0.0;
	      for (int j=0; j<ncomp; j++) acc += RC[u.first](i, j);
	      out << std::right << std::setw(wid) << std::setprecision(6) << acc*disp + off;
	      out << std::right << std::setw(wid) << std::setprecision(6) << data[u.first][i]*disp + off;
	    }
	    out << std::endl;
	  }
	  out.close();
	} else {
	  std::cout << "Could not open <" << filename.str() << ">" << std::endl;
	  exit(-1);
	}
      }
      // END: debug
    }
    // END: not ev only
  }
  // END: M loop



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
