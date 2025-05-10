#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <memory>
#include <random>
#include <string>
#include <cmath>

#include <OrthoFunction.H>	// Orthogonal basis
#include <interp.H>		// Linear interpolation
#include <cxxopts.H>


// Recursion relation evaluation for generalized Laguerre functions
//
double generalized_laguerre(int n, double alpha, double x)
{
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return 1.0 + alpha - x;
  } else {
    return ((2.0 * n - 1.0 + alpha - x) * generalized_laguerre(n - 1, alpha, x) - (n - 1.0 + alpha) * generalized_laguerre(n - 2, alpha, x)) / n;
  }
}


// Class to generate 2d distribution
//
class gen2D
{
protected:

  static constexpr double TOL = 1.0e-12;
  std::shared_ptr<std::mt19937> gen;
  std::shared_ptr<std::uniform_real_distribution<>> uniform;

  double a, phi, alpha;
  int m;

public:

  gen2D(double a, int m, double phi, double alpha) :
    a(a), m(m), phi(phi), alpha(alpha)
  {
    unsigned seed = std::random_device{}();
    gen = std::make_shared<std::mt19937>(seed);
    uniform = std::make_shared<std::uniform_real_distribution<>>(0.0, 1.0);
  }


  virtual std::tuple<double, double, int, int> operator()() = 0;
};

// Generate a 2d exponential distribution
//
class genE : public gen2D
{
private:

  double R;
  double f(double x) { return R - (1.0 - (1.0 + x)*exp(-x)); };
  double df(double x){ return -x*exp(-x); };
  
public:

  genE(double a, int m, double phi, double alpha) : gen2D(a, m, phi, alpha) {}

  std::tuple<double, double, int, int> operator()()
  {
    R = (*uniform)(*gen);
    double x = sqrt(R);

    // Realize a radius
    int n=0;
    for (; n<100; n++) {
      double delF = - f(x)/df(x);
      x += delF;
      if (fabs(delF)<TOL) break;
    }
  
    // Realize an angle
    double P = 2.0*M_PI*(*uniform)(*gen);
    int k=0;
    for (; k<1000; k++) {
      double fp = 0.5*(1.0 + cos((P-phi-alpha*x)*m));
      if (fp > (*uniform)(*gen)) break;
      P = 2.0*M_PI*(*uniform)(*gen);
    }
      
    return {x*a, P, n, k};
  }
};

// Generate a distribution from arrays
//
class genO : public gen2D
{
private:

  std::unique_ptr<Linear1d> interp;
  
public:

  genO(const std::string& filename,
       int m, double phi, double alpha) : gen2D(1.0, m, phi, alpha)
  {
    std::ifstream in (filename);
    std::vector<double> r, d;
    if (in) {
      std::string line;
      while (std::getline(in, line)) {
	std::istringstream iss(line);
	double x, y;
	iss >> x >> y;
	if (iss) {
	  r.push_back(x);
	  d.push_back(y);
	} else break;
      }
    }
    else throw std::runtime_error("Cannot open file " + filename);

    std::vector<double> pp(r.size(), 0.0);
    for (int n=1; n<r.size(); n++) {
      pp[n] = pp[n-1] + 0.5*(d[n-1]*r[n-1] + d[n]*r[n])*(r[n] - r[n-1]);
    }
    for (auto & v : pp) v /= pp.back();

    interp = std::make_unique<Linear1d>(pp, r);
    // Debug
    std::ofstream out("genO.dat");
    if (out) {
      for (int n=0; n<pp.size(); n++) {
	out << std::setw(18) << r[n]
	    << std::setw(18) << d[n]
	    << std::setw(18) << pp[n] << std::endl;
      }
    }
  }

  genO(std::vector<double> r, std::vector<double> d,
       int m, double phi, double alpha) : gen2D(1.0, m, phi, alpha)
  {
    std::vector<double> pp(r.size(), 0.0);
    for (int n=1; n<r.size(); n++) {
      pp[n] = pp[n-1] + 0.5*(d[n-1]*r[n-1] + d[n]*r[n])*(r[n] - r[n-1]);
    }
    for (auto & v : pp) v /= pp.back();

    interp = std::make_unique<Linear1d>(pp, r);
    // Debug
    std::ofstream out("genO.dat");
    if (out) {
      for (int n=0; n<pp.size(); n++) {
	out << std::setw(18) << r[n]
	    << std::setw(18) << d[n]
	    << std::setw(18) << pp[n] << std::endl;
      }
    }
  }

  std::tuple<double, double, int, int> operator()()
  {
    double r = interp->eval((*uniform)(*gen));

    // Realize an angle
    double P = 2.0*M_PI*(*uniform)(*gen);

    int k=0;
    if (m>0) {
      for (; k<1000; k++) {
	double fp = 0.5*(1.0 + cos((P-phi-alpha*r)*m));
	if (fp > (*uniform)(*gen)) break;
	P = 2.0*M_PI*(*uniform)(*gen);
      }
    }
      
    return {r, P, 0, k};
  }
};

// Generate a 2d uniform distribution; a uniform disk with an edge at a.
//
class genU : public gen2D
{
private:

public:

  genU(double a, int m, double phi, double alpha) : gen2D(a, m, phi, alpha) {}

  std::tuple<double, double, int, int> operator()()
  {
    double R = (*uniform)(*gen);
    double x = sqrt(R);

    // Realize an angle
    double P = 2.0*M_PI*(*uniform)(*gen);
    int k=0;
    for (; k<1000; k++) {
      double fp = 0.5*(1.0 + cos((P-phi-alpha*x)*m));
      if (fp > (*uniform)(*gen)) break;
      P = 2.0*M_PI*(*uniform)(*gen);
    }
      
    return {x*a, P, 0, k};
  }
};

int main(int argc, char** argv)
{
  bool logr = false;
  int numr, nmax, nout, N, M, mmax, index, knots;
  double A, rmin, rmax, scale, delta, pitch, phi, Rout;
  unsigned seed;		// Will be inialized by /dev/random if
				// not set on the command line
  std::string filename, modelname, coeffile, densprof, bodyfile;

  // Parse command line
  //
  cxxopts::Options options(argv[0],
			   "Computes an orthogonal set of functions from a giving weighting\n");

  options.add_options()
    ("h,help",     "Print this help message")
    ("logr",       "Plot output grid with logarithmic spacing")
    ("uniform",    "Generate a uniform distribution rather than an exponential")
    ("i,index",    "Use an orthogonal basis of index i for realization",
     cxxopts::value<int>(index)->default_value("0"))
    ("A,length",   "characteristic disk scale length",
     cxxopts::value<double>(A)->default_value("1.0"))
    ("D,delta",    "truncation smoothing scale",
     cxxopts::value<double>(delta)->default_value("0.005"))
    ("P,phi",      "Position angle for M>0 distribution (degrees)",
     cxxopts::value<double>(phi)->default_value("45.0"))
    ("N,mc",       "number of particles for Monte Carlo realization",
     cxxopts::value<int>(N)->default_value("10000"))
    ("M",          "azimuthal order for Monte Carlo realization",
     cxxopts::value<int>(M)->default_value("0"))
    ("k,knots",    "Inner product knots for Stieltjes procedure",
     cxxopts::value<int>(knots)->default_value("400"))
    ("mmax",       "maximum azimuthal order for reconstrution",
     cxxopts::value<int>(mmax)->default_value("4"))
    ("nmax",       "maximum number of radial harmonics in the expansion",
     cxxopts::value<int>(nmax)->default_value("18"))
    ("r,rmin",     "minimum radius for the SL grid",
     cxxopts::value<double>(rmin)->default_value("0.0001"))
    ("R,rmax",     "maximum radius for the SL grid",
     cxxopts::value<double>(rmax)->default_value("20.0"))
    ("s,scale",    "ortho function map scale factor",
     cxxopts::value<double>(scale)->default_value("1.0"))
    ("p,pitch",    "dPhi/dR in scale length units",
     cxxopts::value<double>(pitch)->default_value("0.0"))
    ("nout",       "number of points in the output grid per side",
     cxxopts::value<int>(nout)->default_value("40"))
    ("rout",       "half-size of output grid in scale lenghts",
     cxxopts::value<double>(Rout)->default_value("4.0"))
    ("seed",        "Random number seed. Default: use /dev/random",
     cxxopts::value<unsigned>(seed))
    ("model",      "Input model file",
     cxxopts::value<std::string>(modelname))
    ("o,filename", "Output filename",
     cxxopts::value<std::string>(filename)->default_value("oftest"))
    ("c,coeffile", "Input coefficient file",
     cxxopts::value<std::string>(coeffile))
    ("d,densprof", "Density profile for sampling",
     cxxopts::value<std::string>(densprof))
    ("b,bodyfile", "Body file for sampling",
     cxxopts::value<std::string>(bodyfile))
    ;
  
  //===================
  // Parse options
  //===================

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    return 2;
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    std::cout << options.help() << std::endl << std::endl;
    return 1;
  }

  // Set from /dev/random if not specified
  if (N>0 and vm.count("seed")==0) {
    seed = std::random_device{}();
  }

  phi *= M_PI/180.0;

  // Log spacing?
  //
  if (vm.count("logr")) logr = true;

  // File or expon
  //
  std::function<double(double)> densfunc;
  std::shared_ptr<Linear1d> interp, interp2;

  // Read a model file
  //
  if (vm.count("model")) {
    std::vector<double> r, d;
    std::ifstream in(modelname);
    if (not in) throw std::runtime_error("Error opening file: " + modelname);
    
    std::string line;
    while (std::getline(in, line)) {
      auto pos = line.find_first_of("!#");
      if (pos == std::string::npos) {
	std::istringstream iss(line);
	double x, y;
	iss >> x >> y;
	if (iss) {
	  r.push_back(x);
	  d.push_back(y);
	}
      }
    }
      
    // Compute interpolation functionoid
    //
    double rmin = r.front(), rmax = r.back();
    
    interp = std::make_shared<Linear1d>(r, d);
    densfunc = [&interp](double r)
    {
      return interp->eval(r);
    };
  }
  // Create the weight function for expon
  //
  else {
    if (vm.count("uniform"))
      densfunc = [&](double r)
	{
	  return 1.0;
	};
    else
      densfunc = [&](double r)
	{
	  return exp(-r/A) *
	    0.5*(1.0 + std::erf((rmax - 5.0*delta - r)/delta)) / (A*A);
	};
  }

  // Generate the orthogonal function instance
  //
  OrthoFunction ortho(nmax, densfunc, rmin, rmax, scale, 2);
  ortho.setKnots(knots);

  // Output file for grid
  //
  std::ofstream out(filename + ".dat");
  
  // Define some representative limits
  //
  double Rmax = Rout*A;
  
  // Grid spacing
  //
  if (logr) {
    rmin = log10(rmin);
    rmax = log10(rmax);
  }
  double dR   = (rmax - rmin)/(nout - 1);

  // Print the grid
  //
  double tmass = 0.0;
  for (int i=0; i<nout; i++) {
    double R = rmin + dR*i;
    if (logr) {
      double lR = exp(R - dR);
      R = exp(R);
      tmass += M_PI*(densfunc(R)*R*R + densfunc(lR)*lR*lR)*dR;
    } else {
      double lR = R - dR;
      tmass += M_PI*(densfunc(R)*R + densfunc(lR)*lR)*dR;
    }

    out << std::setw(18) << R;
    for (auto v : ortho(R))
      out << std::setw(18) << v;
    out << std::endl;
  }
  std::cout << "Total mass = " << tmass << std::endl;
	
  // Check orthogonality for a realized weight distribution
  //
  if (N or vm.count("bodyfile") or vm.count("coeffile")) {

    // Normalization for trig fcts
    constexpr double fac = 0.5*M_2_SQRTPI/M_SQRT2;

    // Coefficients
    std::vector<std::vector<std::complex<double>>> coef(mmax+1);
    for (int m=0; m<=mmax; m++) coef[m].resize(nmax, 0.0);

    // Read a coefficient file
    //
    if (vm.count("coeffile")) {
      std::ifstream in(coeffile);
      if (in.good()) {
	std::string line;
        while (std::getline(in, line)) {
	  std::istringstream sout(line);
	  int n; sout >> n;
	  if (sout and n < nmax) {
	    std::vector<std::complex<double>> c;
	    while (sout) {
	      double r, v;
	      sout >> r >> v;
	      c.push_back({r, v});
	    }
	    for (int mm=0; mm<=std::min<int>(mmax, c.size()); mm++) {
	      coef[mm][n] = c[mm];
	    }
	  }
	}
      } else {
	throw std::runtime_error("Error opening coefficient file: " + coeffile);
      }

      N = 0;
    }

    if (N or vm.count("bodyfile")) {
      // Distribution generator
      std::unique_ptr<gen2D> gen;
      if (vm.count("uniform"))
	gen = std::make_unique<genU>(A, M, phi, pitch);
      else if (vm.count("index")) {
	double ratio = 1.0;
	const int ngrid = 1000;
	double lrmin = log(rmin), lrmax = log(rmax);
	double dr = (lrmax - lrmin)/(ngrid - 1);
	for (int n=0; n<ngrid; n++) {
	  double r = exp(lrmin + dr*n);
	  auto p = ortho(r);
	  if (ratio*p(0) + p(index) < 0.0) ratio = -p(index)/p(0);
	}
	
	std::vector<double> rr(ngrid), dd(ngrid);
	for (int n=0; n<ngrid; n++) {
	  rr[n] = exp(lrmin + dr*n);
	  auto p = ortho(rr[n]);
	  dd[n] = ratio*p[0] + p[index];
	}
	std::cout << "Ratio: " << ratio << std::endl;
	gen = std::make_unique<genO>(rr, dd, M, phi, pitch);
      }
      else if (vm.count("densprof"))
	gen = std::make_unique<genO>(densprof, M, phi, pitch);
      else
	gen = std::make_unique<genE>(A, M, phi, pitch);

      if (vm.count("bodyfile")) {
	std::ifstream in(bodyfile);
	if (in) {
	  std::string line;
	  std::getline(in, line);
	  while (std::getline(in, line)) {
	    std::istringstream iss(line);
	    double ms, x, y;
	    iss >> ms >> x >> y;
	    auto p = ortho(sqrt(x*x + y*y));
	    double P = atan2(y, x);
	    for (int m=0; m<=mmax; m++) {
	      auto azi = std::exp(std::complex<double>(0, -P*m));
	      for (int i=0; i<nmax; i++)
		coef[m][i] += fac*ms*p(i)*azi;
	    }
	  }
	}
	else throw std::runtime_error("Error opening body file: " + bodyfile);
      } else {
	// Sample distribution
	for (int n=0; n<N; n++) {
	  // Generate a point
	  auto [R, P, nr, np] = (*gen)();
	  
	  // Get the orthogonal functions
	  auto p = ortho(R);
	  for (int m=0; m<=mmax; m++) {
	    auto azi = std::exp(std::complex<double>(0, -P*m));
	    for (int i=0; i<nmax; i++)
	      coef[m][i] += fac*tmass/N*p(i)*azi;
	  }
	}
      }
    }
    
    std::ofstream cof(filename + ".coef");
    for (int i=0; i<nmax; i++) {
      cof << std::setw(8) << i;
      for (int m=0; m<=mmax; m++)
	cof << std::setw(18) << std::abs(coef[m][i])
	    << std::setw(18) << std::arg(coef[m][i]);
      cof << std::endl;
    }

    int nxy = 100;
    std::ofstream out1(filename + ".mat");
    std::ofstream out2(filename + ".line");
    out1 << std::setw(6) << nxy << std::setw(6) << nxy << std::endl;
    for (int j=0; j<nxy; j++) {
      double y = -Rmax + 2.0*Rmax*j/(nxy-1);
      for (int i=0; i<nxy; i++) {
	double x = -Rmax + 2.0*Rmax*i/(nxy-1);
	    
	auto p = ortho(sqrt(x*x + y*y));	// Get the orthogonal functions
	double phi = atan2(y, x);

	out1 << std::setw(18) << x
	     << std::setw(18) << y;

	std::complex<double> tot = 0.0;
	for (int m=0; m<=mmax; m++) {
	  std::complex<double> c = 0.0;
	  for (int i=0; i<nmax; i++) c += coef[m][i]*p[i]*fac;
	  c *= std::exp(std::complex<double>(0, phi*m));
	  out1 << std::setw(18) << std::real(c);
	  tot += c;
	}
	out1 << std::setw(18) << std::real(tot) << std::endl;
      }

      auto p = ortho(std::abs(y));
      double phi = 0.5*M_PI;
      if (y<0.0) phi = -0.5*M_PI;

      std::complex<double> tot = 0.0;
      for (int m=0; m<=mmax; m++) {
	std::complex<double> c = 0.0;
	for (int i=0; i<nmax; i++) c += coef[m][i]*p[i]*fac;
	c *= std::exp(std::complex<double>(0, phi*m));
	out << std::setw(18) << std::real(c);
	tot += c;
      }
      out2 << std::setw(18) << y << std::setw(18) << std::real(tot) << std::endl;
    }
  }

  std::cout << std::endl
	    << "Orthogonality of the function at the grid points"
	    << std::endl;

  ortho.dumpOrtho("oftest.dump");

  if (vm.count("model")==0 and vm.count("uniform")==0) {
    std::ofstream out("oftest.laguerre");
    const int ngrid = 1000;
    double alpha = 1.0;
    double lrmin = log(rmin), lrmax = log(rmax);
    double dr = (lrmax - lrmin)/(ngrid - 1);
    for (int j=0; j<ngrid; j++) {
    double r = exp(lrmin + dr*j);
      out << std::setw(18) << r;
      for (int n=0; n<nmax; n++)
	out << std::setw(18) << generalized_laguerre(n, alpha, r)
	  * exp(0.5*(std::lgamma(n+1) - std::lgamma(n+alpha+1)));
      out << std::endl;
    }
  }

  std::cout << std::endl
	    << ortho.testOrtho()
	    << std::endl;

  return 0;
}
