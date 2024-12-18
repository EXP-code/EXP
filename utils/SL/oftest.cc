#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <random>
#include <string>
#include <cmath>

#include <OrthoFunction.H>	// Orthogonal basis
#include <cxxopts.H>


// Generate a 2d exponential distribution
//
class genE
{
private:

  static constexpr double TOL = 1.0e-12;
  std::shared_ptr<std::mt19937> gen;
  std::shared_ptr<std::uniform_real_distribution<>> uniform;

  double R, a;
  double f(double x) { return R - (1.0 - (1.0 + x)*exp(-x)); };
  double df(double x){ return -x*exp(-x); };
  
public:

  genE(double a) : a(a) {
    unsigned seed = std::random_device{}();
    gen = std::make_shared<std::mt19937>(seed);
    uniform = std::make_shared<std::uniform_real_distribution<>>(0.0, 1.0);
  }

  std::tuple<double, int> operator()()
  {
    R = (*uniform)(*gen);
    double x = sqrt(R);

    int n=0;
    for (; n<100; n++) {
      double delF = - f(x)/df(x);
      x += delF;
      if (fabs(delF)<TOL) break;
    }
  
    return {x*a, n};
  }
};

int main(int argc, char** argv)
{
  bool logr = false;
  int numr, nmax, nout, N;
  double A, rmin, rmax, scale, delta;
  unsigned seed;		// Will be inialized by /dev/random if
				// not set on the command line
  std::string filename;

  // Parse command line
  //
  cxxopts::Options options(argv[0],
			   "Computes an orthogonal set of functions from a giving weighting\n");

  options.add_options()
    ("h,help",     "Print this help message")
    ("logr",       "Plot output grid with logarithmic spacing")
    ("A,length",   "characteristic disk scale length",
     cxxopts::value<double>(A)->default_value("1.0"))
    ("D,delta",    "truncation smoothing scale",
     cxxopts::value<double>(delta)->default_value("0.005"))
    ("N,mc",       "number of particles for Monte Carlo realization",
     cxxopts::value<int>(N)->default_value("10000"))
    ("nmax",       "maximum number of radial harmonics in the expansion",
     cxxopts::value<int>(nmax)->default_value("10"))
    ("r,rmin",     "minimum radius for the SL grid",
     cxxopts::value<double>(rmin)->default_value("0.0001"))
    ("R,rmax",     "maximum radius for the SL grid",
     cxxopts::value<double>(rmax)->default_value("20.0"))
    ("s,scale",    "ortho function map scale factor",
     cxxopts::value<double>(scale)->default_value("0.01"))
    ("nout",       "number of points in the output grid per side",
     cxxopts::value<int>(nout)->default_value("40"))
    ("seed",    "Random number seed. Default: use /dev/random",
     cxxopts::value<unsigned>(seed))
    ("o,filename", "Output filename",
     cxxopts::value<std::string>(filename)->default_value("testof"))
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

  // Log spacing?
  //
  if (vm.count("logr")) logr = true;

  // Create the weight function instance
  //
  auto densfunc = [&](double r)
  {
    return exp(-r/A) *
      0.5*(1.0 + std::erf((rmax - 5.0*delta - r)/delta)) / (A*A);
  };
  
  // Generate the orthogonal function instance
  //
  OrthoFunction ortho(nmax, densfunc, rmin, rmax, scale, 2);


  // Output file for grid
  //
  std::ofstream out(filename + ".dat");
  
  // Define some representative limits
  //
  double Rmax = 4.0*A;
  
  // Grid spacing
  //
  if (logr) {
    rmin = log10(rmin);
    rmax = log10(rmax);
  }
  double dR   = (rmax - rmin)/(nout - 1);

  // Print the grid
  //
  for (int i=0; i<nout; i++) {
    double R = rmin + dR*i;
    if (logr) R = exp(R);

    out << std::setw(18) << R;
    for (auto v : ortho(R))
      out << std::setw(18) << v;
    out << std::endl;
  }
	
  // Check orthogonality for a realized weight distribution
  //
  if (N) {
    genE gen(A);
    std::vector<double> coef(nmax, 0.0);
  
    for (int n=0; n<N; n++) {
      auto [R, m] = gen();	// Generate a point
      auto p = ortho(R);	// Get the orthogonal functions
      for (int i=0; i<nmax; i++) coef[i] += 1.0/N*p[i];
    }
    
    for (int i=0; i<nmax; i++) {
      std::cout << std::setw( 6) << i
		<< std::setw(18) << coef[i]
		<< std::endl;
    }
  }

  std::cout << std::endl
	    << ortho.testOrtho()
	    << std::endl;

  return 0;
}
