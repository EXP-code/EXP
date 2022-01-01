// Compile string: mpiCC -std=c++11 -O3 -o testQ testQuantile.cc Quantile.cc -lmpi

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <random>

#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/normal.hpp>

#include <cxxopts.H>		// Option parsing
#include <libvars.H>		// EXP library globals
#include "Quantile.H"
 
using namespace NTC;

int main(int ac, char** av)
{
  double mean, stddev, alpha, beta, p1, p2;
  bool weib = false;
  unsigned long N;
  unsigned n;

  cxxopts::Options options(av[0], "Test the quantile class");

  options.add_options()
   ("h,help", "produce help message")
   ("Normal", "use Normal distribution for sampling")
   ("Weibull", "use Weibull distribution for sampling")
   ("mean", "mean for normal distribution",
    cxxopts::value<double>(mean)->default_value("0.5"))
   ("stddev", "stddev for normal distribution",
    cxxopts::value<double>(stddev)->default_value("2.0"))
   ("alpha", "shape parameter for Weibull distribution",
    cxxopts::value<double>(alpha)->default_value("1.5"))
   ("beta", "scale parameter for Weibull distribution",
    cxxopts::value<double>(beta)->default_value("0.5"))
   ("p1", "first test quantile",
    cxxopts::value<double>(p1)->default_value("0.37"))
   ("p2", "first test quantile",
    cxxopts::value<double>(p2)->default_value("0.65"))
   ("N", "number of samples",
    cxxopts::value<unsigned long>(N)->default_value("100000"))
   ("n", "number of bins",
    cxxopts::value<unsigned>(n)->default_value("100"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 1;
  }

  if (vm.count("Normal")) {
    weib = false;
  }

  if (vm.count("Weibull")) {
    weib = true;
  }

  std::random_device rd;
  std::mt19937 mt(rd());

  std::normal_distribution<double> d(mean, stddev);
  std::weibull_distribution<double> w(alpha, beta);
 
  Quantile quant_1(p1), quant_2(p2), quant_n(n);

  double z, exact1, exact2;

  for (unsigned long n=0; n<N; n++) {
    if (weib) z = w(mt);
    else      z = d(mt);
    quant_1.add(z);
    quant_2.add(z);
    quant_n.add(z);
  }

  boost::math::normal_distribution<double>  nm(mean, stddev);
  boost::math::weibull_distribution<double> wd(alpha, beta);

  if (weib) {
    exact1 = boost::math::quantile(wd, p1);
    exact2 = boost::math::quantile(wd, p2);
  } else {
    exact1 = boost::math::quantile(nm, p1);
    exact2 = boost::math::quantile(nm, p2);
  }

  std::cout << std::left
	    << std::setw(10) << "Quant 1"   << p1         << std::endl
	    << std::setw(10) << "Found"     << quant_1()  << std::endl
	    << std::setw(10) << "Expected " << exact1     << std::endl
	    << std::endl
	    << std::setw(10) << "Quant 2"   << p2         << std::endl
	    << std::setw(10) << "Found"     << quant_2()  << std::endl
	    << std::setw(10) << "Expected " << exact2     << std::endl;

  std::cout << std::endl << "Quant n"
	    << " (compare first two and second two columns)"
	    << std::endl << std::endl
	    << std::setw(18) << "p"
	    << std::setw(18) << "P(<z)"
	    << std::setw(18) << "z=Histo(p)"
	    << std::setw(18) << "Q(p)"
	    << std::endl
	    << std::setw(18) << "----------"
	    << std::setw(18) << "----------"
	    << std::setw(18) << "----------"
	    << std::setw(18) << "----------"
	    << std::endl;

  std::vector<double> qv = {0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95};
  for (auto v : qv) {

    double z = quant_n(v), cdf, ext;

    if (weib) {
      cdf = boost::math::cdf(wd, z);
      ext = boost::math::quantile(wd, v);
    } else {
      cdf = boost::math::cdf(nm, z);
      ext = boost::math::quantile(nm, v);
    }
    
    std::cout << std::setw(18) << v
	      << std::setw(18) << cdf
	      << std::setw(18) << z
	      << std::setw(18) << ext
	      << std::endl;
  }

  return EXIT_SUCCESS;
} 
