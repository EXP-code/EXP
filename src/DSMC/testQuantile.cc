// Compile string: mpiCC -std=c++11 -O3 -o testQ testQuantile.cc Quantile.cc -lboost_program_options -lmpi

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <random>

#include <boost/program_options.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/normal.hpp>

#include "Quantile.H"
 
// Unused globals needed in libexputil
//
int myid;
std::mt19937 random_gen;
std::string outdir, runtag;

namespace po = boost::program_options;
using namespace NTC;

int main(int ac, char** av)
{
  double mean, stddev, alpha, beta, p1, p2;
  bool weib = false;
  unsigned long N;
  unsigned n;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("Normal",          "use Normal distribution for sampling")
    ("Weibull",         "use Weibull distribution for sampling")
    ("mean",            po::value<double>(&mean)->default_value(0.5),
     "mean for normal distribution")
    ("stddev",          po::value<double>(&stddev)->default_value(2.0),
     "stddev for normal distribution")
    ("alpha",           po::value<double>(&alpha)->default_value(1.5),
     "shape parameter for Weibull distribution")
    ("beta",            po::value<double>(&beta)->default_value(0.5),
     "scale parameter for Weibull distribution")
    ("p1",            po::value<double>(&p1)->default_value(0.37),
     "first test quantile")
    ("p2",            po::value<double>(&p2)->default_value(0.65),
     "first test quantile")
    ("N",               po::value<unsigned long>(&N)->default_value(100000),
     "number of samples")
    ("n",               po::value<unsigned>(&n)->default_value(100),
     "number of bins")
    ;

  po::variables_map vm;

  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    return -1;
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
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
