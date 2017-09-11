#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>

//
// BOOST stuff
//
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char**argv)
{
  double n0, tol;
  int niter;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce this help message")
    ("density,D", po::value<double>(&n0)->default_value(1.0e-4), 
     "Density in amu/cc. Good for n0<8.5e-2")
    ("tol,e",     po::value<double>(&tol)->default_value(1.0e-10), 
     "error tolerance")
    ("iter,n",    po::value<int>(&niter)->default_value(1000), 
     "maximum number of iterations")
    ;
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  std::array<double, 3> alpha {4.9771e-13, 5.9671e-13, 6.2216e-14};
  std::array<double, 3> beta  {1.6449e-13, 1.8229e-13, 8.3322e-13};
  std::array<double, 3> gamma;
  std::array<double, 3> init  {0.001, 0.001, 0.001};
  std::array<double, 3> curr, last;

  for (int i=0; i<3; i++) gamma[i] = beta[i]/alpha[i];
  double X = 0.76, Y = 0.24, mX = 1.0, mY = 4.0;
  double mu = X/mX + Y/mY;

  curr = init;

  for (int n=0; n<niter; n++) {
    double ne = n0*(X/mX*(1.0 - curr[0]) +
		    Y/mY*(curr[2] + 2.0*(1.0 - curr[0] - curr[2])));
    last = curr;
    curr[0] = ne*(1.0 - last[0])*gamma[0];
    curr[1] = ne*last[2]*gamma[1];
    curr[2] = ne*(1.0 - last[1] - last[2])*gamma[2];

    std::cout << std::setw(8) << n;
    for (int j=0; j<3; j++) std::cout << std::setw(14) << curr[j];
    std::cout << std::endl;

    double tol = 0.0;
    for (int j=0; j<3; j++)
      tol += (curr[j] - last[j]) * (curr[j] - last[j]);
    tol = sqrt(tol);
    if (tol<1.0e-14) break;
  }

  return 0;
}
