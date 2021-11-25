#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <cmath>
#include <array>

// Library variables

#include <libvars.H>
#include <cxxopts.H>

int main(int argc, char**argv)
{
  double n0, tol;
  unsigned T;
  int niter;
  std::string outf;

  cxxopts::Options options(argv[0], "Test of ionization equilibrium");

  options.add_options()
   ("h,help", "produce this help message")
   ("D,density", "Density in amu/cc. Good for n0<8.5e-2",
     cxxopts::value<double>(n0)->default_value("1.0e-4"))
   ("T,temp", "Temperature in Kelvin (integer value)",
     cxxopts::value<unsigned>(T)->default_value("25000"))
   ("e,tol", "error tolerance",
     cxxopts::value<double>(tol)->default_value("1.0e-10"))
   ("n,iter", "maximum number of iterations",
     cxxopts::value<int>(niter)->default_value("1000"))
   ("o,outfile", "data file for makeIon input",
     cxxopts::value<std::string>(outf)->default_value("IonRecombFrac.data"))
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 1;
  }

  std::array<double, 3> alpha {4.9771e-13, 5.9671e-13, 6.2216e-14};

  std::vector<double> Temp;
  std::vector<std::array<double, 3>> barray;
  std::array<double, 3> beta;

  std::ifstream fin("coefficients.dat");
  if (fin) {
    double t;
    std::array<double, 3> v;
    while (fin) {
      fin >> t >> v[0] >> v[1] >> v[2];
      if (fin.good()) {
	Temp.push_back(t);
	barray.push_back(v);
      }
    }
  } else {
    std::cout << "Error opening <coefficients.dat>" << std::endl;
    exit(-1);
  }

  // Interpolate
  if (T<=Temp.back() and T>=Temp.front()) {
    unsigned indx = 0;
    for (indx = 0; indx < Temp.size(); indx++) {
      if (T<Temp[indx]) break;
    }
    if (indx == 0) indx = 1;
    if (indx == Temp.size()) indx = Temp.size()-2;
    double a = (T - Temp[indx-1])/(Temp[indx] - Temp[indx-1]);
    double b = (Temp[indx] - T)  /(Temp[indx] - Temp[indx-1]);
    std::cout << "[a, b] = [" << a << ", " << b << "]" << std::endl;
    std::cout << "[T_1, T_2] = [" << Temp[indx-1] << ", " << Temp[indx] << "]" << std::endl;
    for (int k=0; k<3; k++) beta[k] = a*barray[indx-1][k] + b*barray[indx][k];
  } else {
    std::cout << "T=" << T << " is out of bounds [" << Temp.front()
	      << ", " << Temp.back() << "]" << std::endl;
    exit(-1);
  }

  std::array<double, 3> gamma;
  std::array<double, 3> init  {0.001, 0.001, 0.001};
  std::array<double, 3> curr, last;

  for (int i=0; i<3; i++) gamma[i] = beta[i]/alpha[i];
  double X = 0.76, Y = 0.24, mX = 1.0, mY = 4.0;
  double mu = X/mX + Y/mY;

  curr = init;

  double err = 0.0;

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

    err = 0.0;
    for (int j=0; j<3; j++)
      err += (curr[j] - last[j]) * (curr[j] - last[j]);
    err = sqrt(err);
    if (err < tol) break;
  }

  std::cout << std::endl << std::left
	    << std::setw(24) << "Convergence error"   << err << std::endl
	    << std::setw(24) << "Requested tolerance" << tol << std::endl
	    << std::endl;

  if (err < tol) {
    std::ofstream out(outf);
    if (out) {
      for (int j=0; j<3; j++) out << std::setw(14) << curr[j];
      out << std::endl;
      std::cout << "SUCCESS: "
		<< "file <" << outf << "> written" << std::endl;
    } else {
      std::cout << "FAILURE: "
		<< "error opening <" << outf << "> for output" << std::endl;
    }
  } else {
    std::cout << "FAILURE: no convergence" << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
