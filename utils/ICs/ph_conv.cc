// This method is stable for all densities, unlike ph_ion.
//

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <string>
#include <cmath>
#include <array>

#include <cxxopts.H>

// Library variables

#include <libvars.H>

/**
   Solution for rate balance.  

   n_H  * f_H   * alpha_H   = n_e * n_H  * (1 - f_H)      * beta_H+
   n_He * f_He  * alpha_He  = n_e * n_He * f_He+          * beta_He+
   n_He * f_He+ * alpha_He+ = n_e * n_He * (1-f_He-f_He+) * beta_He++

   LHS=rate of photoionization
   RHS=rate of recombination

   Simplify to fractional rates:

   f_H   * alpha_H   = n_e * (1 - f_H)      * beta_H+
   f_He  * alpha_He  = n_e * f_He+          * beta_He+
   f_He+ * alpha_He+ = n_e * (1-f_He-f_He+) * beta_He++

   Difference in fraction in time interval h is:

   delta_H   = h*[n_e * (1 - f_H)      * beta_H+   - f_H   * alpha_H]
   delta_He  = h*[n_e * f_He+          * beta_He+  - f_He  * alpha_He]
   delta_He+ = h*[n_e * (1-f_He-f_He+) * beta_He++ - f_He+ * alpha_He+]

   Iterative algorithm:

   f_H  (n+1) = f_H  (n) + delta_H
   f_He (n+1) = f_He (n) + delta_He
   f_He+(n+1) = f_He+(n) + delta_He+

*/
int main(int argc, char**argv)
{
  double n0, tol, h, T, z;
  unsigned skip;
  int niter;
  std::string outf;

  cxxopts::Options options(argv[0], "Compute ionization-recombination equilibrium (stable version)");

  options.add_options()
   ("h,help", "produce this help message")
   ("D,density", "Density in amu/cc. Good for n0<8.5e-2",
     cxxopts::value<double>(n0)->default_value("1.0e-4"))
   ("T,temp", "Density in amu/cc. Good for n0<8.5e-2",
     cxxopts::value<double>(T)->default_value("30000.0"))
   ("s,skip", "Iteration step skip for diagnostic output",
     cxxopts::value<unsigned>(skip)->default_value("10"))
   ("H,step", "Time step in years",
     cxxopts::value<double>(h)->default_value("2000.0"))
   ("e,tol", "error tolerance",
     cxxopts::value<double>(tol)->default_value("1.0e-10"))
   ("z,redshift", "redshift",
     cxxopts::value<double>(z)->default_value("0.1"))
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



  std::map< double, std::array<double, 3> > alpha = {
    {0.1, {6.6362e-14, 7.9561e-14, 8.2955e-15}},
    {1.1, {4.9771e-13, 5.9671e-13, 6.2216e-14}}
  };
     
  if (alpha.find(z) == alpha.end()) {
    std::cout << "Could not find <" << z << "> in alpha" << std::endl
	      << "Available values are:";
    for (auto v : alpha) std::cout << " " << v.first;
    std::cout << std::endl;
    exit(-1);
  }

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

  const double year = 365.25*24.0*3600.0;
  h *= year;

  double X = 0.76, Y = 0.24, mX = 1.0, mY = 4.0;
  double mu = X/mX + Y/mY;

  std::array<double, 3> init  {0.1, 0.1, 0.1};
  std::array<double, 3> curr, last, maxd {0, 0, 0}, frac {X, Y, Y};
  std::array<double, 3> delta;

  curr = init;

  double err = 0.0;

  for (int n=0; n<niter; n++) {

    double ne = n0*(X/mX*(1.0 - curr[0]) +
		    Y/mY*(curr[2] + 2.0*(1.0 - curr[0] - curr[2])));
    last = curr;

    delta[0] = h * ( ne*(1.0 - last[0])           * beta[0] - last[0]*alpha[z][0] );
    delta[1] = h * ( ne*last[2]                   * beta[1] - last[1]*alpha[z][1] );
    delta[2] = h * ( ne*(1.0 - last[1] - last[2]) * beta[2] - last[2]*alpha[z][2] );

    for (int j=0; j<3; j++) {
      curr[j] += delta[j];
      maxd[j] = std::max<double>(fabs(delta[j])/curr[j], maxd[j]);
      curr[j] = std::min<double>(1.0, std::max<double>(0.0, curr[j]));
    }

    if (n % skip==0) {
      std::cout << std::setw(8) << n;
      for (int j=0; j<3; j++) std::cout << std::setw(14) << curr[j];
      std::cout << std::setw(14) << ne/n0 << std::endl;
    }

    err = 0.0;
    for (int j=0; j<3; j++) {
      double dif = 0.5 * (curr[j] - last[j]) / (curr[j] + last[j]);
      err += dif*dif;
    }
    err = sqrt(err);
    if (err<tol) break;
  }

  std::cout << std::endl << std::left
	    << std::setw(24) << "Convergence error"   << err << std::endl
	    << std::setw(24) << "Requested tolerance" << tol << std::endl
	    << std::setw(24) << "Max relative error";
  for (auto v : maxd) std::cout << std::setw(14) << v;
  std::cout << std::endl << std::endl;

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
