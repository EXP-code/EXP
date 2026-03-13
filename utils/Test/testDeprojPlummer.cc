#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "Deprojector.H"

using namespace Deproject;

int main()
{
  std::function<double(double)> SigmaFunc, dSigmaFunc, RhoFunc;

  enum class Type { Plummer, Gaussian, Toomre };
  Type which = Type::Plummer;

  switch (which) {
  case Type::Toomre:
    // Test function
    SigmaFunc = [](double R)->double
    { return 1.0 / std::pow(1.0 + R*R, 1.5); };
    // Analytic derivative
    dSigmaFunc = [](double R)->double
    { return -3.0 * R / std::pow(1.0 + R*R, 2.5); };
    // Expected result
    RhoFunc = [](double r)->double
    { return 2.0 / std::pow(1.0 + r*r, 2.0) / M_PI; };
    break;
  case Type::Gaussian:
    // Test function
    SigmaFunc = [](double R)->double
    { return exp(-0.5*R*R); };
    // Analytic derivative
    dSigmaFunc = [](double R)->double
    { return -R*exp(-0.5*R*R); };
    // Expected result
    RhoFunc = [](double r)->double
    { return exp(-0.5*r*r)/sqrt(2.0*M_PI); };
    break;
  default:
  case Type::Plummer:
    // Test function
    SigmaFunc = [](double R)->double
    { return 4.0 / 3.0 / std::pow(1.0 + R*R, 2.0); };
    // Analytic derivative
    dSigmaFunc = [](double R)->double
    { return -16.0 * R / 3.0 / std::pow(1.0 + R*R, 3.0); };
    // Expected result
    RhoFunc = [](double r)->double
    { return 1.0 / std::pow(1.0 + r*r, 2.5); };
    break;
  }
  
  Deprojector D(SigmaFunc, dSigmaFunc, /*R_data_min=*/0.01, /*R_data_max=*/10.0,
		/*R_max_extend=*/50.0, /*tail_power=*/-4.0, /*Ngrid=*/6000);
  
  std::vector<double> r_eval;
  int Nr = 150;
  for (int i = 0; i < Nr; ++i) {
    double t = (double)i / (Nr - 1);
    r_eval.push_back(0.01 + t * 8.0);
  }
  auto rho = D.rho(r_eval);
  
  std::ofstream ofs("rho_test.txt");
  for (size_t i = 0; i < r_eval.size(); ++i) {
    ofs << std::setw(16) << r_eval[i]
	<< std::setw(16) << rho[i]
	<< std::setw(16) << RhoFunc(r_eval[i])
	<< std::endl;
  }
  ofs.close();
  std::cout << "Wrote rho_test.txt\n";
  
  return 0;
}
