#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "Deprojector.H"
#include "EmpDeproj.H"
#include "cxxopts.H"

using namespace Deproject;

int main(int argc, char* argv[])
{

  // Parameters
  //
  std::string type, abel, fname;
  double H, Rmin, Rmax, Rcut, Rwid;
  int Nr, NumR, Nint;

  // Parse command-line options
  //
  cxxopts::Options options("testDonut",
			   "Test the EmpDeproj class for an inner donut-shaped "
			   "density distribution, using the Toomre profile as "
			   "a test case.");

  options.add_options()
    ("h,help", "Print help")
    ("type", "Surface density type (plummer, gaussian, toomre)", cxxopts::value<std::string>()->default_value("toomre"))
    ("abel", "Abel inversion method (derivative, subtraction, ibp)", cxxopts::value<std::string>()->default_value("derivative"))
    ("H", "Scale height for empirical deprojection", cxxopts::value<double>(H)->default_value("0.1"))
    ("Nr", "Number of radial points to evaluate", cxxopts::value<int>(Nr)->default_value("150"))
    ("o,output", "Output file name", cxxopts::value<std::string>(fname)->default_value("rho_test.txt"))
    ("Rmin", "Minimum radius for evaluation", cxxopts::value<double>(Rmin)->default_value("0.01"))
    ("Rmax", "Maximum radius for evaluation", cxxopts::value<double>(Rmax)->default_value("10.0"))
    ("Rcut", "Inner cutoff for donut test", cxxopts::value<double>(Rcut)->default_value("-1.0"))
    ("Rwid", "Width of transition region to inner donut", cxxopts::value<double>(Rwid)->default_value("0.2"))
    ("NumR", "Number of radial points for EmpDeproj", cxxopts::value<int>(NumR)->default_value("1000"))
    ("Nint", "Number of integration points for EmpDeproj", cxxopts::value<int>(Nint)->default_value("800"));

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }

  // Map Abel type string to enum
  //
  EmpDeproj::AbelType type_enum;
  std::map<std::string, EmpDeproj::AbelType> abel_type_map = {
    {"derivative", EmpDeproj::AbelType::Derivative},
    {"subtraction", EmpDeproj::AbelType::Subtraction},
    {"ibp", EmpDeproj::AbelType::IBP}
  };

  // Convert type string to lower case
  //
  std::transform(abel.begin(), abel.end(), abel.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  auto it_abel = abel_type_map.find(result["abel"].as<std::string>());

  if (it_abel != abel_type_map.end()) {
    type_enum = it_abel->second;
  } else {
    throw std::runtime_error("Unknown Abel type: " + result["abel"].as<std::string>());
  }

  std::function<double(double)> SigmaFunc, dSigmaFunc, RhoFunc;
  enum class Type { Plummer, Gaussian, Toomre };
  Type which = Type::Toomre;

  std::map<std::string, Type> type_map = {
    {"plummer", Type::Plummer},
    {"gaussian", Type::Gaussian},
    {"toomre", Type::Toomre}
  };

  // Convert type string to lower case
  //
  std::transform(type.begin(), type.end(), type.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  auto it = type_map.find(result["type"].as<std::string>());

  if (it != type_map.end()) {
    which = it->second;
  } else {
    throw std::runtime_error("Unknown type: " + result["type"].as<std::string>());
  }

  switch (which) {
  case Type::Toomre:
    // Test function
    SigmaFunc = [](double R)->double
    { return 1.0 / std::pow(1.0 + R*R, 1.5); };
    // Analytic derivative
    break;
  case Type::Gaussian:
    // Test function
    SigmaFunc = [](double R)->double
    { return exp(-0.5*R*R); };
    // Analytic derivative
    break;
  default:
  case Type::Plummer:
    // Test function
    SigmaFunc = [](double R)->double
    { return 4.0 / 3.0 / std::pow(1.0 + R*R, 2.0); };
    break;
  }
  
  auto SigmaZFunc = [SigmaFunc, H, Rcut, Rwid](double R, double z)->double
  { double Q = exp(-std::fabs(z)/(2.0*H));
    double sech = 2.0*Q / (1.0 + Q*Q);
    double hole = 1.0;
    if (Rcut > 0.0) {
      double x = (R - Rcut) / Rwid;
      hole = 0.5 * (1.0 + std::tanh(x));
    }
    return SigmaFunc(R)*sech*sech*hole/(4.0*H); };

  EmpDeproj E(H, Rmin, Rmax, NumR, Nint, SigmaZFunc, type_enum);

  std::vector<double> r_eval;
  for (int i = 0; i < Nr; ++i) {
    double t = (double)i / (Nr - 1);
    r_eval.push_back(0.01 + t * 8.0);
  }
  
  std::ofstream ofs(fname);
  for (size_t i = 0; i < r_eval.size(); ++i)
    ofs << std::setw(16) << r_eval[i]
	<< std::setw(16) << E.density(r_eval[i])
	<< std::setw(16) << E.surfaceDensity(r_eval[i])
	<< std::endl;

  ofs.close();
  std::cout << "Wrote " << fname << std::endl;
  
  return 0;
}
