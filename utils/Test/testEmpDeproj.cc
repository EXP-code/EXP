#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <cctype>
#include <functional>
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
  int Nr, Ngrid, NumR, Nint;

  // Define pi in a portable way instead of relying on non-standard M_PI
  const double pi = std::acos(-1.0);

  // Parse command-line options
  //
  cxxopts::Options options("testEmpDeproj",
			   "Test the EmpDeproj class against Deproject "
			   "for various surface density profiles.\n");
  options.add_options()
    ("h,help", "Print help")
    ("type", "Surface density type (plummer, gaussian, toomre)", cxxopts::value<std::string>(type)->default_value("toomre"))
    ("abel", "Abel inversion method (derivative, subtraction, ibp)", cxxopts::value<std::string>(abel)->default_value("derivative"))
    ("H", "Scale height for empirical deprojection", cxxopts::value<double>(H)->default_value("0.1"))
    ("Nr", "Number of radial points to evaluate", cxxopts::value<int>(Nr)->default_value("150"))
    ("o,output", "Output file name", cxxopts::value<std::string>(fname)->default_value("rho_test.txt"))
    ("Rmin", "Minimum radius for evaluation", cxxopts::value<double>(Rmin)->default_value("0.01"))
    ("Rmax", "Maximum radius for evaluation", cxxopts::value<double>(Rmax)->default_value("10.0"))
    ("Rcut", "Inner cutoff for donut test", cxxopts::value<double>(Rcut)->default_value("-1.0"))
    ("Rwid", "Width of transition region to inner donut", cxxopts::value<double>(Rwid)->default_value("0.2"))
    ("Ngrid", "Number of grid points for Deprojector", cxxopts::value<int>(Ngrid)->default_value("6000"))
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

  auto it_abel = abel_type_map.find(abel);

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

  auto it = type_map.find(type);

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
    dSigmaFunc = [](double R)->double
    { return -3.0 * R / std::pow(1.0 + R*R, 2.5); };
    // Expected result
    RhoFunc = [pi](double r)->double
    { return 2.0 / std::pow(1.0 + r*r, 2.0) / pi; };
    break;
  case Type::Gaussian:
    // Test function
    SigmaFunc = [](double R)->double
    { return exp(-0.5*R*R); };
    // Analytic derivative
    dSigmaFunc = [](double R)->double
    { return -R*exp(-0.5*R*R); };
    // Expected result
    RhoFunc = [pi](double r)->double
    { return exp(-0.5*r*r)/sqrt(2.0*pi); };
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
		/*R_max_extend=*/50.0, /*tail_power=*/-4.0, /*Ngrid=*/Ngrid);
  
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
    double t = (Nr > 1) ? static_cast<double>(i) / (Nr - 1) : 0.0;
    r_eval.push_back(Rmin + t * (Rmax - Rmin));
  }
  auto rho = D.rho(r_eval);
  
  std::ofstream ofs(fname);
  for (size_t i = 0; i < r_eval.size(); ++i)
    ofs << std::setw(16) << r_eval[i]
	<< std::setw(16) << rho[i]
	<< std::setw(16) << E.density(r_eval[i])
	<< std::setw(16) << RhoFunc(r_eval[i])
	<< std::setw(16) << SigmaFunc(r_eval[i])
	<< std::setw(16) << E.surfaceDensity(r_eval[i])
	<< std::endl;

  ofs.close();
  std::cout << "Wrote " << fname << std::endl;
  
  return 0;
}
