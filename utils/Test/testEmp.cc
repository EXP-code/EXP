#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <filesystem>
#include <vector>
#include <memory>
#include <cmath>
#include <map>

#include "EmpDeproj.H"
#include "cxxopts.H"

// pybind11 embedding
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
namespace py = pybind11;

int main(int argc, char* argv[])
{
  // Default parameters
  std::string type_opt;
  std::string abel_opt;
  std::string fname;
  double H = 0.1, Rmin = 0.01, Rmax = 10.0, Rcut = -1.0, Rwid = 0.2;
  int Nr = 150, NumR = 1000, Nint = 800;

  // Parse command-line options
  cxxopts::Options options("testDonut",
               "Test the EmpDeproj class for an inner donut-shaped "
               "density distribution, using the Toomre profile as "
               "a test case.");

  options.add_options()
    ("h,help", "Print help")
    ("type", "Surface density type (plummer, gaussian, toomre) - ignored if Python function supplied",
       cxxopts::value<std::string>(type_opt)->default_value("toomre"))
    ("abel", "Abel inversion method (derivative, subtraction, ibp)",
       cxxopts::value<std::string>(abel_opt)->default_value("derivative"))
    ("H", "Scale height for empirical deprojection", cxxopts::value<double>(H)->default_value("0.1"))
    ("Nr", "Number of radial points to evaluate", cxxopts::value<int>(Nr)->default_value("150"))
    ("o,output", "Output file name", cxxopts::value<std::string>(fname)->default_value("rho_test.txt"))
    ("Rmin", "Minimum radius for evaluation", cxxopts::value<double>(Rmin)->default_value("0.01"))
    ("Rmax", "Maximum radius for evaluation", cxxopts::value<double>(Rmax)->default_value("10.0"))
    ("Rcut", "Inner cutoff for donut test", cxxopts::value<double>(Rcut)->default_value("-1.0"))
    ("Rwid", "Width of transition region to inner donut", cxxopts::value<double>(Rwid)->default_value("0.2"))
    ("NumR", "Number of radial points for EmpDeproj", cxxopts::value<int>(NumR)->default_value("1000"))
    ("Nint", "Number of integration points for EmpDeproj", cxxopts::value<int>(Nint)->default_value("800"))
    // Python integration options
    ("pymodule", "Python module name OR path to a .py file containing a function", cxxopts::value<std::string>())
    ("pyfunc", "Function name inside module/file (default: 'Sigma')", cxxopts::value<std::string>()->default_value("Sigma"));

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }

  // Map Abel type string to enum
  EmpDeproj::AbelType type_enum;
  std::map<std::string, EmpDeproj::AbelType> abel_type_map = {
    {"derivative", EmpDeproj::AbelType::Derivative},
    {"subtraction", EmpDeproj::AbelType::Subtraction},
    {"ibp", EmpDeproj::AbelType::IBP}
  };

  std::string abel_l = result["abel"].as<std::string>();
  std::transform(abel_l.begin(), abel_l.end(), abel_l.begin(),
                 [](unsigned char c){ return std::tolower(c); });

  auto it_abel = abel_type_map.find(abel_l);
  if (it_abel != abel_type_map.end()) {
    type_enum = it_abel->second;
  } else {
    throw std::runtime_error("Unknown Abel type: " + result["abel"].as<std::string>());
  }

  // Prepare SigmaZFunc to pass into EmpDeproj
  std::function<double(double,double)> SigmaZFunc;

  // For pybind11 embedding, we need to keep the interpreter alive for
  // the duration of the SigmaZFunc usage. We can use a unique_ptr to
  // manage the lifetime of the interpreter guard. If no Python is
  // used, this will just be an empty guard that does nothing.
  //
  std::unique_ptr<py::scoped_interpreter> pyguard;

  // If user supplied a Python module/file and function, embed Python
  // and load it here
  //
  if (result.count("pymodule")) {
    std::string pymod = result["pymodule"].as<std::string>();
    std::string pyfuncname = result["pyfunc"].as<std::string>();

    // Start the Python interpreter
    pyguard = std::make_unique<py::scoped_interpreter>();

    py::object py_module;

    // If pymod ends with .py, treat it as filepath: insert its
    // directory to sys.path and import by stem
    std::filesystem::path p(pymod);
    try {
      if (p.has_extension() && p.extension() == ".py") {
        std::string dir = p.parent_path().string();
        std::string modname = p.stem().string();
        if (!dir.empty()) {
          py::module_ sys = py::module_::import("sys");
          // insert at front so local dir is found first
          sys.attr("path").attr("insert")(0, dir);
        }
        py_module = py::module_::import(modname.c_str());
      } else {
        // treat as module name
        py_module = py::module_::import(pymod.c_str());
      }
    } catch (const py::error_already_set &e) {
      throw std::runtime_error(std::string("Failed to import Python module '")
			       + pymod + "': " + e.what());
    }

    py::object pyfunc = py_module.attr(pyfuncname.c_str());

    if (!py::isinstance<py::function>(pyfunc) && !py::hasattr(pyfunc, "__call__")) {
      throw std::runtime_error("Python object " + pyfuncname +
			       " is not callable.");
    }

    // Inspect function argument count to decide whether it's Sigma(R)
    // or SigmaZ(R,z).  This is probably overkill and might not be
    // robust for all callables, but it allows some flexibility for
    // users.
    //
    int argcount = 0;
    try {
      // functions have __code__.co_argcount; builtins may not — in
      // that case prefer calling and checking.  I'm still not sure if
      // this is the best way to do it, but it should cover most cases
      // (functions, builtins, callables).  Python experts?
      //
      if (py::hasattr(pyfunc, "__code__")) {
        argcount = pyfunc.attr("__code__").attr("co_argcount").cast<int>();
      } else if (py::hasattr(pyfunc, "__call__") && py::hasattr(pyfunc.attr("__call__"), "__code__")) {
        argcount = pyfunc.attr("__call__").attr("__code__").attr("co_argcount").cast<int>() - 1; // bound method
      } else {
        // fallback, try calling with 2 args first and catch
        argcount = -1;
      }
    } catch (...) {
      argcount = -1;
    }

    if (argcount == 1) {
      // User provided Sigma(R). Wrap it and add vertical profile +
      // hole logic here.  Keep a copy of pyfunc alive by capturing
      // py::object by value.
      std::function<double(double)> Sigma = [pyfunc](double R)->double {
        py::gil_scoped_acquire gil;
        py::object out = pyfunc(R);
        return out.cast<double>();
      };

      SigmaZFunc = [Sigma, H, Rcut, Rwid](double R, double z)->double {
        // vertical profile: sech^2 with scale H (same form as original)
        double Q = std::exp(-std::fabs(z) / (2.0 * H));
        double sech = 2.0 * Q / (1.0 + Q * Q);
        double hole = 1.0;
        if (Rcut > 0.0) {
          double x = (R - Rcut) / Rwid;
          hole = 0.5 * (1.0 + std::tanh(x));
        }
        double s = Sigma(R);
        return s * sech * sech * hole / (4.0 * H);
      };

    } else if (argcount == 2) {
      // User provided SigmaZ(R,z) directly. Use it as-is (no extra
      // vertical/hole logic).
      py::object pyfunc2 = pyfunc;
      SigmaZFunc = [pyfunc2](double R, double z)->double {
        py::gil_scoped_acquire gil;
        py::object out = pyfunc2(R, z);
        return out.cast<double>();
      };

    } else {
      // ambiguous: try calling with 2 args; if that fails try 1 arg
      try {
        // test call with dummy values
        py::gil_scoped_acquire gil;
        pyfunc(1.0, 0.0);
        // succeeded: treat as SigmaZ(R,z)
        py::object pyfunc2 = pyfunc;
        SigmaZFunc = [pyfunc2](double R, double z)->double {
          py::gil_scoped_acquire gil2;
          py::object out = pyfunc2(R, z);
          return out.cast<double>();
        };
      } catch (const py::error_already_set &) {
        // fallback: try as Sigma(R)
        py::object pyfunc1 = pyfunc;
        std::function<double(double)> Sigma = [pyfunc1](double R)->double {
          py::gil_scoped_acquire gil2;
          py::object out = pyfunc1(R);
          return out.cast<double>();
        };
        SigmaZFunc = [Sigma, H, Rcut, Rwid](double R, double z)->double {
          double Q = std::exp(-std::fabs(z) / (2.0 * H));
          double sech = 2.0 * Q / (1.0 + Q * Q);
          double hole = 1.0;
          if (Rcut > 0.0) {
            double x = (R - Rcut) / Rwid;
            hole = 0.5 * (1.0 + std::tanh(x));
          }
          double s = Sigma(R);
          return s * sech * sech * hole / (4.0 * H);
        };
      }
    }
  } // end python handling
  else {
    // No Python supplied: use internal choices (plummer/gaussian/toomre)
    std::function<double(double)> SigmaFunc;
    enum class Type { Plummer, Gaussian, Toomre };
    Type which = Type::Toomre;

    std::map<std::string, Type> type_map = {
      {"plummer", Type::Plummer},
      {"gaussian", Type::Gaussian},
      {"toomre", Type::Toomre}
    };

    std::string type_l = result["type"].as<std::string>();
    std::transform(type_l.begin(), type_l.end(), type_l.begin(),
             [](unsigned char c){ return std::tolower(c); });

    auto it = type_map.find(type_l);
    if (it != type_map.end()) {
      which = it->second;
    } else {
      throw std::runtime_error("Unknown type: " + result["type"].as<std::string>());
    }

    switch (which) {
    case Type::Toomre:
      SigmaFunc = [](double R)->double { return 1.0 / std::pow(1.0 + R*R, 1.5); };
      break;
    case Type::Gaussian:
      SigmaFunc = [](double R)->double { return std::exp(-0.5*R*R); };
      break;
    default:
    case Type::Plummer:
      SigmaFunc = [](double R)->double { return 4.0 / 3.0 / std::pow(1.0 + R*R, 2.0); };
      break;
    }

    // Build SigmaZFunc from SigmaFunc, using the same vertical/hole
    // logic
    SigmaZFunc = [SigmaFunc, H, Rcut, Rwid](double R, double z)->double {
      double Q = std::exp(-std::fabs(z) / (2.0 * H));
      double sech = 2.0 * Q / (1.0 + Q * Q);
      double hole = 1.0;
      if (Rcut > 0.0) {
        double x = (R - Rcut) / Rwid;
        hole = 0.5 * (1.0 + std::tanh(x));
      }
      return SigmaFunc(R) * sech * sech * hole / (4.0 * H);
    };
  } // end else internal choice

  // Construct EmpDeproj and evaluate
  EmpDeproj E(H, Rmin, Rmax, NumR, Nint, SigmaZFunc, type_enum);

  // radial evaluation points
  std::vector<double> r_eval;
  for (int i = 0; i < Nr; ++i) {
    double t = (double)i / (Nr - 1);
    r_eval.push_back(Rmin + t * (Rmax - Rmin));
  }

  std::ofstream ofs(fname);
  if (!ofs) {
    std::cerr << "Failed to open output file: " << fname << std::endl;
    return 1;
  }

  for (size_t i = 0; i < r_eval.size(); ++i) {
    ofs << std::setw(16) << r_eval[i]
        << std::setw(16) << E.density(r_eval[i])
        << std::setw(16) << E.surfaceDensity(r_eval[i])
        << std::endl;
  }

  ofs.close();
  std::cout << "Wrote " << fname << std::endl;

  return 0;
}
