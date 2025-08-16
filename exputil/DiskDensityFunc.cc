#include <DiskDensityFunc.H>

#ifdef HAVE_PYTHON3

namespace py = pybind11;

DiskDensityFunc::DiskDensityFunc(const std::string& modulename,
				 const std::string& funcname)
  : funcname(funcname)
{
  // Check if a Python interpreter exists
  if (Py_IsInitialized() == 0) {
    py::initialize_interpreter();
    // Mark the interpreter as started by this instance
    started = true;
  }

  // Bind the disk_density function from Python
  disk_density =
    py::reinterpret_borrow<py::function>
    ( py::module::import(modulename.c_str()).attr(funcname.c_str()) );
}

DiskDensityFunc::~DiskDensityFunc()
{
  // Only end the interpreter if it was started by this instance
  if (started) py::finalize_interpreter();
}

double DiskDensityFunc::operator() (double R, double z, double phi)
{
  return disk_density(R, z, phi).cast<double>();
}

#else

DiskDensityFunc::DiskDensityFunc(const std::string& modulename,
				 const std::string& funcname)
  : funcname(funcname)
{
  throw std::runtime_error("DiskDensityFunc: you environoment does not have Python3 support.  Use a built-in density target or install Python3 and recompile");
}

DiskDensityFunc::~DiskDensityFunc()
{
}

double DiskDensityFunc::operator() (double R, double z, double phi)
{
  return 0.0;
}

#endif
