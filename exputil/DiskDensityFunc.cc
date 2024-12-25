#include <DiskDensityFunc.H>

namespace py = pybind11;

DiskDensityFunc::DiskDensityFunc(const std::string& modulename)
{
  if (Py_IsInitialized() == 0) {
    py::initialize_interpreter();
    started = true;
  }

  disk_density =
    py::reinterpret_borrow<py::function>
    ( py::module::import(modulename.c_str()).attr("disk_density") );
}

DiskDensityFunc::~DiskDensityFunc()
{
  if (started) py::finalize_interpreter();
}

double DiskDensityFunc::operator() (double R, double z, double phi)
{
  return disk_density(R, z, phi).cast<double>();
}
