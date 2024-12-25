#include <DiskDensityFunc.H>

namespace py = pybind11;	// Convenience


DiskDensityFunc::DiskDensityFunc(const std::string& modulename)
{
  py::initialize_interpreter();

  disk_density =
    py::reinterpret_borrow<py::function>
    ( py::module::import(modulename.c_str()).attr("disk_density") );
}

DiskDensityFunc::~DiskDensityFunc()
{
  py::finalize_interpreter();
}

double DiskDensityFunc::operator() (double R, double z, double phi)
{
  return disk_density(R, z, phi).cast<double>();
}
