#ifndef _DiskDensityFunc_H_
#define _DiskDensityFunc_H_

#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

// Convenience
namespace py = pybind11;

//! Python disk-density function wrapper
class DiskDensityFunc
{
private:

  //! The disk-density function
  py::function disk_density;

public:

  //! Constructor
  DiskDensityFunc(const std::string& modulename);

  //! Destructor
  ~DiskDensityFunc();

  //! Evaluate the density
  double operator() (double R, double z, double phi);

};

#endif
