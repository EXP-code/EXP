#ifndef _DiskDensityFunc_H_
#define _DiskDensityFunc_H_

#include <functional>

#include <config_exp.h>

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

/** Python disk-density function wrapper

    The wrapper class will take a user supplied Python module that
    must supply the 'disk-density' function

    An example disk-density Python module:
    --------------------------cut here--------------------------------
    import math

    def disk_density(R, z, phi):
        """Computes the disk density at a given radius, height, and
        azimuth.  This would be the user's target density for basis
        creation.

        """
        a = 1.0                            # Scale radius
        h = 0.2                            # Scale height
        f = math.exp(-0.5*math.fabs(z)/h)  # Prevent overflows
        sech = 2.0*f / (1.0 + f*f)         # 
        return math.exp(-R/a)*sech*sech/(8*math.pi*h*a*a)
	
    --------------------------cut here--------------------------------

 */
class __attribute__((visibility("default")))
DiskDensityFunc
{
private:

#ifdef HAVE_PYTHON3
  //! The disk-density function
  pybind11::function disk_density;
#endif

  //! Interpreter started?
  bool started = false;

  //! Name of density function
  std::string funcname;

public:

  //! Constructor
  DiskDensityFunc(const std::string& modulename,
		  const std::string& funcname="disk_density");

  //! Destructor
  ~DiskDensityFunc();

  //! Evaluate the density
  double operator() (double R, double z, double phi);

};

#endif
