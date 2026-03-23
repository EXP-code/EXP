#pragma once

#include <functional>
#include <cassert>
#include <cmath>
#include "interp.H"

class EmpDeproj
{
private:
  //! Interpolators for density and mass
  Linear1d densRg, massRg, surfRg;

public:
  //! Abel type
  enum class AbelType { Derivative, Subtraction, IBP };

  //! Construct from analytic Sigma functor
  EmpDeproj(double H, double Rmin, double Rmax, int NUMR, int NINT,
	    std::function<double(double, double)> func,
	    AbelType type = AbelType::Derivative);

  //! Destructor
  ~EmpDeproj() {}

  //! Evaluate deprojected density at a single radius
  double density(double R) const
  {
    assert(R > 0.0);
    return densRg.eval(std::log(R));
  }

  double surfaceDensity(double R) const
  {
    assert(R > 0.0);
    return surfRg.eval(std::log(R));
  }

  //! Evaluate deprojected mass at a single radius
  double mass(double R) const
  {
    assert(R > 0.0);
    return massRg.eval(std::log(R));
  }
};
