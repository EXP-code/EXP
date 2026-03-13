#pragma once

#include <functional>
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
  double density(double R)
  {
    return densRg.eval(log(R));
  }

  double surfaceDensity(double R)
  {
    return surfRg.eval(log(R));
  }

  //! Evaluate deprojected mass at a single radius
  double mass(double R)
  {
    return massRg.eval(log(R));
  }
};
