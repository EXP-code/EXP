#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include <vector>

namespace Deproject
{

  class CubicSpline {
  public:
    CubicSpline() = default;
    CubicSpline(const std::vector<double>& x_in, const std::vector<double>& y_in);
    
    // set data (x must be strictly increasing)
    void set_data(const std::vector<double>& x_in, const std::vector<double>& y_in);
    
    // evaluate spline and its derivative (xx is expected in [xmin(), xmax()]; values outside are extrapolated using the end intervals)
    double eval(double xx) const;
    double deriv(double xx) const;
    
    double xmin() const;
    double xmax() const;
    
  private:
    std::vector<double> x_, y_, y2_;
    int locate(double xx) const;
  };
  
}

#endif // CUBIC_SPLINE_H
