#include <stdexcept>

#include "CubicSpline.H"

namespace Deproject
{
  
  CubicSpline::CubicSpline(const std::vector<double>& x_in, const std::vector<double>& y_in) {
    set_data(x_in, y_in);
  }
  
  void CubicSpline::set_data(const std::vector<double>& x_in, const std::vector<double>& y_in) {
    x_ = x_in;
    y_ = y_in;
    int n = (int)x_.size();
    if (n < 2 || (int)y_.size() != n) throw std::runtime_error("CubicSpline: need at least two points and equal-sized x,y.");
    
    y2_.assign(n, 0.0);
    std::vector<double> u(n - 1, 0.0);
    
    // natural spline boundary conditions (second derivatives at endpoints = 0)
    for (int i = 1; i < n - 1; ++i) {
      double sig = (x_[i] - x_[i-1]) / (x_[i+1] - x_[i-1]);
      double p = sig * y2_[i-1] + 2.0;
      y2_[i] = (sig - 1.0) / p;
      double dY1 = (y_[i+1] - y_[i]) / (x_[i+1] - x_[i]);
      double dY0 = (y_[i]   - y_[i-1]) / (x_[i]   - x_[i-1]);
      u[i] = (6.0 * (dY1 - dY0) / (x_[i+1] - x_[i-1]) - sig * u[i-1]) / p;
    }
    
    for (int k = n - 2; k >= 0; --k) y2_[k] = y2_[k] * y2_[k+1] + u[k];
  }
  
  int CubicSpline::locate(double xx) const {
    int n = (int)x_.size();
    if (xx <= x_.front()) return 0;
    if (xx >= x_.back()) return n - 2;
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      int mid = (lo + hi) >> 1;
      if (x_[mid] > xx) hi = mid; else lo = mid;
    }
    return lo;
  }
  
  double CubicSpline::eval(double xx) const {
    int klo = locate(xx);
    int khi = klo + 1;
    double h = x_[khi] - x_[klo];
    if (h <= 0.0) throw std::runtime_error("CubicSpline::eval: non-increasing x.");
    double A = (x_[khi] - xx) / h;
    double B = (xx - x_[klo]) / h;
    double val = A * y_[klo] + B * y_[khi]
      + ((A*A*A - A) * y2_[klo] + (B*B*B - B) * y2_[khi]) * (h*h) / 6.0;
    return val;
  }
  
  double CubicSpline::deriv(double xx) const {
    int klo = locate(xx);
    int khi = klo + 1;
    double h = x_[khi] - x_[klo];
    if (h <= 0.0) throw std::runtime_error("CubicSpline::deriv: non-increasing x.");
    double A = (x_[khi] - xx) / h;
    double B = (xx - x_[klo]) / h;
    double dy = (y_[khi] - y_[klo]) / h
      + ( - (3.0*A*A - 1.0) * y2_[klo] + (3.0*B*B - 1.0) * y2_[khi] ) * (h / 6.0);
    return dy;
  }
  
  double CubicSpline::xmin() const { return x_.front(); }
  double CubicSpline::xmax() const { return x_.back(); }

}
