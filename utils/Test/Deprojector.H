#ifndef DEPROJECTOR_H
#define DEPROJECTOR_H

#include <functional>
#include <vector>

#include "CubicSpline.H"

namespace Deproject
{
  
  class Deprojector {
  public:
    // Construct from analytic Sigma functor. If dSigmaFunc is empty, numerical derivative is used.
    // R_data_min and R_data_max define where SigmaFunc is trusted for tail matching.
    Deprojector(std::function<double(double)> SigmaFunc,
                std::function<double(double)> dSigmaFunc,
                double R_data_min,
                double R_data_max,
                double R_max_extend = -1.0,
                double tail_power = -3.0,
                int Ngrid = 4000);
    
    // Construct from sampled data arrays (R must be positive and will be sorted)
    Deprojector(const std::vector<double>& R_in,
                const std::vector<double>& Sigma_in,
                double R_max_extend = -1.0,
                double tail_power = -3.0,
                int Ngrid = 4000);
    
    // Evaluate rho at a single radius
    double rho_at(double r) const;
    
    // Evaluate rho for a vector of r
    std::vector<double> rho(const std::vector<double>& r_eval) const;
    
  private:
    void build_grid(int Ngrid);
    
    std::function<double(double)> sigma_func_;
    std::function<double(double)> dsigma_func_; // may be empty
    CubicSpline spline_; // used when constructed from sampled data
    
    double Rdata_min_, Rdata_max_;
    double Rmin_, Rmax_;
    double tail_power_;
    
    std::vector<double> dx_;	// grid spacings
  
    std::vector<double> fineR_;
    std::vector<double> Sigma_f_;
    std::vector<double> dSigma_f_;
  };
  
}

#endif // DEPROJECTOR_H
  
