#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "Deprojector.H"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Deproject
{
  
  Deprojector::Deprojector(std::function<double(double)> SigmaFunc,
			   std::function<double(double)> dSigmaFunc,
			   double R_data_min,
			   double R_data_max,
			   double R_max_extend,
			   double tail_power,
			   int Ngrid)
    : sigma_func_(SigmaFunc), dsigma_func_(dSigmaFunc),
      Rdata_min_(R_data_min), Rdata_max_(R_data_max),
      tail_power_(tail_power)
  {
    if (Rdata_min_ <= 0.0 || Rdata_max_ <= Rdata_min_)
      throw std::runtime_error("Invalid R_data_min/R_data_max.");
    Rmin_ = Rdata_min_;
    Rmax_ = (R_max_extend > Rdata_max_) ? R_max_extend : Rdata_max_;
    build_grid(Ngrid);
  }
  
  Deprojector::Deprojector(const std::vector<double>& R_in,
			   const std::vector<double>& Sigma_in,
			   double R_max_extend,
			   double tail_power,
			   int Ngrid)
    : Rdata_min_(0.0), Rdata_max_(0.0), tail_power_(tail_power)
  {
    if (R_in.size() != Sigma_in.size() || R_in.size() < 2) throw std::runtime_error("Input R and Sigma must be same size and >=2.");
    // copy & sort
    std::vector<std::pair<double,double>> pairs;
    pairs.reserve(R_in.size());
    for (size_t i=0;i<R_in.size();++i) pairs.emplace_back(R_in[i], Sigma_in[i]);
    std::sort(pairs.begin(), pairs.end());
    std::vector<double> R(R_in.size()), S(R_in.size());
    for (size_t i=0;i<pairs.size();++i) { R[i] = pairs[i].first; S[i] = pairs[i].second; }
    
    // ensure positive radii and strictly increasing
    if (R.front() <= 0.0) {
      double eps = 1e-12;
      if (R.front() <= 0.0) R[0] = eps;
      for (size_t i=1;i<R.size();++i) if (R[i] <= R[i-1]) R[i] = R[i-1] + eps;
    }
    
    spline_.set_data(R, S);
    Rdata_min_ = spline_.xmin();
    Rdata_max_ = spline_.xmax();
    Rmin_ = Rdata_min_;
    Rmax_ = (R_max_extend > Rdata_max_) ? R_max_extend : Rdata_max_;
    
    sigma_func_ = [this](double rr){ return this->spline_.eval(rr); };
    dsigma_func_ = [this](double rr){ return this->spline_.deriv(rr); };
    
    build_grid(Ngrid);
  }
  
  
  // --- updated build_grid ---
  void Deprojector::build_grid(int Ngrid) {
    if (Ngrid < 3) Ngrid = 3;
    fineR_.resize(Ngrid);
    for (int i = 0; i < Ngrid; ++i) {
      double t = (double)i / (Ngrid - 1);
      fineR_[i] = Rmin_ + t * (Rmax_ - Rmin_);
    }
    
    // precompute spacing
    dx_.resize(Ngrid - 1);
    for (int i = 0; i < Ngrid - 1; ++i) dx_[i] = fineR_[i+1] - fineR_[i];
    
    Sigma_f_.assign(Ngrid, 0.0);
    dSigma_f_.assign(Ngrid, 0.0);
    
    bool have_dsf = static_cast<bool>(dsigma_func_);
    
    if (have_dsf) {
      for (int i = 0; i < Ngrid; ++i) {
	double rr = fineR_[i];
	if (rr <= Rdata_max_) {
	  Sigma_f_[i] = sigma_func_(rr);
	  dSigma_f_[i] = dsigma_func_(rr);
	} else {
	  double Sig_at = sigma_func_(Rdata_max_);
	  double factor = std::pow(rr / Rdata_max_, tail_power_);
	  Sigma_f_[i] = Sig_at * factor;
	  if (rr > 0.0)
	    dSigma_f_[i] = Sig_at * tail_power_ * std::pow(rr, tail_power_ - 1.0) / std::pow(Rdata_max_, tail_power_);
	  else
	    dSigma_f_[i] = 0.0;
	}
      }
    } else {
      // compute Sigma on grid, then finite-difference derivative using neighbor spacing
      for (int i = 0; i < Ngrid; ++i) {
	double rr = fineR_[i];
	if (rr <= Rdata_max_) Sigma_f_[i] = sigma_func_(rr);
	else {
	  double Sig_at = sigma_func_(Rdata_max_);
	  double factor = std::pow(rr / Rdata_max_, tail_power_);
	  Sigma_f_[i] = Sig_at * factor;
	}
      }
      
      for (int i = 0; i < Ngrid; ++i) {
	if (i > 0 && i < Ngrid - 1) {
	  // centered difference using grid neighbors (robust)
	  double x1 = fineR_[i-1], x2 = fineR_[i+1];
	  double y1 = Sigma_f_[i-1], y2 = Sigma_f_[i+1];
	  dSigma_f_[i] = (y2 - y1) / (x2 - x1);
	} else if (i == 0) {
	  // forward diff
	  double x1 = fineR_[1];
	  dSigma_f_[i] = (Sigma_f_[1] - Sigma_f_[0]) / (x1 - fineR_[0]);
	} else { // i == Ngrid-1
	  double x1 = fineR_[Ngrid-2];
	  dSigma_f_[i] = (Sigma_f_[Ngrid-1] - Sigma_f_[Ngrid-2]) / (fineR_[Ngrid-1] - x1);
	}
      }
    }
  }
  
  // --- updated rho_at ---
  double Deprojector::rho_at(double r) const {
    if (r >= Rmax_) return 0.0;
    
    // find index near r
    // choose local offset delta = 0.5 * local grid spacing to avoid singularity
    auto it0 = std::lower_bound(fineR_.begin(), fineR_.end(), r);
    int idx0 = (int)std::distance(fineR_.begin(), it0);
    // clamp idx0
    if (idx0 <= 0) idx0 = 1;
    if (idx0 >= (int)dx_.size()) idx0 = (int)dx_.size() - 1;
    
    double local_dx = dx_[std::max(0, idx0 - 1)];
    double delta = 0.5 * local_dx;                // half a local cell
    double rstart = r + delta;
    // ensure rstart >= fineR_[0]
    if (rstart < fineR_[0]) rstart = fineR_[0];
    
    // locate starting index after rstart
    auto it = std::lower_bound(fineR_.begin(), fineR_.end(), rstart);
    int idx = (int)std::distance(fineR_.begin(), it);
    if (idx >= (int)fineR_.size()) return 0.0;
    
    double integral = 0.0;
    int N = (int)fineR_.size();
    
    // integrate using trapezoid on intervals [R_{i-1}, R_i] for all i such that R_{i-1} >= rstart or partial first
    if (idx > 0) {
      // partial segment from a = max(fineR_[idx-1], rstart) to b = fineR_[idx]
      int i0 = idx - 1;
      double R0 = fineR_[i0], R1 = fineR_[i0+1];
      double a = std::max(R0, rstart);
      double b = R1;
      if (b > a) {
	// linear interpolation for derivative at 'a'
	double t = (a - R0) / (R1 - R0);
	double dSigma_a = dSigma_f_[i0] + t * (dSigma_f_[i0+1] - dSigma_f_[i0]);
	double f_a = - dSigma_a / std::sqrt(std::max(1e-300, a*a - r*r));
	double f_b = - dSigma_f_[i0+1] / std::sqrt(std::max(1e-300, b*b - r*r));
	integral += 0.5 * (f_a + f_b) * (b - a);
      }
      // full segments from i = idx+1 .. N-1
      for (int i = idx + 1; i < N; ++i) {
	double Rlo = fineR_[i-1], Rhi = fineR_[i];
	double f_lo = - dSigma_f_[i-1] / std::sqrt(std::max(1e-300, Rlo*Rlo - r*r));
	double f_hi = - dSigma_f_[i]   / std::sqrt(std::max(1e-300, Rhi*Rhi - r*r));
	integral += 0.5 * (f_lo + f_hi) * (Rhi - Rlo);
      }
    } else {
      // idx == 0 case: integrate over full segments starting at fineR_[0]
      for (int i = 1; i < N; ++i) {
	double Rlo = fineR_[i-1], Rhi = fineR_[i];
	double f_lo = - dSigma_f_[i-1] / std::sqrt(std::max(1e-300, Rlo*Rlo - r*r));
	double f_hi = - dSigma_f_[i]   / std::sqrt(std::max(1e-300, Rhi*Rhi - r*r));
	integral += 0.5 * (f_lo + f_hi) * (Rhi - Rlo);
      }
    }
    
    return integral / M_PI;
  }
  
  std::vector<double> Deprojector::rho(const std::vector<double>& r_eval) const {
    std::vector<double> out;
    out.reserve(r_eval.size());
    for (double r : r_eval) out.push_back(rho_at(r));
    return out;
  }
  
}
