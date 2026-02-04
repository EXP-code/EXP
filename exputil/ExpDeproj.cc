#include "ExpDeproj.H"

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

void ExpDeproj::initialize()
{
  if (ngrid < 2) {
    throw std::invalid_argument("n must be at least 2");
  }

  rv.resize(ngrid);
  mv.resize(ngrid);

  std::vector<double> dv(ngrid);

  double log_rmin = std::log(rmin);
  double log_rmax = std::log(rmax);
  double dlogr = (log_rmax - log_rmin)/static_cast<double>(ngrid - 1);

  for (int i = 0; i < ngrid; ++i) {
    rv[i] = std::exp(log_rmin + i*dlogr);
    dv[i] = 4.0*M_PI*rv[i]*rv[i]*density(rv[i]);
  }
  
  mv[0] = 0.0;
  for (int i=1; i<ngrid; i++)
    mv[i] = mv[i-1] + 0.5*(dv[i] + dv[i-1])*(rv[i] - rv[i-1]);
}

double ExpDeproj::density(double R)
{
  if (R < 0) {
    throw std::invalid_argument("R must be non-negative");
  }
  return 0.5*std::cyl_bessel_k(0.0, R)/(M_PI*M_PI);
}

double ExpDeproj::mass(double R)
{
  // Initialize precomputed arrays?
  if (mv.size() == 0) {
    initialize();
  }

  if (R < 0) {
    throw std::invalid_argument("R must be non-negative");
  }

  if (R < rmin) return 0.0;
  if (R > rmax) return mv.back();

  auto it = std::lower_bound(rv.begin(), rv.end(), R);
  // If R is slightly larger than the largest grid point due to rounding,
  // lower_bound may return rv.end(); in that case, use the last mass value.
  if (it == rv.end()) {
    return mv.back();
  }

  std::size_t idx = static_cast<std::size_t>(std::distance(rv.begin(), it));
  if (idx >= rv.size()) {
    // Defensive guard; should not happen after the it == rv.end() check.
    return mv.back();
  }
  if (rv[idx] == R) {
    return mv[idx];
  } else {
    // Ensure idx-1 is valid
    if (idx == 0) idx++;
    // Linear interpolation
    double x0 = rv[idx-1];
    double x1 = rv[idx  ];
    double y0 = mv[idx-1];
    double y1 = mv[idx  ];

    return y0 + (y1 - y0)*(R - x0)/(x1 - x0);
  }
}
