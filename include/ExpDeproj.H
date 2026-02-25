#pragma once

#include <vector>
#include <cmath>

class ExpDeproj
{
  const double rmin = 1.0e-4;
  const double rmax = 30.0;

  // Precomputed radial grid
  int ngrid;
  std::vector<double> rv, mv;
  void initialize();

public:
  //! Constructor
  ExpDeproj(int n=4000) : ngrid(n) { initialize(); }

  //! Destructor
  virtual ~ExpDeproj() {}

  //! Density profile at radius R
  double density(double R);

  //! Enclosed mass at radius R
  double mass(double R);
};
