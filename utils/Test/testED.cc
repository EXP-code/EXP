#include <iostream>
#include <iomanip>
#include "ExpDeproj.H"

double projectedDensity(double R, double rmax, int nsteps, ExpDeproj & deproj)
{
  double dz = rmax / double(nsteps-1);
  double sum = 0.0;
  double curr = deproj.density(R), last = 0.0;
  for (int i=1; i<nsteps; i++) {
    last = curr;
    double z = dz * i;
    double r = sqrt(R*R + z*z);
    curr = deproj.density(r);
    sum += (curr + last) * dz;
  }
  return sum;
}

int main()
{
  ExpDeproj deproj;

  double rmin = 1.0e-4;
  double rmax = 30.0;
  int  nsteps = 4000;

  std::cout << "# Testing ExpDeproj class" << std::endl;
  std::cout << "# The density and mass columns are the deprojected model" << std::endl;
  std::cout << "# The projDensity column is the projected deprojected density calculated via numerical integration" << std::endl;
  std::cout << "# The error column is the difference between the projected density and exponential surface density" << std::endl;
  std::cout << "# " << std::endl;
  std::cout << std::setw(18) << "# R"
	    << std::setw(18) << "density"
	    << std::setw(18) << "mass"
	    << std::setw(18) << "projDensity"
	    << std::setw(18) << "error"
	    << std::endl;

  for (int i=0; i<100; i++) {
    double R = rmin + (0.5+i)*(rmax - rmin)/100.0;
    double dens = deproj.density(R);
    double mass = deproj.mass(R);
    double proj = projectedDensity(R, rmax, nsteps, deproj);
    double expt = 0.5 * exp(-R) / M_PI;
    std::cout << std::setw(18) << R
	      << std::setw(18) << dens
	      << std::setw(18) << mass
	      << std::setw(18) << proj
	      << std::setw(18) << (proj - expt)/expt
	      << std::endl;
  }
  
  return 0;
}
