#ifndef _AddSpheres_H
#define _AddSpheres_H

#include <vector>
#include <memory>

#include <massmodel.H>

class AddSpheres
{
 private:
  std::vector<double> r, d, m, p;
  std::shared_ptr<SphericalModelTable> mod;
  double bmass;

 public:
  static int number;			// 4000
  static double Rmin;			// 1.0e-3
  static bool use_mpi;			// false
  static bool logarithmic;		// false

  AddSpheres(std::shared_ptr<AxiSymModel> halo,
	     std::shared_ptr<AxiSymModel> bulge, double bmass);
  ~AddSpheres();

  std::shared_ptr<SphericalModelTable> get_model() { return mod; }
};


#endif
