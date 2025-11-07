#ifndef _AddDisk_H
#define _AddDisk_H

#include <vector>
#include <memory>

#include "massmodel.H"

class AddDisk
{
 private:
  std::vector<double> r, d, m, p;
  std::shared_ptr<SphericalModelTable> mod;
  double dmass;

 public:
  static int number;			// 4000
  static double Rmin;			// 1.0e-3
  static bool use_mpi;			// false
  static bool logarithmic;		// false

  AddDisk(std::shared_ptr<AxiSymModel> halo, std::shared_ptr<AxiSymModel> disk, double dmass);
  ~AddDisk();

  std::shared_ptr<SphericalModelTable> get_model() { return mod; }
};


#endif
