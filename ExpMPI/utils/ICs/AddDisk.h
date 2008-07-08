#ifndef _AddDisk_H
#define _AddDisk_H

#include <vector>

#include <massmodel.h>

class AddDisk
{
 private:
  vector<double> r, d, m, p;
  SphericalModelTable *mod;
  double dmass;

 public:
  static int number;			// 4000
  static double Rmin;			// 1.0e-3
  static bool use_mpi;			// false
  static bool logarithmic;		// false

  AddDisk(AxiSymModel* halo, AxiSymModel* disk, double dmass);
  ~AddDisk();

  SphericalModelTable *get_model() { return mod; }
};


#endif
