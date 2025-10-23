#ifndef _MakeModel_H
#define _MakeModel_H

#include <string>
#include <memory>

#include <Eigen/Eigen>

#include "massmodel.H"

class MakeModel 
{
private:

  int num;
  double rmin, rmax, dr;
  bool logr;
  Eigen::VectorXd r, d, m, p, p1; 
  std::vector<Eigen::VectorXd> dm;


public:
  
  MakeModel(int N, double Rmin, double Rmax, bool logr=true);

  void AddPoint(double r, double mass);

  std::shared_ptr<SphericalModelTable> Compute();

  void WriteModel(string& filename);
};

#endif
