#ifndef EllipForce_H
#define EllipForce_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>

#include "interp.H"
#include "gaussQ.H"

class EllipForce {
private:

  double a, b, c, mass;
  int num, numr;
  std::shared_ptr<LegeQuad> lq;

  std::vector<double> r, m, p;

  void write_cache();
  bool read_cache();

public:

  static string cache_file;

  EllipForce(double a, double b, double c, double mass, 
	     int num=100, int numr=100);
  ~EllipForce();

  double getMass(double x);
  double getPot(double x);
  void PrintTable();
  void PrintModel(const string&);
  void TestTable();

};

#endif // EllipForce
