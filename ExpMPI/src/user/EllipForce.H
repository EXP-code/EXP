#ifndef EllipForce_H
#define EllipForce_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <interp.h>
#include <gaussQ.h>

class EllipForce {
private:

  double a, b, c, mass;
  int num, numr;
  LegeQuad *lq;

  vector<double> r, m, p;

public:

  EllipForce(double a, double b, double c, double mass, 
	     int num=100, int numr=100);
  ~EllipForce();

  double getMass(double x);
  double getPot(double x);
  void PrintTable();
  void TestTable();

};

#endif // EllipForce
