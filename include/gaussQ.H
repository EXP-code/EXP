/*

SYNOPSIS
     #include <gaussQ.H>


DESCRIPTION

     The functions compute the abscissas and  weight  factors  of
     various   Gaussian   quadrature   formulas.

     For the "Jacobi" series, the weight function is
     w(x) = x**alpha * (1-x)**beta
     over the interval 0 <= x <= 1, with alpha, beta > -1.

     For the "Legendre" series, the weight function is
     w(x) = 1
     over the interval 0 <= x <= 1
     (This is computed as a special case of the Jacobi series)

     For the "Laguerre" series, the weight function is
     w(x) = x**alpha * exp(-x)
     over the interval 0 <= x < +infinity, with alpha > -1.

     For the "Hermite" series, the weight function is
     w(x) = abs(x)**alpha * exp(-x**2)
     over the interval -infinity < x < +infinity,
     with alpha > -1.

DIAGNOSTICS
     Aborts with error message on creation if alpha or beta <= -1.

HISTORY
     22-Sep-85  asselin at Carnegie-Mellon University
          Created. (C-code)

     01/18/94 C++ wrapper MDW

*/

#ifndef _gaussQ_H
#define _gaussQ_H

using namespace std;

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>

class GaussQuad 
{
public:

  std::string FunctionID;

  std::vector<double> w;
  std::vector<double> r;
  double alpha;
  double beta;
  int n;

  GaussQuad() : alpha(0), beta(0), n(0) {}

  double weight(const int i) { if (i<0 || i>=n) return bomb("index out of bounds");
  return w[i]; }

  double knot(const int i) { if (i<0 || i>=n) return bomb("index out of bounds");
  return r[i]; }

  std::vector<double>& wV(void) { return w; }
  std::vector<double>& kV(void) { return r; }
  int get_n(void) { return n; }
  double get_alpha(void) { return alpha; }
  double get_beta(void) { return beta; }

  int bomb(const char *s) {
    std::ostringstream sout;
    sout << "ERROR from " << FunctionID << ": " << s;
    throw std::runtime_error(sout.str());
  }
};

class HermQuad : public GaussQuad
{
public:
  HermQuad(int N=10, double ALPHA=0.0);
};

class LaguQuad : public GaussQuad
{
public:
  LaguQuad(int N=10, double ALPHA=0.0);
};

class JacoQuad : public GaussQuad
{
public:
  JacoQuad(int N=10, double ALPHA=0.0, double BETA=0.0);
};

class LegeQuad : public JacoQuad
{
public:
  LegeQuad(int N=10) : JacoQuad(N, 0.0, 0.0) { FunctionID = "LegeQuad"; };
};

#endif
