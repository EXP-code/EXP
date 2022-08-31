
#ifndef _TransformFFT_h
#define _TransformFFT_h

#include <Eigen/Dense>
#include <fftw3.h>

#include <complex>
#include <vector>

class TransformFFT
{
  int N, lo, hi;
  double dr, dk;

  fftw_plan p;
  std::vector<double> in;
  fftw_complex* out;

 public:
  
  TransformFFT(double dr, std::vector<double>& Y);
  TransformFFT(double dr, Eigen::VectorXd& Y);
  ~TransformFFT();

  void Power(Eigen::VectorXd& F, Eigen::VectorXd& P);
  void Power(std::vector<double>& F, std::vector<double>& P);

  void Inverse(Eigen::VectorXd& F, Eigen::VectorXcd& W);
  void Inverse(std::vector<double>& F, 
	       std::vector<double>& Wr, std::vector<double>& Wi);
};

#endif
