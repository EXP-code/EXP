#include <cmath>
#include <iostream>
#include <iomanip>

#include <TransformFFT.H>

TransformFFT::TransformFFT(double DR, std::vector<double>& Y)
{
  N = Y.size();

  dr = DR;
  dk = 2.0*M_PI/dr/N;

  in = Y;
  out.resize(N/2+1);
  
  p = fftw_plan_dft_r2c_1d(N, in.data(),
			   reinterpret_cast<fftw_complex*>(out.data()),
			   FFTW_ESTIMATE);

  fftw_execute(p);
}


TransformFFT::TransformFFT(double DR, Eigen::VectorXd& Y)
{
  N = Y.size();

  dr = DR;
  dk = 2.0*M_PI/dr/N;

  in.resize(N);
  for (int j=0; j<N; j++) in[j] = Y[j];

  out.resize(N/2+1);
  
  p = fftw_plan_dft_r2c_1d(N, in.data(),
			   reinterpret_cast<fftw_complex*>(out.data()),
			   FFTW_ESTIMATE);

  
  fftw_execute(p);
}


TransformFFT::~TransformFFT()
{
  fftw_destroy_plan(p);
}

void TransformFFT::Power(Eigen::VectorXd& F, Eigen::VectorXd& P)
{
  F.resize(N/2+1);
  P.resize(N/2+1);
  
  double d2 = dr*dr;

  F(0) = 0.0;
  P(0) = d2 * std::norm(out[0]);

  for (int j=1; j<N/2; j++) {
    F(j) = dk * j;
    P(j) = d2 * std::norm(out[j]) * 2.0;
  }
  
  F(N/2) = dk * N/2;
  P(N/2) = d2 * std::norm(out[N/2]);
}

void TransformFFT::Power(std::vector<double>& F, std::vector<double>& P)
{
  F.resize(N/2+1);
  P.resize(N/2+1);

  double d2 = dr*dr;

  F[0] = 0.0;
  P[0] = d2 * std::norm(out[0]);

  for (int j=1; j<N/2; j++) {
    F[j] = dk * j;
    P[j] = d2 * std::norm(out[j]) * 2.0;
  }
  
  F[N/2] = dk * N/2;
  P[N/2] = d2 * std::norm(out[N/2]);
}


void TransformFFT::Inverse(Eigen::VectorXd& F, Eigen::VectorXcd& W)
{
  F.resize(N/2+1);
  W.resize(N/2+1);

  for (int j=0; j<N/2+1; j++) {
    F[j] = dk*j;
    W[j] = out[j] * dr;
  }
  
}

void TransformFFT::Inverse(std::vector<double>& F, 
			   std::vector<double>& Wr, std::vector<double>& Wi)
{
  F .resize(N/2+1);
  Wr.resize(N/2+1);
  Wi.resize(N/2+1);

  for (int j=0; j<N/2+1; j++) {
    F[j] = dk*j;
    Wr[j] = std::real(out[j]) * dr;
    Wi[j] = std::imag(out[j]) * dr;
  }
}

