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
  out = new fftw_complex [N/2+1];
  
  p = fftw_plan_dft_r2c_1d(N, in.data(), out, FFTW_ESTIMATE);

  fftw_execute(p);
}


TransformFFT::TransformFFT(double DR, Eigen::VectorXd& Y)
{
  N = Y.size();

  dr = DR;
  dk = 2.0*M_PI/dr/N;

  in.resize(N);
  for (int j=0; j<N; j++) in[j] = Y[j];

  out = new fftw_complex [N/2+1];
  
  p = fftw_plan_dft_r2c_1d(N, in.data(), out, FFTW_ESTIMATE);

  
  fftw_execute(p);
}


TransformFFT::~TransformFFT()
{
  delete [] out;
  fftw_destroy_plan(p);
}

void TransformFFT::Power(Eigen::VectorXd& F, Eigen::VectorXd& P)
{
  F.resize(N/2+1);
  P.resize(N/2+1);
  
  double d2 = dr*dr;

  F(0) = 0.0;
  P(0) = d2 * (out[0][0]*out[0][0] + out[0][1]*out[0][1]);

  for (int j=1; j<N/2; j++) {
    F(j) = dk * j;
    P(j) = d2 * (out[j][0]*out[j][0] + out[j][1]*out[j][1]) * 2.0;
  }
  
  F(N/2) = dk * N/2;
  P(N/2) = d2 * (out[N/2][0]*out[N/2][0] + out[N/2][1]*out[N/2][1]);
}

void TransformFFT::Power(std::vector<double>& F, std::vector<double>& P)
{
  F.resize(N/2+1);
  P.resize(N/2+1);

  double d2 = dr*dr;

  F[0] = 0.0;
  P[0] = d2 * ( out[0][0]*out[0][0] + out[0][1]*out[0][1] );

  for (int j=1; j<N/2; j++) {
    F[j] = dk * j;
    P[j] = d2 * out[j][0]*out[j][0] + out[j][1]*out[j][1] * 2.0;
  }
  
  F[N/2] = dk * N/2;
  P[N/2] = d2 * (out[N/2][0]*out[N/2][0] + out[N/2][1]*out[N/2][1]);
}


void TransformFFT::Inverse(Eigen::VectorXd& F, Eigen::VectorXcd& W)
{
  F.resize(N/2+1);
  W.resize(N/2+1);

  for (int j=0; j<N/2+1; j++) {
    F[j] = dk*j;
    W[j] = std::complex<double>(out[j][0], out[j][1]) * dr;
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
    Wr[j] = out[j][0] * dr;
    Wi[j] = out[j][1] * dr;
  }
}

