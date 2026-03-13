#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "Deprojector.H"

using namespace Deproject;

int main()
{
  // Example A: construct from sampled data
  {
    std::vector<double> Rdata, Sigma;
    int Ndata = 2000;
    double Rmin = 0.01, Rmax = 10.0;
    for (int i = 0; i < Ndata; ++i) {
      double t = (double)i / (Ndata - 1);
      double r = Rmin + t * (Rmax - Rmin);
      Rdata.push_back(r);
      Sigma.push_back(std::pow(1.0 + r*r, -1.5));
    }

    Deprojector D(Rdata, Sigma, /*R_max_extend=*/50.0, /*tail_power=*/-4.0, /*Ngrid=*/6000);

    std::vector<double> r_eval;
    int Nr = 150;
    for (int i = 0; i < Nr; ++i) {
      double t = (double)i / (Nr - 1);
      r_eval.push_back(0.01 + t * 8.0);
    }
    auto rho = D.rho(r_eval);
    
    std::ofstream ofs("rho_from_sampled.txt");
    for (size_t i = 0; i < r_eval.size(); ++i) ofs << r_eval[i] << " " << rho[i] << "\n";
    ofs.close();
    std::cout << "Wrote rho_from_sampled.txt\n";
  }
  
  // Example B: construct from analytic functor + analytic derivative
  {
    auto SigmaFunc = [](double R)->double { return std::pow(1.0 + R*R, -1.5); };
    auto dSigmaFunc = [](double R)->double { return -3.0 * R / std::pow(1.0 + R*R, -2.5); }; // analytic derivative
    
    Deprojector D(SigmaFunc, dSigmaFunc, /*R_data_min=*/0.01, /*R_data_max=*/10.0,
		  /*R_max_extend=*/50.0, /*tail_power=*/-4.0, /*Ngrid=*/6000);

    std::vector<double> r_eval;
    int Nr = 150;
    for (int i = 0; i < Nr; ++i) {
      double t = (double)i / (Nr - 1);
      r_eval.push_back(0.01 + t * 8.0);
    }
    auto rho = D.rho(r_eval);

    std::ofstream ofs("rho_from_functor.txt");
    for (size_t i = 0; i < r_eval.size(); ++i) ofs << r_eval[i] << " " << rho[i] << "\n";
    ofs.close();
    std::cout << "Wrote rho_from_functor.txt\n";
  }
  
  return 0;
}
