#ifndef _king_H
#define _king_H

#include "massmodel.H"

class KingSphere : public SphericalModelTable
{
private:
  double kingK;
  double betaK;
  double E0K;

public:

  KingSphere(string filename) : SphericalModelTable(filename) {
    ModelID = "KingSphere"; }
      
  void setup_df(void) {
    if (get_num_param()<5) bomb("wrong number of parameters");
    kingK = get_param(2)*exp(-get_param(1));
    betaK = -2.0*get_param(3);
    E0K = -get_param(4)/get_param(5);
  }
  double distf(double E) { return kingK*(exp(betaK*(E-E0K)) - 1.0); }
  double dfde(double E)  { return betaK*kingK*exp(betaK*(E-E0K)); }
  double d2fde2(double E) { return betaK*betaK*kingK*exp(betaK*(E-E0K)); }
};

#endif
