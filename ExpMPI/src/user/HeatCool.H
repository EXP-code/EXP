
class HeatCool
{
  typedef double(*func_1d)(double);
  
private:
  double n_h, g0, g1, g2, r;
  double crate, hrate, compcrate, trate;

  static double gp0_H, gp0_He, gp0_Hep;
  static double eps_H, eps_He, eps_Hep;

  double g_H(double t);
  double g_He(double t);
  double g_Hep(double t);
  double a_Hp(double t);
  double a_Hep(double t);
  double a_Hepp(double t);
  double a_p(double t);
  double rate_Hp(double t);
  double rate_Hep(double t);
  double rate_Hepp(double t);
  double rate_H(double t);
  double rate_H_cex(double t);
  double rate_He(double t);  
  double rate_br(double t);
  double radrate_Hp(double t);
  double radrate_Hep(double t);
  double radrate_Hepp(double t);

  void xion(double t, double *x, double *x_1, double *x_2, double *x_3, double *p_ne);
  static double fac1(double);
  static double fac2(double);
  static double romberg_o(func_1d func, double a, double b, double eps);
  static double spline_int(double* xa, double* ya, double* y2a, int n, double x);
  static void ionize();


public:


  HeatCool(double n, double T);

  static void initialize();

  double CoolRate() { return crate; }
  double HeatRate() { return hrate; }
  double CompCRate() { return compcrate; }
  double TotalRate() { return trate; }
};
