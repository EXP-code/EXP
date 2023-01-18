#include <functional>
#include <vector>

class HeatCool
{
  struct Rates 
  {
    double crate;
    double hrate;
    double cmpcrate;
    double trate;
  };
  
private:
  double n_h, g0, g1, g2, r;
  double crate, hrate, cmpcrate, trate;

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
  static double romberg_o(std::function<double(double)> func,
			  double a, double b, double eps);
  static double spline_int(double* xa, double* ya, double* y2a, int n, double x);
  static void ionize();


  double Nmin, Nmax, Tmin, Tmax, dN, dT;
  unsigned Nnum, Tnum;
  string Cache;
  bool table;

  vector< vector<Rates> > data;

  void compute(double n, double T);

  bool readCache();
  void writeCache();
  void debugCache();

public:

  static bool dbgCache;
  static bool useCache;

  HeatCool() : table(false) {}

  HeatCool(double n, double T) : table(false) { compute(n, T); }

  HeatCool(double nmin, double nmax, double Tmin, double Tmax,
	   unsigned nnum, unsigned Tnum, string cache="");

  static void initialize();

  void setPoint(double n, double T);

  double CoolRate() { return crate; }
  double HeatRate() { return hrate; }
  double Compcrate() { return cmpcrate; }
  double TotalRate() { return trate; }

  double CoolRate(double N, double T);
};
