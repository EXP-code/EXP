#ifndef SIM_ANN_HPP_
#define SIM_ANN_HPP_ 1.4

#define _USE_MATH_DEFINES

#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <limits>
#include <cmath>

class SimAnn
{
  typedef std::function<double(std::vector<double>&)> Func1d;
  typedef std::function<std::vector<double>(std::vector<double>&)> FuncMap;

 private:

  Func1d func;			// Objective function
  FuncMap pmap;			// Parameter mapping for output

  int ddwell;
  double rrange;
  double t0;
  double K;
  double rho;
  double dt;		  // temperature increment to use when melting
  double tscale	;
  int maxit;
  double c_jump;		// phase transition step size
  int fsave;

  std::mt19937 gen;	       // Random number generator
  std::uniform_real_distribution<double> number_range, number_01;

  int err;
  std::vector<double> x, xnew, xbest;
  double y, dy, ybest;

  int equilibrate(const double t, const int n);
  std::string fname;
  void log_state(int);
  std::ostream *str;

 public:

  SimAnn() : ddwell(20), rrange(M_PI_2),
	     t0(0.0), K(1.0), rho(0.5), dt(0.1), tscale(0.1),
	     maxit(400), c_jump(100.0), fsave(0)
  {
    number_range = std::uniform_real_distribution<double>(-rrange, rrange);
    number_01    = std::uniform_real_distribution<double>(0.0, 1.0);
  }

  SimAnn(std::function<double(std::vector<double>&)> f, const int d = 1);

  ~SimAnn() 
  { 
    // Nothing so far
  }
  
  //! Set the objective function
  void set_up(Func1d f, const int seed=10);
  
  //! Set map function
  void set_map(FuncMap f) { pmap = f; }
  
  //! Look for the melting temprature
  double melt  (const int iters = -1);

  //! Perform the annealing
  double anneal(const int iters = -1);
  
  int iterations(const int m = -1) { if ( m > 0 ) maxit = m;
				     return maxit; }

  int dwell(const int d = -1)      { if ( d > 0 ) ddwell = d;
				     return ddwell; }

  double Boltzmann(const double k = -1.0)
  { if ( k > 0.0 ) K = k; return K; }
  
  double learning_rate(const double r = -1.0)
  { if ( r > 0.0 ) rho = r; return rho; }

  double temperature(const double t = -1.0)
  { if ( t > 0.0 ) t0 = t; return t0; }

  double jump(const double j = -1.0)
  { if ( j > 0.0 ) c_jump = j; return c_jump; }

  double range(const double r = -1.0)
  {
    if ( r > 0.0 ) {
      rrange = r;
      number_range = std::uniform_real_distribution<double>(-r,r);
    }
    return rrange; 
  }

  void initial (std::vector<double>& xi) { x = xi; }
  std::vector<double> current() { return x;     }
  std::vector<double> optimum() { return xbest; }

  void save_states(const std::string& name)
  {
    fsave = 1;
    fname = name;
  }

  void report_states(std::ostream& out)
  {
    fsave = 2;
    str   = &out;
  }

};

#endif

