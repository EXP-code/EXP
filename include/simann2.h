// A general purpose Simulated Annealing Class
// This version allows vector data

#ifndef SIM_ANNEAL_HPP_
#define SIM_ANNEAL_HPP_ 1.3

#include <functional>
#include <string>
#include <random>

class SimAnneal
{
 private:

  typedef std::function<double(std::vector<double>&)> Func1d;

  Func1d func;                  // 
  int dimension;		// 
  int ddwell;			// 
  double rrange;		// 
  double t0;			// 
  double K;			// 
  double rho;			// 
  double dt;			// temperature increment to use when melting
  double tscale	;		// 
  int maxit;			// 
  double c_jump;		// phase transition step size
  int fsave;			// 

  int err;
  std::vector<double> x, xnew, xbest;
  double y, dy, ybest;

				// Random number distributions
  std::uniform_real_distribution<> number_range, number_01;

  int equilibrate(const double t, const int n);
  std::string fname;
  void log_state(int);

 public:

  SimAnneal() :	 func(NULL), dimension(1), ddwell(20), rrange(M_PI/2.0), 
    t0(0.0), K(1.0), rho(0.5), dt(0.1), tscale(0.1), maxit(400), c_jump(100.0),
    fsave(0) {

    std::uniform_real_distribution<>::param_type params1(-rrange, rrange);
    std::uniform_real_distribution<>::param_type params2(0.0, 1.0);

    number_range.param(params1);
    number_01.param(params2);
  }

  SimAnneal(Func1d f, const int d = 1);

  ~SimAnneal() 
  { 
    // NONE
  }
  
  int set_up(Func1d f, const int d = 1, const uint32_t seed=10);
  
  const int operator!() const { return err; }
  
  double melt(const int iters = -1);
  double anneal(const int iters = -1);
  
  int iterations(const int m = -1) { if ( m > 0 ) maxit = m;
				     return maxit; }
  int dwell(const int d = -1)      { if ( d > 0 ) ddwell = d;
				     return ddwell; }
  double Boltzmann(const double k = -1.0)
    { if ( k > 0.0 ) K = k;
      return K; }
  
  double learning_rate(const double r = -1.0)
    { if ( r > 0.0 ) rho = r;
      return rho; }
  double temperature(const double t = -1.0)
    { if ( t > 0.0 ) t0 = t;
      return t0; }
  double jump(const double j = -1.0)
    { if ( j > 0.0 ) c_jump = j;
      return c_jump; }
  double range(const double r = -1.0)
    { 
      if ( r > 0.0 ) 
	{ rrange = r;
	  std::uniform_real_distribution<>::param_type params1(-rrange, rrange);
	  number_range.param(params1);}
      return rrange; 
    }
  void initial(std::vector<double>& xinit) { x = xinit;    }
  void current(std::vector<double>& xcur)  { xcur = x;     }
  void optimum(std::vector<double>& xopt)  { xopt = xbest; }
  void save_states(const char *name) {
    fsave = 1;
    fname = std::string(name);
  }
};

#endif

