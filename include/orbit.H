#ifndef _orbit_H
#define _orbit_H

#include <memory>
#include <Eigen/Eigen>
#include <gaussQ.H>

class AxiSymModel;
class AxiSymBiorth;

class RegularOrbit
{
public:

  int dof;
  std::string OrbitID;

  bool model_defined;
  bool action_defined;
  bool freq_defined;
  bool angle_defined;
  bool biorth_defined;

  Eigen::Vector3d action, freq;
  
  virtual ~RegularOrbit() {}

  virtual double get_action(const int) = 0;
  virtual double get_freq(const int) = 0;
  virtual double get_angle(const int, const double) = 0;

				// Error function
  void bomb(const char*   s);
  void bomb(const string& s);
  void warn(const string& fct, const string& msg);
};
  
struct ANGLE_GRID {
  Eigen::MatrixXd t;  // RECS
  Eigen::MatrixXd w1;
  Eigen::MatrixXd dw1dt;
  Eigen::MatrixXd f;
  Eigen::MatrixXd r;
  Eigen::MatrixXd fr; // NMAX x RECS
  int num;
};



class SphericalOrbit : public RegularOrbit
{

private:
  double energy;
  double kappa;
  double beta;

  double jmax;
  double r_peri;
  double r_apo;
  double r_circ;

  std::shared_ptr<AxiSymModel>  model;
  std::shared_ptr<AxiSymBiorth> biorth;
  int RECUR;

  void compute_action(void) { compute_freq(); }
  void compute_freq(void);
  void compute_angles(void);
  void compute_freq_epi(void);
  void compute_angles_epi(void);
  void compute_biorth(void);

  struct ANGLE_GRID angle_grid;
  int recs;
  int nmax;
  int nbsct;

  int l;
  int l1s, l2s;
  Eigen::VectorXd cosvec;

  GaussQuad Gkn;

public:

  static double RMAXF;		// Factors of model rmax for upper 
				// turning point search [default=3.0]

  static double tol;		// Root finder tolerance [default=1.0e-8]

  static double ZFRAC;		// Factor for minimum radius
				// turning point search [default=1.0e-3]

  static double tolnr;		// r_apo location using Newton-Raphson 
				// refinement				

				// Error if Energy < Phi(Rmin) [default=true]
  static bool guard;


  // Constructors


  SphericalOrbit(void);
  SphericalOrbit(std::shared_ptr<AxiSymModel> model);
  SphericalOrbit(std::shared_ptr<AxiSymModel> model,
		 double Energy, double kappa, double beta = 0.0);

  virtual ~SphericalOrbit() {}

  SphericalOrbit(const SphericalOrbit &);
  SphericalOrbit &operator=(const SphericalOrbit &);

				// Required functions

  double get_action(const int i) {
    if (i>dof) bomb("actions not defined!");
    if (!action_defined) compute_action();
    return action[i]; }

  double get_freq(const int i) {
    if (i>dof) bomb("frequencies not defined!");
    if (!freq_defined) compute_freq();
    return freq[i]; }

  /*
     i     val
     -     ---
     1      w1
     2      w2
     3      t
     4      dw1dt
     5      f
     6      r
     7      phi
   */

  double get_angle(const int i, const double time);

				// Specific functions

  void new_orbit(double Energy, double kappa, double beta = 0.0);

  void set_biorth(std::shared_ptr<AxiSymBiorth> type,
		  int ll, int max, int recur = 0) 
    { nmax = max; l = ll; biorth = type; RECUR = recur; 
    biorth_defined = false; }

  void set_numerical_params(int RECS=64, int NMAX=40, int NBSCT=3)
    { recs=RECS; nmax=NMAX; nbsct=NBSCT;}

  double pot_trans(int l1, int l2, double (*func)(double) );
  double pot_trans(int l1, int l2, int n );
  void   pot_trans(int l1, int l2, Eigen::VectorXd& t);

				// Access to underlying grid for pot_trans
  struct ANGLE_GRID * get_angle_grid(void) { 
    if (!freq_defined) compute_freq();
    if (!angle_defined) compute_angles();
    return &angle_grid; }
  
				// Safe access

  double Energy(void) { return energy; }

  double AngMom(void)
  { 
    if (freq_defined == false) compute_freq();
    return jmax*kappa;
  }

  double Jmax(void)
  {
    if (freq_defined == false) compute_freq();
    return jmax;
  }

  double Kappa(void) { return kappa; }

  double Beta(void)  { return beta; }
  
  double peri(void)
  { 
    if (freq_defined == false) compute_freq();
    return r_peri;
  }

  double apo(void)
  {
    if (freq_defined == false) compute_freq();
    return r_apo;
  }

  double circ(void)
  {
    if (freq_defined == false) compute_freq();
    return r_circ;
  }

  AxiSymModel& modl(void) { return *model; }

  AxiSymBiorth& orth(void) { return *biorth; }

};

#endif
