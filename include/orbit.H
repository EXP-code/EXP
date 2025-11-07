#ifndef _orbit_H
#define _orbit_H

#include <functional>
#include <memory>
#include <Eigen/Eigen>
#include "gaussQ.H"

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
  double ap, am, sp, sm;

  std::shared_ptr<AxiSymModel>  model;
  std::shared_ptr<AxiSymBiorth> biorth;
  int RECUR;

  void compute_action(void) { compute_freq(); }
  void compute_freq(void);
  void compute_angles(void);
  void compute_angles_old(void);
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

  std::tuple<double, double, bool> search
  (std::function<double(double)> func, double rmin, double rmax);
  
public:

  static int    Nseg;           // Number of trial segments for root finding 

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

  //@{
  //! Required functions for accessing phase-space parameters

  /** Return the action where:
      index 0 is radial action, I_r
      index 1 is angular action, J
      index 2 should be J_z/J, current set to zero
  */
  double get_action(const int i)
  {
    if (i>dof)           bomb("actions not defined!");
    if (i<0 or i>2)      bomb("action index must be [0,1,2]");
    if (!action_defined) compute_action();
  
  return action[i];
  }

  //! Get radial angle for a given radius
  double get_w1(double r, double vr);

  /** Return the orbital frequency where:
      index 0 is radial frequency, I_r
      index 1 is angular frequency, J
      index 2 is frequency of line of nodes (zero for a spherical system)
  */
  double get_freq(const int i)
  {
    if (i>dof)         bomb("frequencies not defined!");
    if (i<0 or i>2)    bomb("frequency index must be [0,1,2]");
    if (!freq_defined) compute_freq();
    
    return freq[i];
  }

  /** Get angle dependent quantities for the orbit where the
     index signifies:

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
  //@}

  //! Set the energy and kappa for the orbit
  void new_orbit(double Energy, double kappa, double beta = 0.0);

  //! Set the biorthogonal ffunction type
  void set_biorth(std::shared_ptr<AxiSymBiorth> type,
		  int ll, int max, int recur = 0) 
    { nmax = max; l = ll; biorth = type; RECUR = recur; 
    biorth_defined = false; }

  //! Set the numerical parameters for angle transform computation
  void set_numerical_params(int RECS=64, int NMAX=40, int NBSCT=3)
    { recs=RECS; nmax=NMAX; nbsct=NBSCT;}

  double pot_trans(int l1, int l2, std::function<double(double)> func );
  double pot_trans(int l1, int l2, int n );
  void   pot_trans(int l1, int l2, Eigen::VectorXd& t);

  //! Access to underlying grid for pot_trans
  struct ANGLE_GRID * get_angle_grid(void)
  { 
    if (!freq_defined) compute_freq();
    if (!angle_defined) compute_angles();
    return &angle_grid;
  }
  
  //@{
  //! Safe access to internal quantities

  //! The value of the orbital energy
  double Energy(void) { return energy; }

  //! The value of the orbital angular momentum
  double AngMom(void)
  { 
    if (freq_defined == false) compute_freq();
    return jmax*kappa;
  }

  //! The value of the current J_{max}(E)
  double Jmax(void)
  {
    if (freq_defined == false) compute_freq();
    return jmax;
  }

  //! The value of the current J/J_{max}(E)
  double Kappa(void) { return kappa; }

  //! The specified altitude of the orbital plane
  double Beta(void)  { return beta; }
  
  //! Radius of the pericenter of this orbit
  double peri(void)
  { 
    if (freq_defined == false) compute_freq();
    return r_peri;
  }

  //! Radius of the apocenter of this orbit
  double apo(void)
  {
    if (freq_defined == false) compute_freq();
    return r_apo;
  }

  //! Radius of the circular orbit with energy E
  double circ(void)
  {
    if (freq_defined == false) compute_freq();
    return r_circ;
  }

  //! Get the underlying axisymmetric model object
  AxiSymModel& modl(void) { return *model; }

  //! Get the underlying biorthogonal function object
  AxiSymBiorth& orth(void) { return *biorth; }

  //@}

};

#endif
