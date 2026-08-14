#ifndef _orbit1d_H
#define _orbit1d_H

#include <functional>
#include <complex>
#include <memory>

#include <Eigen/Eigen>

#include "orbit.H"
// #include "massmodel1d.H"

class OneDModel;
class OneDBiorth;

class OneDOrbit : public RegularOrbit
{

  using Vector  = Eigen::VectorXd;
  using Matrix  = Eigen::MatrixXd;
  using CVector = Eigen::VectorXcd;
  using CMatrix = Eigen::MatrixXcd;
  using Complex = std::complex<double>;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  double energy;
  double vp;			// keep track orbit information in in-plane
				// direction

				// computed in compute_freq
  double zmax;
  double freq;
  double action;

  std::shared_ptr<OneDModel>  model;
  std::shared_ptr<OneDBiorth> biorth;

  void compute_action(void) { compute_freq(); }
  void compute_freq(void);
  void compute_angles(void);
  void compute_biorth(void);

  struct OneD_Angle_Grid {
    int size;
    Vector time;
    Matrix angl;
    Matrix hght;
  } angle_grid;

  constexpr static double TOL = 1.0e-10;

				// Number of points in angle grid
  int NGRID;
				// Number of knots for angle integrals
  int NINT;
				// Number of points for frequency integral
  int NFREQ;
				// Number of points for potential transform
  int NREC;

  std::shared_ptr<LegeQuad> lkw; // Cache this for computing angles

  bool pot_cache = false;	 // Cache potential transforms for a
				 // given biorthogonal set

  std::vector<Vector> pot_eval;	 // Cache potential evaluations for pot_trans


public:

  enum class Type { Angle, Height, Velocity };

  //@{
  //! Constructors
  OneDOrbit(void);
  OneDOrbit(std::shared_ptr<OneDModel> model);
  OneDOrbit(std::shared_ptr<OneDModel> model, double Energy, double vparallel=0.0);
  //@}

  virtual ~OneDOrbit() = default;

  //! Copy constructor and assignment operator (defaulted)
  OneDOrbit(OneDOrbit &orb);
  //  OneDOrbit &operator=(OneDOrbit &);

				// Required functions

  double get_action(const int i=0) {
    if (action_defined == false) compute_action();
    return action; }

  double get_freq(const int i=0) {
    if (freq_defined == false) compute_freq();
    return freq; }

  double get_angle(const int i, const double time);

  double get_angle(const Type f, const double time);

				// Specific functions

  void new_orbit(double Energy, double vparallel=0.0);

  void set_biorth(std::shared_ptr<OneDBiorth> type) {
    biorth = type;
    biorth_defined = true;
  }

  void set_numerical_params(int nrec, int ngrid=100, int nfreq=64, 
			    int nint=20) 
  { NREC = nrec; NGRID = ngrid; NFREQ = nfreq; NINT = nint;}

  double pot_trans(const int k, double (*func)(double) );
  double pot_trans(const int k, const int n);
  void   pot_trans(const int k, Eigen::VectorXd& t);
  
  Complex pot_trans_complex(const int k, std::function<double(double)> );
  Complex pot_trans_complex(const int k, const int n);
  Complex pot_trans_complex_check(const int k, const int n);
  Complex pot_trans_complex_check(const int k, const int nx, const int n);
  void    pot_trans_complex(const int k, Eigen::VectorXcd& t);
  void    set_pot_cache(bool val, int n=0);
  
  double pot_trans(const int kx, const int ky, double (*func)(double) );
  double pot_trans(const int kx, const int ky, const int n);
  void   pot_trans(const int kx, const int ky, Eigen::VectorXd& t);
  
  Complex pot_trans_complex(const int kx, const int ky, const int kz,
			    double (*func)(double) );
  Complex pot_trans_complex(const int kx, const int ky, const int kz,
			    const int n);
  void    pot_trans_complex(const int kx, const int ky, const int kz,
			    Eigen::VectorXcd& t);
  void    pot_trans_complex_check(const int kx, const int ky, const int kz,
				  Eigen::VectorXcd& t);
  
  //@{
  // Safe access
  const double Energy(void) { return energy; }
  const double Vpara(void) { return vp; }

  const double apo(void) {
    if (freq_defined == false) compute_freq();
    return zmax;
  }

  OneDModel& modl(void) { return *model; }

  OneDBiorth& orth(void) { return *biorth; }
  //@}
};

#endif
