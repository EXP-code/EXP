#ifndef _phase_H
#define _phase_H

#include <ostream>
#include <vector>

#include <Eigen/Eigen>

#include <numerical.H>

class Ensemble;

typedef Eigen::Vector3d (*force_func_ptr)(double, Eigen::Vector3d &, Eigen::Vector3d &);
typedef double (*pot_func_ptr)(double, Eigen::Vector3d &);

class Phase
{
protected:
  Eigen::Vector3d x;		// position vector
  Eigen::Vector3d v;		// velocity vector
  double m;			// mass
  double t;			// current time
  double work;			// total work done so far
  
  static double vscale;		    // typical length scale
  static double xscale;		    // typical speed scale
  static double tolerance;	    // desired integration accuracy
  static int potential_defined;	    // has a potential been set?
  static ode_integrator integrator; // type of integrator
  static symp_integrator mapping;   // type of sympletic integrator
  static char integrator_name[8];   // name of ODE integrator
  static char symp_integrator_name[8]; // name of SYMP integrator
public:
  // constructors
  
  Phase(void);
  Phase(const Phase &);
  
  
  // specification of model
  
  static void Register_Potential(pot_func_ptr p, force_func_ptr f)
  {potential_defined=1; Potential = p; Force = f;}
  static Eigen::Vector3d Default_Force(double, Eigen::Vector3d &, Eigen::Vector3d &);
  static double Default_Potential(double, Eigen::Vector3d &);
  static pot_func_ptr Potential;
  static force_func_ptr Force;
  
  
  
  // assignment operator
  
  Phase &operator=(const Phase &);
  
  
  
  // access to protected data
  
  Eigen::Vector3d &Position(void)	{return x;}
  Eigen::Vector3d &Velocity(void)	{return v;}
  double &Position(int i)	{return x[i];}
  double &Velocity(int i)	{return v[i];}
  double &Time(void)	{return t;}
  double &Work(void) 	{return work;}
  double &Mass(void)	{return m;}
  
  
  
  // set and get the scale dimensions and tolerance
  
  void set_Tolerance(double tol)	{tolerance = tol;}
  void set_Xscale(double xx) 	{xscale = xx;}
  void set_Vscale(double vv) 	{vscale = vv;}
  
  double get_Tolerance(void) 	{return tolerance;}
  double get_Xscale(void) 	{return xscale;}
  double get_Vscale(void) 	{return vscale;}
  
  
  // set or retreive the ode integrator type
  
  void set_Integrator(char *);
  char *get_Integrator(void);
  
  
  
  
  // print a phase-space point
  
  void print(std::ostream& out);
  
  
  // rotate frame or point about z axis
  
  Phase rotate_frame(double);
  Phase rotate_view(double);
  
  
  
  // various functions of the coordinates
  
  double Speed(void) 	{return sqrt(v.dot(v));}
  double Energy(void)	{return 0.5*m*v.dot(v) + m*(*Potential)(t, x);}
  double Jacobi(void);
  double Accuracy(Phase &);
  Eigen::Vector3d Angular_Momentum(void)	{return m*x.cross(v);}
  
  
  
  // numerical integration
  
  Phase integrate_to(double);
  Phase advance(double&);
  Phase x_return_map(void);
  
  // symplectic mapping
  
  Phase integrate_to_symplectic(double, double);
  Phase advance_symplectic(double &);
  
  // predictor-corrector stuff
  
  void set_PC(double, int);
  Phase PC_advance(double);
  
  
  void binread(std::istream& fp);
  void binwrite(std::ostream& fp);
  
  friend class Ensemble;
};


void Newton(double, Eigen::VectorXd&, Eigen::VectorXd&);
double v_circ(double);
double vertical_kappa(double);
double epicyclic_kappa(double);
double Oort_A(double);
double Oort_B(double);

class Frame_Rotation
{
public:
  double omega, omegasq, omega2;
  void set_omega(double w) 	
  {
    omega = w;
    omegasq = omega*omega;
    omega2 = 2.0*omega;
  }
  void set_corotation(double rc)
  {
    if (rc <= 0.0) omega=0.0;
    else omega = v_circ(rc)/rc;
    omegasq = omega*omega;
    omega2 = 2.0*omega;
  }
};


Eigen::Matrix3d z_Rotation(double);
typedef int (*freezing_function)(Phase &);



class Ensemble
{
protected:
  int Nstars;
  std::vector<Phase> stars;
  double t;
  
public:
  Ensemble(void);
  ~Ensemble(void);
  Ensemble(int);
  Ensemble(Ensemble &);
  
  void setsize(int);
  int get_Nstars(void) 	{return Nstars;}
  
  Ensemble &operator=(Ensemble &);
  
  Phase &operator[](int);
  
  Ensemble integrate_to(double);
  Ensemble Step_Positions(double);
  Ensemble Step_Velocities(double, Eigen::Vector3d *);
  Ensemble integrate_to(double, freezing_function);
  
  void write_log(std::ostream&);
  void write_snapshot(std::ostream&);
  void read_snapshot(std::istream&);
  void read_tidesfile(std::istream&);
  void write_orbits(std::ostream&);
  void read_orbits(std::istream&);
  
  void settime(double);
  void set_masses(double);
  void scale_masses(double);
  void scale_positions(double);
  void scale_speeds(double);
  void rotate_frame(double);
  void rotate_view(double);
  
  void translate(Eigen::Vector3d &);
  
  double total_Energy(void);
  double total_Kinetic(void);
  double total_Potential(void);
  double Virial(void);
  Eigen::Vector3d CM_Position(void);
  Eigen::Vector3d CM_Velocity(void);
  Eigen::Vector3d total_Angular_Momentum(void);
  
  Eigen::Matrix3d Moment_Tensor(double);
  Eigen::Matrix3d Inertia_Tensor(double);
  Eigen::Vector3d Principal_Axes(double, Eigen::Matrix3d &);
  Eigen::Vector3d Solid_Body_Frequency(double);
  
  void binread(std::istream& fp);
  void binwrite(std::ostream& fp);
  
};


double inverse_interp_zero(double, double, double, double,
			   double, double );

void forward_interp(double, double, double, double,
		    double, double, double, double&, double&);

#endif
