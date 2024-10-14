#ifndef _Orient_H
#define _Orient_H

#include <vector>
#include <set>
#include <deque>
#include <algorithm>

#include <euler.H>
#include <Particle.H>

#include <config_exp.h>

class Component;

extern bool use_cuda;

/** Class to hold energy and angular momentum for use in centering
    a component.  Helper class for class Orient. */
class EL3 
{
public:
  /// Time
  double T;
  /// Mass
  double M;
  /// Binding energy
  double E;			
  /// Angular momentum
  Eigen::Vector3d L;			
  /// Position
  Eigen::Vector3d R;
  
  //! Constructor
  EL3() {}
  
  //! Comparison member for ordering energy
  bool operator<(const EL3& x) const {
    return E < x.E;
  }
  
  //! Debug output
  void debug() const;
};


/** Comparsion function for set<>
 */
struct ltEL3
{
  bool operator()(const EL3& s1, const EL3& s2) const
  {
    return (s1.E < s2.E);
  }
};

/** Class to keep track of orientation 

Log file contents:

<ol>
<li> Time </li>
<li> Energy </li>
<li> Particles used </li>
<li> Current Axis(x) from regression </li>
<li> Current Axis(y) </li>
<li> Current Axis(z) </li>
<li> Last Axis(x) from particles </li>
<li> Last Axis(y) </li>
<li> Last Axis(z) </li>
<li> Current Center(x) from regression </li>
<li> Current Center(y) </li>
<li> Current Center(z) </li>
<li> Current Center(x) from analytic approx </li>
<li> Current Center(y) </li>
<li> Current Center(z) </li>
<li> Last Center(x) from particles </li>
<li> Last Center(y) </li>
<li> Last Center(z) </li>
<li> COM(x) </li>
<li> COM(y) </li>
<li> COM(z) </li>
<li> System COM(x) </li>
<li> System COM(y) </li>
<li> System COM(z) </li>
</ol>
*/
class Orient
{
  typedef pair<double, Eigen::VectorXd> DV;
  
private:
  EL3 t;
  set<EL3, ltEL3> angm;
  deque<DV> sumsA, sumsC;
  int keep, current;
  double damp;
  Eigen::Matrix3d body, orig;
  Eigen::Vector3d axis, center, axis1, center1, center0, cenvel0;
  double sumX, sumX2;
  Eigen::Vector3d sumY, sumY2, sumXY, slope, intercept;
  double sigA, sigC, sigCz;
  
  int many, used;
  double lasttime, deltaT;
  double Ecurr, Elast;
  int Nlast;
  unsigned int oflags;
  unsigned int cflags;
  string logfile;
  bool linear;
  
  vector<double> pos, psa, vel;
  
  void accumulate_cpu(double time, Component* c);
#if HAVE_LIBCUDA==1
  void accumulate_gpu(double time, Component* c);
#endif

public:
  //! Mask flags for Orient type
  enum OrientFlags {AXIS=1, CENTER=2};
  
  //! Mask flags for Orient behavior
  enum ControlFlags {DIAG=1, KE=2, EXTERNAL=4};
  
  //! Constructor
  Orient(int number_to_keep, int target,
	 unsigned orient_flags, unsigned control_flags, string logfile,
	 double dt=0.0, double damping=1.0);
  
  //! Destructor
  ~Orient();
  
  //! Set initial center (zero by default)
  void set_center(double x, double y, double z) {
    center[0] = center0[0] = x;
    center[1] = center0[1] = y;
    center[2] = center0[2] = z;
  }
  
  //! Set initial center velocity (zero by default)
  void set_cenvel(double u, double v, double w) {
    cenvel0[0] = u;
    cenvel0[1] = v;
    cenvel0[2] = w;
  }
  
  //! Write log entry
  void logEntry(double time, Component *c);
  
  //! Set center to move at constant velocity
  void set_linear() { linear = true; }
  
  //! Register phase space of _num_ particles and store angular momentum vector for the lowest _many_ binding energies
  void accumulate(double time, Component* c);
  
  //! Return transformation to new body axes
  Eigen::Matrix3d& transformBody(void) {return body;};
  
  //! Return transformation to original coordinates
  Eigen::Matrix3d& transformOrig(void) {return orig;};
  
  //! Return current axis
  Eigen::Vector3d& currentAxis(void) {return axis;};
  
  //! Return current center
  Eigen::Vector3d& currentCenter(void) {return center;};
  
  //! Return variance of axis determination (linear least squares solution)
  double currentAxisVar(void) {return sigA;}
  
  //! Return variance for center determination (linear least squares solution)
  double currentCenterVar(void) {return sigC;}
  
  //! Return variance for center determination for z component alone
  double currentCenterVarZ(void) {return sigCz;}
  
  //! Return number of particles used
  int currentUsed(void) {return used;};
  
  //! Return energy for disk ang mom
  double currentE(void) {return Ecurr;};
  
};

#endif
