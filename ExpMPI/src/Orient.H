// This may look like C code, but it is really -*- C++ -*-

#ifndef _Orient_H
#define _Orient_H


#include <vector>
#include <deque>
#include <algorithm>

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>

#include <Vector.h>
#include <Particle.H>
#include <Component.H>

/** Class to hold energy and angular momentum for use in centering
    a component.  Helper class for class Orient. */
class EL3 
{
public:
  
  /// Mass
  double M;
  /// Binding energy
  double E;			
  /// Angular momentum
  Vector L;			
  /// Position
  Vector R;

  //! Constructor
  EL3() {
    L.setsize(1, 3);
    R.setsize(1, 3);
  }

  //! Comparison member for ordering energy
  bool operator<(const EL3& x) const {
    return E < x.E;
  }
};

/** Class to keep track of orientation 

Log file contents:

-#	Time
-#	Energy
-#	Particles used
-#	Current Axis(x) from regression
-#	Current Axis(y)
-#	Current Axis(z)
-#	Last Axis(x) from particles
-#	Last Axis(y)
-#	Last Axis(z)
-#	Current Center(x) from analytic approx
-#	Current Center(y)
-#	Current Center(z)
-#	Current Center(x) from regression
-#	Current Center(y)
-#	Current Center(z)
-#	Last Center(x) from particles
-#	Last Center(y)
-#	Last Center(z)
-#	Effective COM(x)
-#	Effective COM(y)
-#	Effective COM(z)
-#	E grad
-#	Delta E
.
*/
class Orient
{
private:
  EL3 t;
  vector<EL3> angm;
  deque<Vector> sumsA, sumsC;
  int keep, current;
  Matrix body, orig;
  Vector axis, center, axis1, center1, center0, cenvel0;
  double sumX, sumX2;
  Vector sumY, sumY2, sumXY, slope;
  double sigA, sigC, sigCz;

  int many, used;
  double lasttime;
  double Ecurr, Elast, Egrad;
  int Nlast;
  unsigned int oflags;
  unsigned int cflags;
  string logfile;
  bool linear;

  vector<double> pos, psa, vel;

  ACG *gen;
  Normal *gauss;

  void write_log(double time, double Egrad, double dE, Component *c);

public:
  //! Mask flags for Orient type
  enum OrientFlags {AXIS=1, CENTER=2};

  //! Mask flags for Orient behavior
  enum ControlFlags {DIAG=1, KE=2, EXTERNAL=4};

  //! Constructor
  Orient(int number_to_keep, int target, double rinit, 
	 unsigned orient_flags, unsigned control_flags, string logfile);

  //! Destructor
  ~Orient();

  //! Set initial center (zero by default)
  void set_center(double x, double y, double z) {
    center[1] = center0[1] = x;
    center[2] = center0[2] = y;
    center[3] = center0[3] = z;
  }

  //! Set initial center velocity (zero by default)
  void set_cenvel(double u, double v, double w) {
    cenvel0[1] = u;
    cenvel0[2] = v;
    cenvel0[3] = w;
  }

  //! Write log entry
  void logEntry(double time, Component *c) { write_log(time, 0.0, 0.0, c); }

  //! Set center to move at constant velocity
  void set_linear() { linear = true; }

  //! Register phase space of _num_ particles and store angular momentum vector for the lowest _many_ binding energies
  void accumulate(double time, Component* c);

  //! Return transformation to new body axes
  Matrix& transformBody(void) {return body;};

  //! Return transformation to original coordinates
  Matrix& transformOrig(void) {return orig;};

  //! Return current axis
  Vector& currentAxis(void) {return axis;};

  //! Return current center
  Vector& currentCenter(void) {return center;};

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
