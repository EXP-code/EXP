// This may look like C code, but it is really -*- C++ -*-

#ifndef _Orient_H
#define _Orient_H


#include <vector>
#include <deque>
#include <algorithm>

#include <Vector.h>

#include <Particle.H>

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

/** Class to keep track of orientation */
class Orient
{
private:
  EL3 t;
  vector<EL3> angm;
  deque<Vector> sumsA, sumsC;
  int keep, current;
  Matrix body, orig;
  Vector axis, center, axis1, center1, center0;
  double sumX, sumX2;
  Vector sumY, sumY2, sumXY, slope;
  double sigA, sigC, sigCz;

  int many, used;
  double Ecurr, Elast, Egrad;
  int Nlast;
  int flags;
  string logfile;
  bool verbose;

public:
  //! Mask flags for Orient behavior
  enum OrientFlags {AXIS=1, CENTER=2};

  //! Constructor
  Orient(int number_to_keep, int target, double rinit, int flags, 
	 string logfile, bool verbose=false);

  //! Register phase space of <num> particles and store angular momentum vector for the lowest <many> binding energies
  void accumulate(double time, vector<Particle> *p, double* com);

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
