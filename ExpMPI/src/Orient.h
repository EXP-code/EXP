// This may look like C code, but it is really -*- C++ -*-

#include <vector>
#include <deque>
#include <algorithm>

// #include <function.h>

#include <Vector.h>

/** Class to hold energy and angular momentum */
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

  EL3() {
    L.setsize(1, 3);
    R.setsize(1, 3);
  }

  /// For ordering energy
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
  Vector axis, center, axis1, center1;
  double sumX, sumX2;
  Vector sumY, sumY2, sumXY, slope;
  double sigA, sigC, sigCz;

  int many, used;
  double Ecurr, Elast, Egrad;
  int Nlast;

public:

  /// Constructor
  Orient(int number_to_keep, int target, double rinit);

  /// Register phase space of <num> particles and store angular momentum vector for the lowest <many> binding energies
  void accumulate(int* c, double* mass, double* pot,
		  double* x, double* y, double* z,
		  double* u, double* v, double* w,
		  double* com, int num);

  /// Return transformation to new body axes
  Matrix& transformBody(void) {return body;};

  /// Return transformation to original coordinates
  Matrix& transformOrig(void) {return orig;};

  /// Return current center
  Vector& currentAxis(void) {return axis;};

  /// Return current center
  Vector& currentCenter(void) {return center;};

  /// Return variances
  double currentAxisVar(void) {return sigA;}
  double currentCenterVar(void) {return sigC;}
  double currentCenterVarZ(void) {return sigCz;}

  /// Return number of particles used
  int currentUsed(void) {return used;};

  /// Return energy for disk ang mom
  double currentE(void) {return Ecurr;};

};
