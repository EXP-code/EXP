
#include <iostream>
#include <cstdlib>
#include <assert.h>
#include <cmath>
#include "numerical.H"
#include "phase.H"


/* 
   Predictor-Corrector integrator for phase-space class
*/


/* class for storage of previous values for predictor-corrector */

class PC_table
{
public:
  Phase p1, p2, p3;
  Eigen::Vector3d f1, f2, f3;
};

static PC_table PC;



/* 
   Build table of previous points for predictor-corrector. 
   Uses the normal integrator.
*/


void Phase::set_PC(double dt, int n_subdivide)
{
  int i;
  double h;
  
  h = dt/((double) n_subdivide);
  PC.p1 = *this;
  
  for (i=1; i<=n_subdivide; i++)
    {
      PC.p1 = PC.p1.integrate_to(t - i*h);
    }
  PC.f1 = Force(PC.p1.t, PC.p1.x, PC.p1.v);
  
  PC.p2 = PC.p1;
  
  for (i=1; i<=n_subdivide; i++)
    {
      PC.p2 = PC.p2.integrate_to(t-dt-i*h );
    }
  PC.f2 = Force(PC.p2.t, PC.p2.x, PC.p2.v);
  
  PC.p3 = PC.p2;
  
  for (i=1; i<=n_subdivide; i++)
    {
      PC.p3 = PC.p3.integrate_to(t-2.0*dt-i*h);
    }
  PC.f3 = Force(PC.p3.t, PC.p3.x, PC.p3.v);
  
  
}


Phase Phase::PC_advance(double dt)
{
  static Phase pred, next;
  static Eigen::Vector3d f, fp;
  
  
  /* calculate force at current position */
  
  f = Force(t, x, v);
  
  
  
  /* do predictor step */
  
  pred.x = x + (55.0*v - 59.0*PC.p1.v + 37.0*PC.p2.v 
		- 9.0*PC.p3.v)*dt/24.0;
  pred.v = v + (55.0*f - 59.0*PC.f1 + 37.0*PC.f2 - 9.0*PC.f3)*dt/24.0;
  pred.t = t + dt;
  
  
  
  
  /* calculate force at predicted position */
  
  fp = Force(pred.t, pred.x, pred.v);
  
  
  
  
  /* now do corrector step */
  
  next.x = x + (9.0*pred.v + 19.0*v - 5.0*PC.p1.v + PC.p2.v)*dt/24.0;
  next.v = v + (9.0*fp + 19.0*f - 5.0*PC.f1 + PC.f2)*dt/24.0;
  next.t = t+ dt;
  
  
  
  
  /* update table of previous values */
  
  PC.p3 = PC.p2;
  PC.p2 = PC.p1;
  PC.p1 = *this;
  PC.f3 = PC.f2;
  PC.f2 = PC.f1;
  PC.f1 = f;
  
  return next;
}















