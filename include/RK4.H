#ifndef _RK4_H_
#define _RK4_H_

#include <Eigen/Eigen>

namespace ODE
{
  //! RK4 takes one 4th-order Runge-Kutta step
  class RK4
  {

  public:

    //! Reference to an Eigen vector
    using Vref = Eigen::Ref<Eigen::VectorXd>;

    //! Functoid defintion for the force
    using Force = std::function<Vref(double, Vref)>;

  private:

    //! Force instance
    Force force;
  
  public:

    //! Constructor defining the force function
    RK4(Force force) : force(force) {}

    //! Take a single RK4 step
    double step(Vref x, double t, double dt)
    {
      Eigen::VectorXd f0 = force(t, x), xt;
      Eigen::VectorXd f1 = force(t + 0.5*dt, (xt=x + 0.5*dt*f0));
      Eigen::VectorXd f2 = force(t + 0.5*dt, (xt=x + 0.5*dt*f1));
      Eigen::VectorXd f3 = force(t + dt, (xt=x + dt*f2));
      
      x += dt * ( f0 + 2.0*f1 + 2.0*f2 + f3 ) / 6.0;
      
      return t + dt;
    }
  };
}

#endif // _RK4_H_
