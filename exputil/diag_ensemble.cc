#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <time.h>
#include "numerical.H"

#include "phase.H"


double Ensemble::total_Energy(void)
{
  return total_Potential() + total_Kinetic();
}


double Ensemble::total_Kinetic(void)
{
  int i;
  double ktot=0.0;
  
  for (i=0; i<Nstars; i++) 
    ktot += 0.5*stars[i].Mass()*stars[i].Velocity().dot(stars[i].Velocity());
  
  return ktot;
}


double Ensemble::total_Potential(void)
{
  int i;
  double ptot=0.0;
  
  /* factor 1/2 accounts for self-energy */
  
  for (i=0; i<Nstars; i++) 
    ptot += 0.5*stars[i].Mass()*
      (*Phase::Potential)(stars[i].Time(), 
			  stars[i].Position());
  
  return ptot;
}





double Ensemble::Virial(void)
{
  return 2.0*total_Kinetic() + total_Potential();
}




Eigen::Vector3d Ensemble::CM_Position(void)
{
  static Eigen::Vector3d x_cm;
  double mtot;
  int i;
  
  x_cm.setZero();
  for (i=0; i<Nstars; i++) 
    {
      x_cm += stars[i].Mass()*stars[i].Position();
      mtot += stars[i].Mass();
    }
  x_cm /= mtot;
  
  return x_cm;
}





Eigen::Vector3d Ensemble::CM_Velocity(void)
{
  static Eigen::Vector3d v_cm;
  double mtot;
  int i;
  
  v_cm.setZero();
  for (i=0; i<Nstars; i++) 
    {
      v_cm += stars[i].Mass()*stars[i].Velocity();
      mtot += stars[i].Mass();
    }
  v_cm /= mtot;
  
  return v_cm;
}


Eigen::Vector3d Ensemble::total_Angular_Momentum(void)
{
  int i;
  static Eigen::Vector3d Jtot;
  
  Jtot.setZero();
  
  for (i=0; i<Nstars; i++) Jtot += stars[i].Angular_Momentum();
  
  return Jtot;
}






/* 
   get second moments of all stars within a limiting 
   radius rmax.
   
   M_{i,j} = \sum_{r<rmax} m_{i,j} x_i x_j
   
   
*/

Eigen::Matrix3d Ensemble::Moment_Tensor(double rmax)
{
  Eigen::Matrix3d M;
  double r, mtot;
  Eigen::Vector3d x, x0;
  
  M.setZero();
  x0.setZero();
  mtot = 0.0;
  
  // find CM position of stars within rmax of origin
  
  for (int n=0; n<Nstars; n++)
    {
      x = stars[n].x;
      if (sqrt(x.dot(x))>rmax) continue;
      x0 += stars[n].m*x;
      mtot += stars[n].m;
    }
  
  if (mtot==0.0) return M;
  x0 /= mtot;
  
  
  // find moment tensor of stars within rmax of CM position
  
  mtot = 0.0;
  
  for (int n=0; n<Nstars; n++)
    {
      x = stars[n].x-x0;
      r = sqrt(x.dot(x));
      if (r > rmax) continue;
      mtot += stars[n].m;
      for (int i=0; i<3; i++)
	{
	  for (int j=0; j<3; j++)
	    {
	      M(i, j) += stars[n].m*x[i]*x[j];
	    }
	}
    }
  
  
  if (mtot>0.0) M /= mtot;
  
  return M;
}


/*
  get the principal axes of the moment tensor
*/


Eigen::Vector3d Ensemble::Principal_Axes(double rmax,
					 Eigen::Matrix3d &directions)
{
  Eigen::Matrix3d M;
  static Eigen::Vector3d PA;
  Eigen::Vector3d V;
  V.setZero();
  
  M = Moment_Tensor(rmax);
  
  if (M.trace() == 0.0) 
    {
      PA.setZero();
      return PA;
    }
  
  Eigen::EigenSolver<Eigen::Matrix3d> es(directions);

  V = es.eigenvalues().real();
  PA = V;
  return PA;
  
  
}


Eigen::Matrix3d Ensemble::Inertia_Tensor(double rmax)
{
  double r2;
  Eigen::Matrix3d M, I;
  
  M = Moment_Tensor(rmax);
  
  r2 = M.trace();
  
  I = -M;
  for (int k=0; k<3; k++) I(k, k) += r2;
  
  return I;
}



Eigen::Vector3d Ensemble::Solid_Body_Frequency(double rmax)
{
  Eigen::Matrix3d I, Iinv;
  Eigen::Vector3d J;
  
  I = Inertia_Tensor(rmax);
  J.setZero();
  double mtot = 0.0;
  
  for (int i=0; i<Nstars; i++)
    {
      if (stars[i].x.dot(stars[i].x) < rmax*rmax)
	{
	  J += stars[i].m * stars[i].x.cross(stars[i].v);
	  mtot += stars[i].m;
	}
    }
  J /= mtot;
  
  Iinv = I.inverse();
  
  return Iinv*J;
}


