#include <iostream>
#include <cmath>
#include <cstdlib>

#include <time.h>
#include <numerical.H>

#include "phase.H"


/*
  default constructor. Sets everything to zero
*/


Ensemble::Ensemble(void)
{
  Nstars = 0;
  t = 0.0;
}


/*
  construct an ensemble of n stars.
*/


Ensemble::Ensemble(int n)
{
  Nstars = n;
  t = 0.0;
  stars.resize(Nstars);
}



/*
  copy constructor
*/

Ensemble::Ensemble(Ensemble &e)
{
  Nstars = e.Nstars;
  t = e.t;
  stars.resize(Nstars);
  for (int i=0; i<Nstars; i++) stars[i] = e.stars[i];
}


/*
  destructor: free the stars if they exist. If they don't
  exist, print a little warning.
*/

Ensemble::~Ensemble(void)
{
  // Nothing
}


/*
  set the ensemble size. This can be used to initialise
  or to change the size.
*/

void Ensemble::setsize(int n)
{
  Nstars = n;
  std::cout << "Nstars = " << Nstars << std::endl;
  stars.resize(Nstars);
  puts("done allocating stars");
}






/*
  assignment operator. Copy everything from one 
  ensemble to another.
*/


Ensemble &Ensemble::operator=(Ensemble &e)
{
  int i;
  
  if (Nstars && Nstars!=e.Nstars)
    {
      puts("error in Ensemble=");
      exit(0);
    }
  
  Nstars = e.Nstars;
  t = e.t;
  stars.resize(Nstars);

  for (i=0; i<Nstars; i++) stars[i] = e.stars[i];
  
  return *this;
}

Phase &Ensemble::operator[](int i)
{
  if (i<0 || i>=Nstars) 
    {
      puts("subscript out of range in Ensemble[]");
      exit(0);
    }
  
  return stars[i];
}








/* 
   various ways of manipulating an ensemble: scaling, 
   rotation, etc.	
*/




void Ensemble::scale_masses(double ms)
{
  int i;
  
  for (i=0; i<Nstars; i++) stars[i].Mass() *= ms;
}




void Ensemble::scale_positions(double xs)
{
  int i;
  
  for (i=0; i<Nstars; i++) stars[i].Position() *= xs;
}


void Ensemble::scale_speeds(double vs)
{
  int i;
  
  for (i=0; i<Nstars; i++) stars[i].Velocity() *= vs;
}








/*
  set the system time
*/

void Ensemble::settime(double t)
{
  int i;
  
  t = t;
  for (i=0; i<Nstars; i++) stars[i].Time() = t;
}	


void Ensemble::set_masses(double m)
{
  int i;
  
  for (i=0; i<Nstars; i++) stars[i].Mass() = m;
}	






/*
  rotate positions through an angle, theta, about the
  Z axis.
*/


void Ensemble::rotate_view(double theta)
{
  int i;
  
  for (i=0; i<Nstars; i++) stars[i] = stars[i].rotate_view(theta);
}









/* 
   convert velocities to a rotating frame.
*/


void Ensemble::rotate_frame(double omega)
{
  int i;
  
  for (i=0; i<Nstars; i++) stars[i] = stars[i].rotate_frame(omega);
}

void Ensemble::translate(Eigen::Vector3d &dx)
{
  for (int i=0; i<Nstars; i++) stars[i].Position() -= dx;
}
