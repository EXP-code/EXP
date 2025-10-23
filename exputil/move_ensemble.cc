#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include <time.h>
#include "numerical.H"

#include "phase.H"

Ensemble Ensemble::Step_Positions(double dt)
{
  int i;
  Ensemble newstars(Nstars);
  
  for (i=0; i<Nstars; i++)
    {
      newstars[i] = stars[i];
      newstars[i].x += stars[i].v*dt;
    }
  
  return newstars;
}


Ensemble Ensemble::Step_Velocities(double dt, Eigen::Vector3d *f)
{
  int i;
  Ensemble newstars(Nstars);
  
  for (i=0; i<Nstars; i++)
    {
      newstars[i] = stars[i];
      newstars[i].v +=  f[i]*dt;
    }
  
  return newstars;
}	



/*
  allow certain stars to be held fixed
*/

Ensemble Ensemble::integrate_to(double tfin, freezing_function frozen)
{
  int i;
  Ensemble newstars(Nstars);
  
  for (i=0; i<Nstars; i++)
    {
      stars[i].work = 0.0;
      if (i%(Nstars/10) == 1)
	{
	  time_t realtime;
	  time(&realtime);
	  std::cerr << ctime(&realtime);
	  std::cerr << "integrating star #" << i
		    << " of " << Nstars <<  "    ";
	}
      if (! (int) frozen(stars[i]))     // is the star frozen?
	{
	  newstars.stars[i] = stars[i].integrate_to(tfin);
	}
      else
	{
	  newstars.stars[i] = stars[i];
	}
      if (i%(Nstars/10) == 1) 
	std::cerr << "accuracy = "
		  << newstars.stars[i].Accuracy(stars[i]) << std::endl;
    }
  newstars.t = tfin;
  
  return newstars;
}




/*
  integrate all stars
*/


Ensemble Ensemble::integrate_to(double tfin)
{
  int i;
  Ensemble newstars(Nstars);
  
  for (i=0; i<Nstars; i++)
    {
      stars[i].work = 0.0;
      if (Nstars < 200 || i%(Nstars/100) == 1)
	{
	  time_t realtime;
	  time(&realtime);
	  std::cerr << ctime(&realtime);
	  std::cerr << "integrating star #" << i 
		    << " of " << Nstars << "       ";
	}
      newstars.stars[i] = stars[i].integrate_to(tfin);
      if (Nstars < 20 || i%(Nstars/10) == 1)
	std::cerr << "accuracy = " <<
	  newstars.stars[i].Accuracy(stars[i]) << std::endl;
    }
  newstars.t = tfin;
  
  return newstars;
}
