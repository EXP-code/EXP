#include "expand.h"

#include <ExternalForce.H>

ExternalForce::ExternalForce(string& line) : PotAccel(line)
{
				// Do nothing
}

void ExternalForce::get_acceleration_and_potential(vector<Particle>* P)
{
  particles = P;		// "Register" particles
  nbodies = (*particles).size(); // And compute number of bodies

  
  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  determine_acceleration_and_potential();

  MPL_stop_timer();
}


void ExternalForce::determine_coefficients(void)
{
  exp_thread_fork(true);
}

void * ExternalForce::determine_coefficients_thread(void * arg)
{
				// Do nothing
}

void ExternalForce::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);
}



