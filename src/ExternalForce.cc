#include "expand.H"

#include <ExternalForce.H>

ExternalForce::ExternalForce(const YAML::Node& conf) : PotAccel(0, conf)
{
				// Do nothing
}

void ExternalForce::get_acceleration_and_potential(Component *C)
{
  cC = C;
  nbodies = cC->Number();	// And compute number of bodies

  
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
  print_timings(id + ": coefficient timings");
}

void * ExternalForce::determine_coefficients_thread(void * arg)
{
				// Do nothing
  return (NULL);
}

void ExternalForce::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif
  exp_thread_fork(false);
  print_timings(id + ": acceleration timings");
}

void ExternalForce::print_divider(void)
{
  if (myid) return;

  char c = cout.fill('-');
  cout << setw(80) << "-" << endl;
  cout.fill(c);
}

#if HAVE_LIBCUDA==1
void ExternalForce::getParticlesCuda(Component *c)
{
  if (use_cuda and not cudaAware()) {
    if (not comp->fetched[c]) {
      c->CudaToParticles();
      comp->fetched[c] = true;
    }
  }
}
#endif

