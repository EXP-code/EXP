
#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#define DENSITY
// #define SELECTOR

#include <values.h>

#include "expand.h"

#include <gaussQ.h>
#include <Sphere.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

Sphere::Sphere(string&line) : SphericalBasis(line)
{
				// Defaults
  rmin = 1.0e-3;
  numr = 1000;

				// Get initialization info
  initialize();

  SLGridSph::mpi = 1;		// Turn on MPI

				// Generate Sturm-Liouville grid
  ortho = new SLGridSph(Lmax, nmax, numr, rmin, rmax);

  setup();
}


void Sphere::initialize()
{
  string val;

  if (get_value("rmin", val)) rmin = atof(val.c_str());
  if (get_value("numr", val)) numr = atoi(val.c_str());

}

Sphere::~Sphere(void)
{
  delete ortho;
}


void Sphere::get_dpotl(int lmax, int nmax, double r, Matrix& p, Matrix& dp,
		       int tid)
{
  ortho->get_pot(p, r);
  ortho->get_force(dp, r);
}

void Sphere::get_potl(int lmax, int nmax, double r, Matrix& p, int tid)
{
  ortho->get_pot(p, r);
}

void Sphere::get_dens(int lmax, int nmax, double r, Matrix& p, int tid)
{
  ortho->get_dens(p, r);
}

void Sphere::get_potl_dens(int lmax, int nmax, double r, Matrix& p, Matrix& d,
			   int tid)
{
  ortho->get_pot(p, r);
  ortho->get_dens(d, r);
}
