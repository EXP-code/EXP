
#include <values.h>

#include "expand.h"

#include <gaussQ.h>
#include <Sphere.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

Sphere::Sphere(string& line, MixtureSL* m) : SphericalBasis(line, m)
{
  id = "Sphere SL";
				// Defaults
  rmin = 1.0e-3;
  rs = 0.067*rmax;
  numr = 2000;
  cmap = 1;
  diverge = 0;
  dfac = 1.0;
  model_file = "SLGridSph.model";
  cache_file = "SLGridSph.cache";

				// Get initialization info
  initialize();

  SLGridSph::mpi = 1;		// Turn on MPI
  SLGridSph::model_file_name = model_file;
  SLGridSph::sph_cache_name = cache_file + "." + runtag;
  

				// Generate Sturm-Liouville grid
  ortho = new SLGridSph(Lmax, nmax, numr, rmin, rmax, cmap, rs, diverge, dfac);

  setup();
}


void Sphere::initialize()
{
  string val;

  if (get_value("rmin", val))      rmin = atof(val.c_str());
  if (get_value("rs", val))        rs = atof(val.c_str());
  if (get_value("numr", val))      numr = atoi(val.c_str());
  if (get_value("cmap", val))      cmap = atoi(val.c_str());
  if (get_value("diverge", val))   diverge = atoi(val.c_str());
  if (get_value("dfac", val))      dfac = atof(val.c_str());
  if (get_value("modelname", val)) model_file = val.c_str();
  if (get_value("cachename", val)) cache_file = val.c_str();

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
