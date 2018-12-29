
#include <values.h>

#include "expand.h"

#include <gaussQ.h>
#include <Sphere.H>

Sphere::Sphere(const YAML::Node& conf, MixtureBasis* m) : SphericalBasis(conf, m)
{
  id = "Sphere SL";
				// Defaults
  rmin = min<double>(1.0e-6, 1.0e-6*rmax);
  rs = 0.067*rmax;
  numr = 2000;
  cmap = 1;
  diverge = 0;
  dfac = 1.0;
  model_file = "SLGridSph.model";
  cache_file = "SLGridSph.cache";

				// Get initialization info
  initialize();

				// Enable MPI code for more than one node
  if (numprocs>1) SLGridSph::mpi = 1;

  SLGridSph::model_file_name = homedir + model_file;
  SLGridSph::sph_cache_name  = outdir  + cache_file + "." + runtag;
  

				// Generate Sturm-Liouville grid
  ortho = new SLGridSph(Lmax, nmax, numr, rmin, rmax, true,
			cmap, rs, diverge, dfac);

  setup();
}


void Sphere::initialize()
{
  try {
    if (conf["rmin"])      rmin       = conf["rmin"].as<double>();
    if (conf["rs"])        rs         = conf["rs"].as<double>();
    if (conf["numr"])      numr       = conf["numr"].as<int>();
    if (conf["cmap"])      cmap       = conf["cmap"].as<int>();
    if (conf["diverge"])   diverge    = conf["diverge"].as<int>();
    if (conf["dfac"])      dfac       = conf["dfac"].as<double>();
    if (conf["modelname"]) model_file = conf["modelname"].as<std::string>();
    if (conf["cachename"]) cache_file = conf["cachename"].as<std::string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Sphere: "
			   << error.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }
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
