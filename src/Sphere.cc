
#include <limits.h>

#include <boost/make_shared.hpp>

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
  tnext = dtime = 0.0;
  recompute  = false;

				// Get initialization info
  initialize();

				// Basis computation logic
  if (dtime>0.0) {
    recompute = true;
    tnext = tnow + dtime;
  }
				// Enable MPI code for more than one node
  if (numprocs>1) SLGridSph::mpi = 1;

  SLGridSph::model_file_name = homedir + model_file;
  SLGridSph::sph_cache_name  = outdir  + cache_file + "." + runtag;
  

				// Generate Sturm-Liouville grid
  ortho = boost::make_shared<SLGridSph>(Lmax, nmax, numr, rmin, rmax, true,
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
    if (conf["dtime"])     dtime      = conf["dtime"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Sphere: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}

Sphere::~Sphere(void)
{
  // NADA
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

void Sphere::make_model()
{
  // Made a mass model
  //
  double Rmin = rmin;
  double Rmax = component->rtrunc;
  double dr   = (Rmax - Rmin)/numr;

  std::vector<double> histo(numr, 0.0);

  bool not_done = true;

  while (not_done) {

    for (int i=0; i<component->Particles().size(); i++) {

      double rr = 0.0;
      for (int k=0; k<3; k++) {
	double pos = component->Pos(i, k, Component::Local | Component::Centered);
	rr += pos*pos;
      }
      rr = sqrt(rr);
      
      if (rr < Rmax) {
	int id = (rr - Rmin)/dr;
	histo[id] += component->Part(i)->mass;
      }
    }
    
    // All processes get complete mass histogram
    //
    MPI_Allreduce(MPI_IN_PLACE, &histo[0], numr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Check for zero mass bins
    //
    int n = 0;
    for (; n<numr; n++) {
      if (histo[n]<=0.0) break;
    }

    // No zero mass bins
    //
    if (n==numr) {
      not_done = false;
    }
    // Found a zero mass bin
    //
    else {
      double lastR = Rmax;
      Rmax = Rmin + dr*(n-1);
      dr = (Rmax - Rmin)/numr;

      // Diagnostic warning
      //
      if (myid==0 and Rmax/lastR > 2.0) {
	std::cout << "Sphere::make_model at T=" << tnow << " has small radius "
		  << Rmax << "/" << lastR << std::endl;
      }
    }
  }

  std::vector<double> R(numr), D(numr), M(numr), P(numr), P1(numr);

  for (int i=0; i<numr; i++) {

    double r1 = Rmin + dr*i;	// Bin edges
    double r2 = Rmin + dr*(i+1);
    R[i] = rmin + dr*(0.5+i);	// Bin center

				// Bin density
    D[i] = histo[i] / (4.0*M_PI/3.0 * (pow(r2, 3.0) - pow(r1, 3.0)));

    if (i==0) M[i] = 0.0;	// Consistent cumulative mass at edges
    else M[i] = M[i-1] + 2.0*M_PI*dr*(D[i]*r2*r2 + D[i-1]*r1*r1);
    
				// Outer potential
    if (i==0) P1[i] = 4.0*M_PI*dr*D[i]*r2;
    else P1[i] = P1[i-1] + 2.0*M_PI*dr*(D[i]*r2 + D[i-1]*r1);
  }

				// Full potential
  for (int i=numr-1; i>0; i--) 
    P[i] = M[i]/R[i] + P1[numr-1] - 0.5*(P1[i] + P1[i-1]);
  P[0] = M[0]/R[0] + P1[numr-1] - 0.5*P1[0];

  for (int i=0; i<numr; i++) P[i] *= -1.0;

  // Create a new spherical model
  //
  SphModTblPtr mod = boost::make_shared<SphericalModelTable>(numr, &R[0]-1, &D[0]-1, &M[0]-1, &P[0]-1);

  // Regenerate Sturm-Liouville grid
  //
  ortho = boost::make_shared<SLGridSph>(Lmax, nmax, numr, Rmin, Rmax, mod, false);

  // Update time trigger
  //
  tnext = tnow + dtime;
}
