
#include <limits>

#include "expand.H"

#include <gaussQ.H>
#include <Cylinder2d.H>
#include <interp.H>

const std::set<std::string>
Cylinder2d::valid_keys = {
  "knots",
  "numr",
  "acyl",
  "cmap",
  "logr",
};

Cylinder2d::Cylinder2d(Component* c0, const YAML::Node& conf, MixtureBasis* m) :
  PolarBasis(c0, conf, m)
{
  id = "Cylinder2d EOF";
				// Defaults
  knots      = 40;
  numr       = 2000;
  cmap       = 1;
  acyl       = 0.01;
  cache_file = "Cyl2d.cache";
  tnext      = 0.0;
  logr       = false;

				// Get initialization info
  initialize();

				// Basis computation logic
  if (dtime>0.0) {
    recompute = true;
    tnext = tnow + dtime;
  }
				// Enable MPI code for more than one node
  if (numprocs>1) SLGridSph::mpi = 1;

  std::string modelname = homedir + model_file;
  std::string cachename = outdir  + cache_file + "." + runtag;

  id += ", model=" + modelname;

				// Generate Sturm-Liouville grid
  ortho = std::make_shared<EmpCyl2D>(Mmax, nmax, knots, numr,
				     rmin, rmax, A, scale, cmap, logr,
				     model, biorth, cachename);

				// Get the min and max expansion radii
  rmin  = ortho->getRmin();
  rmax  = ortho->getRmax();

  if (myid==0) {
    std::string sep("----    ");
    std::cout << "---- Cylinder2d parameters: "
	      << std::endl << sep << "lmax="        << Lmax
	      << std::endl << sep << "nmax="        << nmax
	      << std::endl << sep << "cmap="        << cmap
	      << std::endl << sep << "rmin="        << rmin
	      << std::endl << sep << "rmax="        << rmax
	      << std::endl << sep << "acyl="        << acyl
	      << std::endl << sep << "logr="        << std::boolalpha << logr
	      << std::endl << sep << "NO_L0="       << std::boolalpha << NO_L0
	      << std::endl << sep << "NO_L1="       << std::boolalpha << NO_L1
	      << std::endl << sep << "EVEN_L="      << std::boolalpha << EVEN_L
	      << std::endl << sep << "EVEN_M="      << std::boolalpha << EVEN_M
	      << std::endl << sep << "M0_ONLY="     << std::boolalpha << M0_only
	      << std::endl << sep << "selfgrav="    << std::boolalpha << self_consistent
	      << std::endl;
  }


  setup();
}


void Cylinder2d::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["acyl"])      acyl       = conf["acyl"].as<double>();
    if (conf["numr"])      numr       = conf["numr"].as<int>();
    if (conf["knots"])     knots      = conf["knots"].as<int>();
    if (conf["cmap"])      cmap       = conf["cmap"].as<int>();
    if (conf["diverge"])   diverge    = conf["diverge"].as<int>();
    if (conf["dfac"])      dfac       = conf["dfac"].as<double>();
    if (conf["cachename"]) cache_file = conf["cachename"].as<std::string>();
    if (conf["dtime"])     dtime      = conf["dtime"].as<double>();
    if (conf["logr"])      logr       = conf["logr"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Cylinder2d: "
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

Cylinder2d::~Cylinder2d(void)
{
  // NADA
}


void Cylinder2d::get_dpotl(int lmax, int nmax, double r, Eigen::MatrixXd& p,
		       Eigen::MatrixXd& dp, int tid)
{
  ortho->get_pot(p, r);
  ortho->get_force(dp, r);
}

void Cylinder2d::get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid)
{
  ortho->get_pot(p, r);
}

void Cylinder2d::get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid)
{
  ortho->get_dens(p, r);
}

void Cylinder2d::get_potl_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p,
			   Eigen::MatrixXd& d, int tid)
{
  ortho->get_pot(p, r);
  ortho->get_dens(d, r);
}
