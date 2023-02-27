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
  cachename  = "Cyl2d.cache";
  tnext      = 0.0;
  logr       = false;

				// Get initialization info
  initialize();

				// Enable MPI code for more than one node
  // if (numprocs>1) BiorthCyl::mpi = 1;

  std::string cachename = outdir  + cachename + "." + runtag;

  id += ", model=" + model;

				// Generate Sturm-Liouville grid
  ortho = std::make_shared<BiorthCyl>(conf);

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
	      << std::endl << sep << "NO_L0="       << std::boolalpha << NO_M0
	      << std::endl << sep << "NO_L1="       << std::boolalpha << NO_M1
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
    if (conf["cachename"]) cachename  = conf["cachename"].as<std::string>();
    if (conf["dtime"])     dtime      = conf["dtime"].as<double>();
    if (conf["logr"])      logr       = conf["logr"].as<bool>();
    if (conf["model"])     model      = conf["model"].as<std::string>();
    if (conf["biorth"])    biorth     = conf["biorth"].as<std::string>();
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


void Cylinder2d::get_dpotl(double r, double z,
			   Eigen::MatrixXd& p, Eigen::MatrixXd& dpr, Eigen::MatrixXd& dpz, int tid)
{
  ortho->get_pot(p, r, z);
  ortho->get_rforce(dpr, r, z);
  ortho->get_zforce(dpz, r, z);
}

void Cylinder2d::get_potl(double r, double z, Eigen::MatrixXd& p, int tid)
{
  ortho->get_pot(p, r, z);
}

void Cylinder2d::get_dens(double r, double z, Eigen::MatrixXd& p, int tid)
{
  ortho->get_dens(p, r, z);
}

void Cylinder2d::get_potl_dens(double r, double z, Eigen::MatrixXd& p,
			       Eigen::MatrixXd& d, int tid)
{
  ortho->get_pot(p, r, z);
  ortho->get_dens(d, r, z);
}


