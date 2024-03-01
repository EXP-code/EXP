#include <limits>

#include "expand.H"

#include <gaussQ.H>
#include <FlatDisk.H>
#include <interp.H>

const std::set<std::string>
FlatDisk::valid_keys = {
  "nmaxfid",
  "numr",
  "rcylmin",
  "rcylmax",
  "mmax",
  "numx",
  "numy",
  "NQDHT",
  "knots",
  "logr",
  "model",
  "biorth",
  "diskconf",
  "cachename"
};

FlatDisk::FlatDisk(Component* c0, const YAML::Node& conf, MixtureBasis* m) :
  PolarBasis(c0, conf, m)
{
  id = "FlatDisk EOF";
				// Defaults
  knots      = 40;
  numr       = 2000;
  rcylmin    = 0.0;
  rcylmax    = 10.0;
  logr       = false;
  is_flat    = true;

				// Get initialization info
  initialize();

  id += "[" + model + "][" + biorth + "]";

  if (myid==0) {
    std::string sep("----    ");
    std::cout << "---- FlatDisk parameters: "
	      << std::endl << sep << "lmax="        << Lmax
	      << std::endl << sep << "nmax="        << nmax
	      << std::endl << sep << "rcylmin="     << rcylmin
	      << std::endl << sep << "rcylmax="     << rcylmax
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


void FlatDisk::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["rcylmin"])   rcylmin    = conf["rcylmin"].as<double>();
    if (conf["rcylmax"])   rcylmax    = conf["rcylmax"].as<double>();
    if (conf["mmax"])      mmax       = conf["mmax"].as<int>();
    if (conf["numr"])      numr       = conf["numr"].as<int>();
    if (conf["knots"])     knots      = conf["knots"].as<int>();
    if (conf["logr"])      logr       = conf["logr"].as<bool>();
    if (conf["model"])     model      = conf["model"].as<std::string>();
    if (conf["biorth"])    biorth     = conf["biorth"].as<std::string>();

    if (conf["mmax"]) Lmax = Mmax = mmax; // Override base-class values
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in FlatDisk: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // Create the BiorthCyl instance
  ortho = std::make_shared<BiorthCyl>(conf);

}

FlatDisk::~FlatDisk(void)
{
  // NADA
}


void FlatDisk::get_dpotl(double r, double z,
			 Eigen::MatrixXd& p,
			 Eigen::MatrixXd& dpr,
			 Eigen::MatrixXd& dpz, int tid)
{
  ortho->get_pot   (p,   r, z);
  ortho->get_rforce(dpr, r, z);
  ortho->get_zforce(dpz, r, z);
}

void FlatDisk::get_potl(double r, double z, Eigen::MatrixXd& p, int tid)
{
  ortho->get_pot(p, r, z);
}

void FlatDisk::get_dens(double r, double z, Eigen::MatrixXd& p, int tid)
{
  ortho->get_dens(p, r, z);
}

void FlatDisk::get_potl_dens(double r, double z, Eigen::MatrixXd& p,
			     Eigen::MatrixXd& d, int tid)
{
  ortho->get_pot (p, r, z);
  ortho->get_dens(d, r, z);
}

