#include <BasisClient.H>

BasicBasis::BasicBasis(const YAML::Node& CONF) : conf(CONF)
{
  // Nothing yet
}

SphericalSL::SphericalSL(const YAML::Node& CONF) : BasicBasis(CONF)
{
  try {
    if (conf["rs"])        rs         = conf["rs"].as<double>();
    if (conf["numr"])      numr       = conf["numr"].as<int>();
    if (conf["nums"])      nums       = conf["nums"].as<int>();
    if (conf["cmap"])      cmap       = conf["cmap"].as<int>();
    if (conf["diverge"])   diverge    = conf["diverge"].as<int>();
    if (conf["dfac"])      dfac       = conf["dfac"].as<double>();
    if (conf["modelname"]) model_file = conf["modelname"].as<std::string>();
    if (conf["cachename"]) cache_file = conf["cachename"].as<std::string>();
    if (conf["logr"])      logr       = conf["logr"].as<bool>();
    if (conf["plummer"])   plummer    = conf["plummer"].as<bool>();

    if (conf["scale"]) 
      scale = conf["scale"].as<double>();
    else
      scale = 1.0;

    if (conf["rmin"]) 
      rmin = conf["rmin"].as<double>();
    else
      rmin = 0.0;

    if (conf["rmax"]) 
      rmax = conf["rmax"].as<double>();
    else
      rmax = std::numeric_limits<double>::max();

    if (conf["NO_L0"])   NO_L0   = conf["NO_L0"].as<bool>();
    if (conf["NO_L1"])   NO_L1   = conf["NO_L1"].as<bool>();
    if (conf["EVEN_L"])  EVEN_L  = conf["EVEN_L"].as<bool>();
    if (conf["EVEN_M"])  EVEN_M  = conf["EVEN_M"].as<bool>();
    if (conf["M0_ONLY"]) M0_only = conf["M0_ONLY"].as<bool>();
  } 
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing 'force' for Component <"
			   << name << ">: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(235);
  }
}


CylindricalSL::CylindricalSL(const YAML::Node& CONF) : BasicBasis(CONF)
{
  try {
    if (cconf["com"     ])  com_system = cconf["com"     ].as<bool>();
    if (cconf["comlog"  ])     com_log = cconf["comlog"  ].as<bool>();
    if (cconf["timers"  ])      timers = cconf["timers"  ].as<bool>();
  } 
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing 'force' for Component <"
			   << name << ">: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(236);
  }
}


BasisClient::BasisClient(const YAML::Node& conf)
{
  // Load parameters from YAML configuration node
  try {
    id = conf["id"].as<std::string>();
  } 
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing force id in Basis Client"
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(233);
  }

  // Parameters for force
  YAML::Node fconf;

  try {
    fconf = conf["parameters"];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing force 'parameters' for <"
			   << id << ">: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << cforce                << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(234);
  }

  if ( !id.compare("sphereSL") ) {
    basis = std::make_shared<SphericalSL>(fconf);
  }
  else if ( !id.compare("cylinder") ) {
    basis = std::make_shared<CylindricalSL>(fconf);
  }
  else {
    std::string msg("I don't know about the basis: ");
    msg += id;
    throw std::runtime_error(msg);
  }

}
