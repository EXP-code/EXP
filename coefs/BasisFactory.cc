#include <algorithm>

#include <YamlCheck.H>
#include <EXPException.H>
#include <BasisFactory.H>
#include <BiorthBasis.H>
#include <VelocityBasis.H>
#include <exputils.H>

#ifdef HAVE_FE_ENABLE
#include <cfenv>
#endif

namespace BasisClasses
{
  std::map<Basis::Coord, std::string> Basis::coordLabels =
    { {Basis::Coord::Spherical, "Spherical"}, {Basis::Coord::Cylindrical, "Cylindrical"},
      {Basis::Coord::Cartesian, "Cartesian"}, {Basis::Coord::None, "None"} };

  Basis::Basis(const YAML::Node& CONF)
  {
    // Copy the YAML config
    //
    node = CONF;

    // Complete the initialization
    //
    initialize();
  }

  Basis::Basis(const std::string& confstr)
  {
    try {
      // Read the YAML from a string
      //
      node = YAML::Load(confstr);
    }
    catch (const std::runtime_error& error) {
      std::cout << "Basis constructor: found a problem in the YAML config"
		<< std::endl;
      throw;
    }

    // Complete the initialization
    //
    initialize();
  }

  void Basis::initialize()
  {
#ifdef HAVE_FE_ENABLE
    // Flag invalid FP results only, such as 0/0 or infinity - infinity
    // or sqrt(-1).
    //
    // feenableexcept(FE_INVALID);
#endif
  
    // Check whether MPI is initialized
    //
    int flag;
    MPI_Initialized(&flag);
    if (flag) use_mpi = true;
    else      use_mpi = false;
    
    // Fall back sanity (works for me but this needs to be fixed
    // generally)
    //
    if (use_mpi) {
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    }

    // Parameters for force
    //
    try {
      // First get the id name; we know it exists because it has been
      // parsed by the factory
      name = node["id"].as<std::string>();
      // Then . . . 
      conf = node["parameters"];
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing force 'parameters' for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << node                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Basis: error parsing YAML");
    }
    
    // Set coefficient center to zero by default
    //
    coefctr = {0.0, 0.0, 0.0};

    // Default null coordinate type
    //
    coordinates = Coord::None;
  }
  
  Basis::Coord Basis::parseFieldType(std::string coord_type)
  {
    // Find coordinate type
    //
    std::transform(coord_type.begin(), coord_type.end(), coord_type.begin(),
    [](unsigned char c){ return std::tolower(c); });

    Coord ctype;

    if (coord_type.find("cyl") != std::string::npos) {
      ctype = Coord::Cylindrical;
    } else if (coord_type.find("cart") != std::string::npos) {
      ctype = Coord::Cartesian;
    } else if (coord_type.find("none") != std::string::npos) {
      ctype = Coord::None;
    } else {
      ctype = Coord::Spherical;
    }

    return ctype;
  }

  std::shared_ptr<Basis> Basis::factory_string(const std::string& conf)
  {
    YAML::Node node;
    
    try {
      // Read the YAML from a string
      //
      node = YAML::Load(conf);
    }
    catch (const std::runtime_error& error) {
      std::cout << "Basis constructor: found a problem in the YAML config"
		<< std::endl;
      throw;
    }

    // Complete the initialization
    //
    return factory_initialize(node);
  }

  std::shared_ptr<Basis> Basis::factory(const YAML::Node& conf)
  {
    return factory_initialize(conf);
  }


  std::shared_ptr<Basis> Basis::factory_initialize(const YAML::Node& conf)
  {
    std::shared_ptr<Basis> basis;
    std::string name;
    
    // Load parameters from YAML configuration node
    try {
      name = conf["id"].as<std::string>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing force id in BasisFactory"
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("BasisFactory: error parsing YAML");
    }
    
    try {
      if ( !name.compare("sphereSL") ) {
	basis = std::make_shared<SphericalSL>(conf);
      }
      else if ( !name.compare("cylinder") ) {
	basis = std::make_shared<Cylindrical>(conf);
      }
      else if ( !name.compare("flatdisk") ) {
	basis = std::make_shared<FlatDisk>(conf);
      }
      else if ( !name.compare("cube") ) {
	basis = std::make_shared<Cube>(conf);
      }
      else if ( !name.compare("velocity") ) {
	basis = std::make_shared<VelocityBasis>(conf);
      }
      else {
	std::string msg("I don't know about the basis named: ");
	msg += name;
	msg += ". Known types are currently 'sphereSL', 'cylinder' and 'flatdisk'";
	throw std::runtime_error(msg);
      }
    }
    catch (std::exception& e) {
      std::cout << "Error in BasisFactory constructor: " << e.what() << std::endl;
      throw;			// Rethrow the exception?
    }
    
    return basis;
  }

  std::vector<double>
  Basis::operator()(double x1, double x2, double x3, const Coord ctype)
  {
    if (ctype==Coord::Spherical)
      return sph_eval(x1, x2, x3);
    else if (ctype==Coord::Cylindrical)
      return cyl_eval(x1, x2, x3);
    else if (ctype==Coord::Cartesian)
      return crt_eval(x1, x2, x3);
    else {
      return sph_eval(x1, x2, x3);
    };
  }
    
  std::vector<double> Basis::getFields(double x, double y, double z)
  {
    return crt_eval(x, y, z);
  }
    
}
// END namespace BasisClasses
