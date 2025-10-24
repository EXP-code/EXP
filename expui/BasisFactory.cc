#include <algorithm>

#include "YamlCheck.H"
#include "EXPException.H"
#include "BasisFactory.H"
#include "BiorthBasis.H"
#include "FieldBasis.H"
#include "exputils.H"

#ifdef HAVE_FE_ENABLE
#include <cfenv>
#endif

namespace BasisClasses
{
  std::map<Basis::Coord, std::string> Basis::coordLabels =
    { {Basis::Coord::Spherical,   "Spherical"  },
      {Basis::Coord::Cylindrical, "Cylindrical"},
      {Basis::Coord::Cartesian,   "Cartesian"  },
      {Basis::Coord::None,        "None"       } };

  Basis::Basis(const YAML::Node& CONF, const std::string& id)
  {
    // Class name
    //
    name = id;

    // Copy the YAML config
    //
    node = CONF;

    // Complete the initialization
    //
    initialize();
  }

  Basis::Basis(const std::string& confstr, const std::string& id)
  {
    // Assign class name
    //
    name = id;

    // Read the YAML from a string
    //
    try {
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
      conf = node["parameters"];
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing Basis parameters for <"
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
      std::cout << "Basis::factory constructor: found a problem in the YAML config"
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
      if (myid==0) std::cout << "Error parsing force id in Basis::factory"
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Basis::factory: error parsing YAML");
    }
    
    try {
      if ( !name.compare("sphereSL") ) {
	basis = std::make_shared<SphericalSL>(conf);
      }
      else if ( !name.compare("bessel") ) {
	basis = std::make_shared<Bessel>(conf);
      }
      else if ( !name.compare("cylinder") ) {
	basis = std::make_shared<Cylindrical>(conf);
      }
      else if ( !name.compare("flatdisk") ) {
	basis = std::make_shared<FlatDisk>(conf);
      }
      else if ( !name.compare("CBDisk") ) {
	basis = std::make_shared<CBDisk>(conf);
      }
      else if ( !name.compare("slabSL") ) {
	basis = std::make_shared<Slab>(conf);
      }
      else if ( !name.compare("cube") ) {
	basis = std::make_shared<Cube>(conf);
      }
      else if ( !name.compare("field") ) {
	basis = std::make_shared<FieldBasis>(conf);
      }
      else if ( !name.compare("velocity") ) {
	basis = std::make_shared<VelocityBasis>(conf);
      }
      else {
	std::string msg("I don't know about the basis named: ");
	msg += name;
	msg += ". Known types are currently 'sphereSL', 'cylinder', 'flatdisk', 'CBDisk', 'slabSL', 'cube', 'field', and 'velocity'";
	throw std::runtime_error(msg);
      }
    }
    catch (std::exception& e) {
      std::cout << "Error in Basis::factory constructor: " << e.what() << std::endl;
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
    
  std::tuple<std::map<std::string, Eigen::VectorXd>, Eigen::VectorXd>
  Basis::getFieldsCoefs
  (double x, double y, double z, std::shared_ptr<CoefClasses::Coefs> coefs)
  {
    // Python dictonary for return
    std::map<std::string, Eigen::VectorXd> ret;

    // Times for the coefficients
    auto times  = coefs->Times();

    // Initialize the dictionary/map
    auto fields = getFieldLabels(coordinates);
    for (auto s : fields) ret[s].resize(times.size());

    // Make the return dictionary of arrays
    for (int i=0; i<times.size(); i++) {
      set_coefs(coefs->getCoefStruct(times[i]));
      // The field evaluation
      auto v = crt_eval(x, y, z); 
      // Pack the fields into the dictionary
      for (int j=0; j<fields.size(); j++) ret[fields[j]][i] = v[j];
    }

    // An attempt at an efficient return type for the time array
    Eigen::VectorXd T =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(times.data(), times.size());

    // Return the dictionary and the time array
    return {ret, T};
  }

  // Generate coefficients from the accumulated array values
  CoefClasses::CoefStrPtr Basis::makeFromArray(double time)
  {
    make_coefs();
    load_coefs(coefret, time);
    return coefret;
  }

  // Generate coefficients from a phase-space table
  //
  CoefClasses::CoefStrPtr Basis::createFromArray
  (Eigen::VectorXd& m, RowMatrixXd& p, double time, Eigen::Vector3d ctr,
   RowMatrix3d rot, bool roundrobin, bool posvelrows)
  {
    initFromArray(ctr, rot);
    addFromArray(m, p, roundrobin, posvelrows);
    return makeFromArray(time);
  }

  void Basis::setNonInertial(int N, const Eigen::VectorXd& t, const Eigen::MatrixXd& pos)
  {
    // Sanity checks
    if (t.size() < 1)
      throw std::runtime_error("Basis: setNonInertial: no times in time array");

    if (t.size() != pos.rows())
      throw std::runtime_error("Basis::setNonInertial: size mismatch in time and position arrays");

    // Set the data
    Naccel  = N;
    t_accel = t;
    p_accel = pos;
  }

  void Basis::setNonInertial(int N, const std::string& orient)
  {
    std::ifstream in(orient);

    if (not in) {
      throw std::runtime_error("Cannot open Orient file with centering data: " + orient);
    }

    const int cbufsiz = 16384;
    std::unique_ptr<char[]> cbuf(new char [cbufsiz]);
	
    // Look for data and write it while
    // accumlating data for averaging
    Eigen::Vector3d testread;
    double time, dummy;

    std::vector<double> times;
    std::vector<Eigen::Vector3d> centers;

    while (in) {

      in.getline(cbuf.get(), cbufsiz);
      if (in.rdstate() & (ios::failbit | ios::eofbit)) break;

      // Skip comment lines
      //
      if (cbuf[0] == '#') continue;

      std::istringstream line(cbuf.get());

      // Read until current time is reached
      line >> time;		// 
      line >> dummy;
      line >> dummy;
      
      bool allRead = true;
      for (int i=0; i<8; i++) {
	if (line.eof()) allRead = false;
	for (int k; k<3; k++) line >> testread(k);
      }
      if (allRead) {
	times.push_back(time);
	centers.push_back(testread);
      }
    }

    // Repack data
    Naccel = N;
    t_accel.resize(times.size());
    p_accel.resize(times.size(), 3);
    for (int i=0; i<times.size(); i++) {
      t_accel(i) = times[i];
      for (int k=0; k<3; k++) p_accel(i, k) = centers[i](k);
    }
  }

  
  Eigen::Vector3d Basis::currentAccel(double time)
  {
    Eigen::Vector3d ret;

    auto n = t_accel.size();

    // Allow a little bit of buffer in the allowable on-grid range but
    // otherwise force termination
    //
    if ( time < t_accel(0  ) - 0.5*(t_accel(1  ) - t_accel(0  ))  ||
	 time > t_accel(n-1) + 0.5*(t_accel(n-1) - t_accel(n-2)) ) {
      
      std::ostringstream msg;
      msg << "Basis::currentAccel: " << time
	  << " is outside the range of the non-inertial DB ["
	  << t_accel(0) << ", " << t_accel(n-1) << "]";

      throw std::runtime_error(msg.str());
    }
    // Do the quadratic interpolation
    //
    else {
      int imin = 0;
      int imax = std::lower_bound(t_accel.data(), t_accel.data()+n, time) - t_accel.data();

      // Get a range around the current time of approx size Naccel
      //
      imax = std::min<int>(n-1, imax + Naccel/2);
      imin = std::max<int>(imax - Naccel, 0);

      int num = imax - imin + 1;
      Eigen::VectorXd t(num);
      Eigen::MatrixXd p(num, 3);
      for (int i=imin; i<=imax; i++) {
	t(i-imin) = t_accel(i);
	for (int k=0; k<3; k++)	p(i-imin, k) = p_accel(i, k);
      }
      
      for (int k=0; k<3; k++)
	ret(k) = 2.0*std::get<0>(QuadLS<Eigen::VectorXd>(t, p.col(k)).coefs());
    }
    return ret;
  }

}
// END namespace BasisClasses
