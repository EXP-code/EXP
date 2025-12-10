#include <iostream>
#include <iomanip>

#include "CenterFile.H"
#include "localmpi.H"

// CenterFile contructor
//
CenterFile::CenterFile(const YAML::Node& conf)
{
  std::string name, type;

  try {
    name = conf["file"].as<std::string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0)
      std::cout << __FILE__ << ": " << __LINE__ << std::endl
		<< "CenterFile error finding center file name: "
		<< error.what() << std::endl
		<< std::string(60, '-') << std::endl
		<< "Config node"        << std::endl
		<< std::string(60, '-') << std::endl
		<< conf                 << std::endl
		<< std::string(60, '-') << std::endl;
    throw std::runtime_error("CenterFile: missing <file> parameter");
  }

  try {
    type = conf["type"].as<std::string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0)
      std::cout << "CenterFile error finding center type: "
		<< error.what() << std::endl
		<< std::string(60, '-') << std::endl
		<< "Config node"        << std::endl
		<< std::string(60, '-') << std::endl
		<< conf
		<< std::string(60, '-') << std::endl;
    throw std::runtime_error("CenterFile: bad center <type>");
  }
  
  // Convert type to upper case
  //
  std::for_each(type.begin(), type.end(), [](char & c){ c = std::toupper(c); });

  // Select start column for center values
  //
  int start = 0;
  if      (type.find("EJ" ) == 0) { start = 16; }
  else if (type.find("COM") == 0) { start = 19; }
  else {
    if (myid==0)
      std::cout << "CenterFile error parsing center 'type': "
		<< "fouund <" << type << "> but expected either "
		<< "'COM' or 'EJ'" << std::endl
		<< std::string(60, '-') << std::endl;
    throw std::runtime_error("CenterFile: error parsing type");
  }

  // Open input file
  //
  std::ifstream in(name);

  // Check for open file
  //
  if (in.is_open()) {

    std::string line;

    // Get a line as a string
    //
    while (std::getline(in, line)) {
      // Put line into a stringstream
      //
      std::istringstream iss(line); 

      // Collect the fields into a vector
      //
      std::vector<double> fields;

      double x;
      while(iss >> x) fields.push_back(x);

      // Assign the time and center data
      //
      if (fields.size() >= start+3) {
	time.push_back(fields[0]);
	data.push_back({fields[start-1], fields[start], fields[start+1]});
      }
    }
  } else {
    if (myid==0)
      std::cout << "CenterFile error opening center file <" << name << ">"
		<< std::endl
		<< std::string(60, '-') << std::endl;
    throw std::runtime_error("CenterFile: error opening file");
  }

  // Some chatty output
  //
  if (myid==0) {
    std::cout << "---- CenterFile: opened and read file <" << name << "> "
	      << "and read " << type << " center columns between times "
	      << "[" << time.front() << ", " << time.back() << "]"
	      << std::endl;
  }
}

// CenterFile interpolator
//
std::array<double, 3> CenterFile::operator()(double T)
{
  // Prints error message and quits if out of bounds
  //
  if (T < time.front() or T > time.back()) {
    if (myid==0)
      std::cout << "CenterFile range error: T="
		<< T << " is off grid"<< std::endl
		<< std::string(60, '-') << std::endl;
    throw std::runtime_error("CenterFile::operator(): range error");
  }
  
  // Check limiting cases
  //
  if (T==time.front()) return data.front();
  if (T==time.back() ) return data.back();

  // Otherwise: binary search for location
  //
  size_t i = lower_bound(time.begin(), time.end(), T) - time.begin();

  // Set up for linear interpolation
  //
  double lo = time[i-1];
  double hi = time[i];

  double a = (hi - T)/(hi - lo);
  double b = (T - lo)/(hi - lo);

  // Perform the interpolation
  //
  std::array<double, 3> ret;
  for (int k=0; k<3; k++) ret[k] = a * data[i-1][k] + b * data[i][k];

  return ret;
}
