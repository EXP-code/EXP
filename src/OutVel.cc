#include <filesystem>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "expand.H"

#include "OutVel.H"

const std::set<std::string>
OutVel::valid_keys = {
  "modelname",
  "nint",
  "nintsub",
  "name",
  "dof",
  "rmapping",
  "rmin",
  "rmax",
  "ascl",
  "delta",
  "lmax",
  "nmax",
  "model"
};

OutVel::OutVel(const YAML::Node& conf) : Output(conf)
{
  // Defaults
  //
  nint    = 10;
  nintsub = std::numeric_limits<int>::max();
  tcomp   = NULL;
  dof     = 3;

  // Retrieve parameters
  //
  initialize();

  // Target output file
  //
  outfile = outdir + "velcoef." +  tcomp->name + "." + runtag;

  // Check for valid model type
  //
  std::vector<std::string> modelTypes { "file", "expon" };
  auto it = std::find(modelTypes.begin(), modelTypes.end(), model);
  if (it == modelTypes.end()) {
    std::ostringstream sout;
    sout << "OutVel: found type <" << model << ">.  Must be one of ";
    for (auto v : modelTypes) sout << " " << v;
    throw std::runtime_error(sout.str());
  }

  // Check dof value
  //
  if (dof!=2 and dof!=3) {
    std::ostringstream sout;
    sout << "OutVel: found " << dof << " for dof.  Must be 2 or 3.";
    throw std::runtime_error(sout.str());
  }

  // Create the basis
  //
  YAML::Node node;
  node["parameters"] = conf;

  // Remove OutVel only keys
  //
  for (auto s : {"nint", "nintsub", "name"}) {
    if (node["parameters"][s]) node["parameters"].remove(s);
  }
  
  basis = std::make_shared<BasisClasses::VelocityBasis>(node);

  // Create the coefficient container based on the dimensionality.
  // Currently, these are spherical and polar grids.
  //
  if (dof==2)
    coefs = std::make_shared<CoefClasses::CylFldCoefs>();
  else
    coefs = std::make_shared<CoefClasses::SphFldCoefs>();
}

void OutVel::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["dof"])          dof      = conf["dof"  ].as<int>();
    if (conf["nint"])         nint     = conf["nint" ].as<int>();

    if (conf["nintsub"]) {
      nintsub  = conf["nintsub"].as<int>();
      if (nintsub <= 0) nintsub = 1;
    }

    if (conf["model"])		// Model type is a required parameter
      model    = conf["model"].as<std::string>();
    else {
      std::string message = "OutVel: no model specified. Please specify "
	"either 'file' with the model 'modelname' or 'expon' for the\n"
	"exponential disk model (i.e. Laguerre polynomials)";
      throw std::runtime_error(message);
    }

    // Search for desired component
    //
    if (conf["name"]) {
      std::string tmp = conf["name"].as<std::string>();
      for (auto c : comp->components) {
	if (!(c->name.compare(tmp))) tcomp  = c;
      }
    }

    // Success is required!
    //
    if (!tcomp) {
      std::string message = "OutVel: no component to trace. Please specify "
	"the component name using the 'name' parameter.";
      throw std::runtime_error(message);
    }

    if (conf["modelname"]) modelname = conf["modelname"].as<std::string>();

  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutVel: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("OutVel::initialize: error parsing YAML");
  }
}

void OutVel::Run(int n, int mstep, bool last)
{
  // Skip this master step
  //
  if (n % nint != 0 && !last) return;

  // Skip this sub step
  //
  if (mstep < std::numeric_limits<int>::max() and mstep % nintsub != 0) return;

  // Don't duplicate times
  //
  if (tnow <= prev) return;

  // Record current time
  //
  prev = tnow;

  // Zero coefficients
  //
  basis->reset_coefs();
    
  // Compute coefficients
  //
  PartMapItr ibeg = tcomp->Particles().begin();
  PartMapItr iend = tcomp->Particles().end();

  constexpr std::complex<double> I(0.0, 1.0);

#pragma omp parallel
  {
    for (PartMapItr it=ibeg; it!=iend; ++it) {
      
      // Execute block for one particle at a time...
#pragma omp single nowait
      {
	double M = it->second->mass;

	double x = it->second->pos[0];
	double y = it->second->pos[1];
	double z = it->second->pos[2];
	
	double u = it->second->vel[0];
	double v = it->second->vel[1];
	double w = it->second->vel[2];
	
	basis->accumulate(M, x, y, z, u, v, w);
      }
    }
  }
  
  // Make coefficients and enter in coefficient DB
  //
  coefs->add(basis->makeFromArray(tnow));

  // Only root node writes the coefficient file
  //
  if (myid==0) {

    // Check if file exists and extend the existing HDF5 file
    //
    if (std::filesystem::exists(outfile)) {
      if (dof==2)
	std::dynamic_pointer_cast<CoefClasses::CylFldCoefs>(coefs)->
	  ExtendH5Coefs(outfile);
      else
	std::dynamic_pointer_cast<CoefClasses::SphFldCoefs>(coefs)->
	  ExtendH5Coefs(outfile);
    }
    // Otherwise, write a new HDF5 coefficient file
    //
    else {
      if (dof==2)
	std::dynamic_pointer_cast<CoefClasses::CylFldCoefs>(coefs)->
	  WriteH5Coefs(outfile);
      else
	std::dynamic_pointer_cast<CoefClasses::SphFldCoefs>(coefs)->
	  WriteH5Coefs(outfile);
    }
  }
  // END: write to file
}

