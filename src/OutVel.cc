#include <filesystem>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <expand.H>

#include <OutVel.H>

const std::set<std::string>
OutVel::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "name",
  "dof",
  "scale",
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
  nint    = 10;
  nintsub = std::numeric_limits<int>::max();
  tcomp   = NULL;
  dof     = 3;

  initialize();

  // Target output file
  //
  outfile = "velcoef." + tcomp->name;

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
  if (dof!=2 or dof!=3) {
    std::ostringstream sout;
    sout << "OutVel: found " << dof << " for dof.  Must be 2 or 3.";
    throw std::runtime_error(sout.str());
  }

  // Create the basis
  //
  basis = std::make_shared<BasisClasses::VelocityBasis>(conf);

  // Create the coefficient container
  //
  if (dof==2)
    coefs = std::make_shared<CoefClasses::PolarVelCoefs>();
  else
    coefs = std::make_shared<CoefClasses::SphVelCoefs>();
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
    if (conf["name"])
      {				// Search for desired component
	std::string tmp = conf["name"].as<std::string>();
	for (auto c : comp->components) {
	  if (!(c->name.compare(tmp))) tcomp  = c;
	}
      }

    if (!tcomp) {
      std::string message = "OutVel: no component to trace. Please specify "
	"the component name using the 'name' parameter.";
      throw std::runtime_error(message);
    }

    if (conf["filename"])
      {
	filename = conf["filename"].as<std::string>();
      }
    else
      {
	filename = outdir + "velcoef." + tcomp->name + "." + runtag;
      }

  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutVel: "
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

void OutVel::Run(int n, int mstep, bool last)
{
  // Skip this master step
  //
  if (n % nint != 0 && !last) return;

  // Skip this sub step
  //
  if (mstep < std::numeric_limits<int>::max() and mstep % nintsub != 0) return;

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
  
  // Make coefficients
  //
  basis->make_coefs();

  // Check if file exists and extend the existing HDF5 file
  //
  if (std::filesystem::exists(outfile)) {
    if (dof==2)
      std::dynamic_pointer_cast<CoefClasses::PolarVelCoefs>(coefs)->
	ExtendH5Coefs(outfile);
    else
      std::dynamic_pointer_cast<CoefClasses::SphVelCoefs>(coefs)->
	ExtendH5Coefs(outfile);
  }
  // Otherwise, write a new HDF5 coefficient file
  //
  else {
    if (dof==2)
      std::dynamic_pointer_cast<CoefClasses::PolarVelCoefs>(coefs)->
	WriteH5Coefs(outfile);
    else
      std::dynamic_pointer_cast<CoefClasses::SphVelCoefs>(coefs)->
	WriteH5Coefs(outfile);
  }
}

