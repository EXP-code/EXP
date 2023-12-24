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
  lmax    = 4;
  nmax    = 10;
  rmin    = 1.0e-4;
  rmax    = 2.0;
  ascl    = 0.01;
  delta   = 0.005;
  scale   = 0.05;
  model   = "file";

  initialize();

  if (!(tcomp->force->HaveCoefDump())) {
    if (myid==0) {
      cerr << "OutVel: no coefficients for this force\n";
    }
  }

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

  // Allocate storage for coefficients
  //
  coefs.resize(dof);

  if (dof==2) {
    for (auto & v : coefs) v.resize(lmax+1, nmax);
  } else {
    for (auto & v : coefs) v.resize((lmax+1)*(lmax+2)/2, nmax);
  }

  // Create model needed for density prefactor in OrthoFunction
  //
  if (model == "file") {
    std::vector<double> r, d;
    std::ifstream in(filename);
    if (not in) throw std::runtime_error("Error opening file: " + filename);

    std::string line;
    while (std::getline(in, line)) {
      auto pos = line.find_first_of("!#");
      if (pos == std::string::npos) {
	std::istringstream iss(line);
	double x, y;
	iss >> x >> y;
	if (iss) {
	  r.push_back(x);
	  d.push_back(y);
	}
      }
    }
    
    
    double rmin = r.front(), rmax = r.back();
  
    interp = std::make_shared<Linear1d>(r, d);
    densfunc = [this](double r)
    {
      return this->interp->eval(r);
    };
    
  } else if (model == "expon") {

    densfunc = [&](double r)
    {
      return exp(-r/ascl) * 0.5*(1.0 + std::erf((rmax - 5.0*delta - r)/delta)) / ascl;
    };

  } else {
    throw InternalError("OutVel: model logic failure?! "
			"You should not be here...",
			__FILE__, __LINE__);
  }

  // Generate the orthogonal function instance
  //
  ortho = std::make_shared<OrthoFunction>(nmax, densfunc, rmin, rmax, scale);

}

void OutVel::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["nint"])         nint     = conf["nint" ].as<int>();
    if (conf["nintsub"]) {
      nintsub  = conf["nintsub"].as<int>();
      if (nintsub <= 0) nintsub = 1;
    }

    if (conf["lmax"])         lmax     = conf["lmax" ].as<int>();
    if (conf["nmax"])         nmax     = conf["nmax" ].as<int>();
    if (conf["rmin"])         rmin     = conf["rmin" ].as<double>();
    if (conf["rmax"])         rmax     = conf["rmax" ].as<double>();
    if (conf["ascl"])         ascl     = conf["ascl" ].as<double>();
    if (conf["delta"])        delta    = conf["delta"].as<double>();
    if (conf["scale"])        scale    = conf["scale"].as<double>();
    if (conf["model"])        model    = conf["model"].as<std::string>();

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
	filename = outdir + "outcoef." + tcomp->name + "." + runtag;
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

  if (myid==0) {
    // Notes.  To be implemented:
    // 1. Zero coefficients
    // 2. Run through phase space and compute expansion of flow field
    // 3. Save and append to HDF5 file
  }

}
