#include <iostream>
#include <fstream>
#include <iomanip>

#include "expand.H"
#include "Basis.H"
#include "OutSample.H"

const std::set<std::string>
OutSample::valid_keys = {
  "floatType",
  "filename",
  "name",
  "level",
  "chunksize",
  "compress",
  "szip"
};

OutSample::OutSample(const YAML::Node& conf) : Output(conf)
{
  nintsub   = std::numeric_limits<int>::max();
  tcomp     = NULL;

  initialize();

  if (!(tcomp->force->hasSubsample())) {
    throw std::runtime_error("OutSample: no subsampling for this force");
  }

  covarStore = std::make_shared<BasisClasses::SubsampleCovariance>
    ([this](HighFive::File& file){this->tcomp->force->writeCovarH5Params(file);},
     this->tcomp->force->id, floatType, this->tcomp->force->FullCovar());

  covarStore->setCovarH5Compress(level, chunksize, shuffle, szip);

  if (myid==0)
    std::cout << "---- OutSample: writing subsample coefficients for component '"
	      << tcomp->name << "' to file " << filename << " for force '"
	      << tcomp->force->id << "'" << std::endl;
}

void OutSample::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    // Search for desired component by name
    if (conf["name"]) {
      std::string tmp = conf["name"].as<std::string>();
      for (auto c : comp->components) {
	if (!(c->name.compare(tmp))) tcomp  = c;
      }
    }

    if (!tcomp) {
      std::string message = "OutSample: no component to trace. Please specify "
	"the component name using the 'name' parameter.";
      throw std::runtime_error(message);
    }

    if (conf["filename"]) {
      filename = outdir + conf["filename"].as<std::string>();
    } else {
      filename = outdir + "outcoef." + tcomp->name + "." + runtag;
    }

    if (conf["floatType"]) {
      floatType = conf["floatType"].as<bool>();
    }

    if (conf["level"]) {
      level = conf["level"].as<unsigned>();
    }

    if (conf["chunksize"]) {
      chunksize = conf["chunksize"].as<unsigned>();
    }

    if (conf["shuffle"]) {
      shuffle = conf["shuffle"].as<bool>();
    }

    if (conf["szip"]) {
      szip = conf["szip"].as<bool>();
    }

  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutSample: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("OutSample::initialize: error parsing YAML");
  }
}

void OutSample::Run(int n, int mstep, bool last)
{
  // Subsampling is not available for output
  //
  if (not tcomp->force->subsampleReady()) return;
  
  // Check for repeat time
  //
  if (tnow <= prev) return;

  // Cache the current simulation time
  prev = tnow;
    
  // Write the subsample data
  auto elem = tcomp->force->getSubsample();
  covarStore->writeCoefCovariance(tcomp->name, runtag, elem, tnow);

}
