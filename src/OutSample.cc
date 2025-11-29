#include <iostream>
#include <fstream>
#include <iomanip>

#include "expand.H"
#include "Basis.H"
#include "OutSample.H"

const std::set<std::string>
OutSample::valid_keys = {
  "filename",
  "nint",
  "name"
};

OutSample::OutSample(const YAML::Node& conf) : Output(conf)
{
  nint      = 10;
  nintsub   = std::numeric_limits<int>::max();
  tcomp     = NULL;

  initialize();

  if (!(tcomp->force->hasSubsample())) {
    throw std::runtime_error("OutSample: no subsampling for this force");
  }

}

void OutSample::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["nint"])         nint     = conf["nint"].as<int>();
    if (conf["name"])
      {				// Search for desired component
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

    if (conf["filename"])
      {
	filename = outdir + conf["filename"].as<std::string>();
      }
    else
      {
	filename = outdir + "outcoef." + tcomp->name + "." + runtag;
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

  // Skip this master step
  //
  if (n % nint != 0 && !last) return;

  // Check for repeat time
  //
  if (tnow <= prev) return;

  prev = tnow;

  // Write the subsample data
  tcomp->force->writeSubsample();
}
