#include <limits>

#include "expand.H"

#include <EJcom.H>

const std::set<std::string>
EJcom::valid_keys = {
  "cfac",
  "alpha"
};

EJcom::EJcom(Component* c0, const YAML::Node& conf) : TwoCenter(c0, conf)
{
  id = "EJcom";

  // Defaults
  //
  cfac  = 1.0;
  alpha = 1.0;

  // Get initialization info
  //
  initialize();
}

void EJcom::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("TwoCenter", "parameter", unmatched, __FILE__, __LINE__, 1020);

  // Assign values from YAML
  //
  if (conf["cfac"])  cfac  = conf["cfac"].as<double>();
  if (conf["alpha"]) alpha = conf["alpha"].as<double>();
}


double EJcom::mixture(double* pos)
{
  double del=0.0, dif=0.0;

  for (int k=0; k<3; k++) {
    del += (pos[k]   - inner[k]) * (pos[k]   - inner[k]);
    dif += (outer[k] - inner[k]) * (outer[k] - inner[k]);
  }

  double value = erf(cfac*pow(del/(dif+1.0e-10), 0.5*alpha));
  
  if (multistep==0 || mstep==0) accum_histo(value);

  return value;
}
