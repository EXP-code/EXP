#include "expand.H"
#include "Output.H"


Output::Output(const YAML::Node& CONF) : conf(CONF)
{
  // Default interval
  nint = 50;
  nintsub = std::numeric_limits<int>::max();

  // Add keys
  for (YAML::const_iterator it=conf.begin(); it!=conf.end(); ++it) {
    current_keys.insert(it->first.as<std::string>());
  }

}
