#include "YamlCheck.H"

std::set<std::string> YamlCheck(const YAML::Node& node, const std::set<std::string>& valid)
{
  std::set<std::string> umatch;
  for (YAML::const_iterator it=node.begin(); it!=node.end(); ++it) {
    auto st = it->first.as<std::string>();
    if (valid.find(st) == valid.end()) umatch.insert(st);
  }
  return umatch;
}
