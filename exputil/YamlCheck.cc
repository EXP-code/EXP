#include <YamlCheck.H>

std::vector<std::string> YamlCheck(const YAML::Node& node, const std::set<std::string>& valid)
{
  std::vector<std::string> umatch;
  for (YAML::const_iterator it=node.begin(); it!=node.end(); ++it) {
    auto st = it->first.as<std::string>();
    if (valid.find(st) == valid.end()) umatch.push_back(st);
  }
  return umatch;
}
