#ifndef _YAMLCHECK_H
#define _YAMLCHECK_H

#include <vector>
#include <string>
#include <set>

#include <yaml-cpp/yaml.h>

/**
   Check for unmatched configuration keys in a YAML::node

   Returns an std::vector of unmatched keys which is empty if all keys are matched
 */
std::vector<std::string> YamlCheck(const YAML::Node& node, const std::set<std::string>& valid);

#endif
