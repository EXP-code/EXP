#include "parseVersionString.H"

std::vector<int> parseVersionString(const std::string& version_str)
{
  std::vector<int> version;
  size_t pos = 0;

  while (pos<version_str.size()) {

    size_t dot_pos = version_str.find('.', pos);
    size_t n = dot_pos;
    if (dot_pos != std::string::npos) n = dot_pos - pos;

    version.push_back(atoi(version_str.substr(pos, n).c_str()));

    if (dot_pos == std::string::npos) return version;

    pos = dot_pos + 1;
  }

  return version;
}
