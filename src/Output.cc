#include "expand.H"
#include <Output.H>


Output::Output(const YAML::Node& CONF) : conf(CONF)
{
  nint = 50;			// Default interval
  nintsub = std::numeric_limits<int>::max();
}
