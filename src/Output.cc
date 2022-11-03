#include "expand.H"
#include <Output.H>


Output::Output(const YAML::Node& CONF) : conf(CONF)
{
  nint = 50;			// Default interval
  nintsub = std::numeric_limits<int>::max();
}


void Output::bomb(const string& msg)
{
  std::ostringstream sout;
  sout << "Output [" << id << ": " << msg;
  throw GenericError(sout.str(), __FILE__, __LINE__);
}
