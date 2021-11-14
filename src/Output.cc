#include "expand.H"
#include <Output.H>


Output::Output(const YAML::Node& CONF) : conf(CONF)
{
  nint = 50;			// Default interval
  nintsub = std::numeric_limits<int>::max();
}


void Output::bomb(const string& msg)
{
  std::cerr << "Output [" << id << ": " << msg << endl;
  MPI_Abort(MPI_COMM_WORLD, 499);
}
