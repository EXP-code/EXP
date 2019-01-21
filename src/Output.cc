#include <boost/lexical_cast.hpp>

#include "expand.h"
#include <Output.H>


Output::Output(const YAML::Node& CONF) : conf(CONF)
{
  nint = 50;			// Default interval
}


void Output::bomb(const string& msg)
{
  cerr << "Output [" << id << ": " << msg << endl;
  MPI_Abort(MPI_COMM_WORLD, 499);
}
