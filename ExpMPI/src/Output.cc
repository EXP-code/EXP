
#include "expand.h"
#include <Output.H>


Output::Output(string& line)
{
  StringTok<string> tokens(line);

				// Default interval
  nint = 50;

				// comma separated tokens
  string token = tokens(",");
  
  while (token.size()) {

    StringTok<string> parse(token);
    pair<string, string> spair;
    spair.first = trimLeft(trimRight(parse("=")));
    spair.second = trimLeft(trimRight(parse("=")));
    namevalue.push_back(spair);

				// Next parameter
    token = tokens(",");
  }

}

int Output::get_value(const string& name, string& value)
{
  list< pair<string, string> >::iterator it;
  for (it=namevalue.begin(); it!=namevalue.end(); it++) {
    if (it->first.compare(name) == 0) {
      value = it->second;
      return 1;
    }
  }
  return 0;
}


void Output::bomb(const string& msg)
{
  cerr << "Output [" << id << ": " << msg << endl;
  MPI_Abort(MPI_COMM_WORLD, 499);
}
