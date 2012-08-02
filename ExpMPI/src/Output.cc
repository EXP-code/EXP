#include <boost/lexical_cast.hpp>

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


std::map<int, std::string> Output::get_value_array(const string& name)
{
  std::map<int, string> values;
  int indx;

  list< pair<string, string> >::iterator it;
  for (it=namevalue.begin(); it!=namevalue.end(); it++) {
    string key = name + "(";
    if (it->first.compare(0, key.size(), key) == 0) {
      string sindx = it->first.substr(key.size(), it->first.find(")"));
      try {
	indx = boost::lexical_cast<int>(sindx);
      } 
      catch( boost::bad_lexical_cast const& ) {
	std::cout << "Output::get_value_array: input string <"
		  << it->first << "> is not valid" << std::endl;
      }
      values[indx] = it->second;
    }
  }
  return values;
}

void Output::bomb(const string& msg)
{
  cerr << "Output [" << id << ": " << msg << endl;
  MPI_Abort(MPI_COMM_WORLD, 499);
}
