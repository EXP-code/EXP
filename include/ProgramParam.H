#ifndef _ProgramParam_H
#define _ProgramParam_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include <ParamParse.H>

using namespace std;

class PP2_record 
{

public:
  string name, help, type;
  void *value;

  void parse(const PP2_record *p, const string& s);

  PP2_record() : value(0) {};
  PP2_record(const PP2_record &s)
  {
    name = s.name;
    help = s.help;
    type = s.type;
    parse(&s, "");
  }
  // ~PP2_record() { delete value; }
};

struct program_option 
{
  char name[20];
  char type[20];
  char dflt[40];
  char desc[160];
};

typedef pair<string, string> sspair;

class ProgramParam : public ParamParse
{
protected:
  list<PP2_record> database2;
  const char *desc;

public:

  //! Constructors
  ProgramParam() : ParamParse("="), desc("") {}
  ProgramParam(const char *desc, program_option *init=0);

  //! Usage
  void usage(char *prog);

  //! Print default values
  void print_default(char *prog);

  //! Print current values
  void print_params(char *prog, string& file);

  //! Help on parameters
  void help();

  //! Add an entry
  void add_entry(const string &name, const string& type, const string& deflt, const string &help);

  //! Set an entry list
  void set_entry_list(vector<sspair> &params);

  //! Set an entry
  void set_entry(const string &name, const string &value);

  //! Get an entry
  template <typename T> 
  T get(const string &name);

  //! Set a value
  template <typename T> 
  void set(const string &name, const T& value);

  //! Set the description string
  void set_description(char* text) { desc = text; }

  //! Parse command line
  int parse_args(int argc, char **argv);

};


template<typename T>
void ProgramParam::set(const string& name, const T& val)
{
  list<PP2_record>::iterator it;
  for (it=database2.begin(); it!=database2.end(); it++) 
    {
      if (it->name.compare(name)==0) { // Found it!

	if (typeid(T)==typeid(bool)) {
	  if (it->type.compare("bool")) {
	    ostringstream msg;
	    msg << "ProgramParam::set: expected <bool> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(int)) {
	  if (it->type.compare("int")) {
	    ostringstream msg;
	    msg << "ProgramParam::set: expected <int> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(float)) {
	  if (it->type.compare("float")) {
	    ostringstream msg;
	    msg << "ProgramParam::set: expected <float> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(double)) {
	  if (it->type.compare("double")) {
	    ostringstream msg;
	    msg << "ProgramParam::set: expected <double> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(string)) {
	  if (it->type.compare("string")) {
	    ostringstream msg;
	    msg << "ProgramParam::set: expected <string> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(char)) {
	  if (it->type.compare("char")) {
	    ostringstream msg;
	    msg << "ProgramParam::set: expected <char> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else {

	  ostringstream msg;
	  msg << "ProgramParam::set: type mismatch on <"
	      << it->type << "> type!";
	  throw msg.str().c_str();
	}

	T* p = new T;
	*p = val;
	it->value = (void *)p;

	return;
      }

    }

  ostringstream msg;
  msg << "ProgramParam::set: no database entry for variable <"
      << name << ">";
  throw msg.str().c_str();
}



template<typename T>
T ProgramParam::get(const string &name)
{
  list<PP2_record>::iterator it;
  for (it=database2.begin(); it!=database2.end(); it++) 
    {
      if (it->name.compare(name)==0) { // Found it!

	if (typeid(T)==typeid(bool)) {
	  if (it->type.compare("bool")) {
	    ostringstream msg;
	    msg << "ProgramParam::get: expected <bool> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(int)) {
	  if (it->type.compare("int")) {
	    ostringstream msg;
	    msg << "ProgramParam::get: expected <int> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(float)) {
	  if (it->type.compare("float")) {
	    ostringstream msg;
	    msg << "ProgramParam::get: expected <float> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(double)) {
	  if (it->type.compare("double")) {
	    ostringstream msg;
	    msg << "ProgramParam::get: expected <double> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(string)) {
	  if (it->type.compare("string")) {
	    ostringstream msg;
	    msg << "ProgramParam::get: expected <string> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(char)) {
	  if (it->type.compare("char")) {
	    ostringstream msg;
	    msg << "ProgramParam::get: expected <char> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

      }

    }

  ostringstream msg;
  msg << "ProgramParam::get: no database entry for variable <"
      << name << ">";
  cerr << msg.str() << endl;
  throw msg.str().c_str();

  return (T)(0);
}

#endif
