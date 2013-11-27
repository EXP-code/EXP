#ifndef _ParamDatabase_H
#define _ParamDatabase_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>
#include <typeinfo>
#include <utility>

#include <ParamParse.H>

using namespace std;

class PP1_record 
{
public:
  string name, type;
  void *value;

  void parse(const PP1_record *p, const string& s);

  PP1_record() { value = 0;};
  PP1_record(const PP1_record &s)
  {
    name = s.name;
    type = s.type;
    parse(&s, "");
  }

  // ~PP1_record() {}
};

struct database_record
{
  char name[20];
  char type[20];
  char dflt[20];
};

class ParamDatabase : public ParamParse
{
protected:
  list<PP1_record> database2;

public:

  //! Constructors
  ParamDatabase() : ParamParse("=") {}
  ParamDatabase(database_record *init);

  //! Add an entry
  void add_entry(const string &name, const string& type, const string& deflt);

  //! Set an entry
  void set_entry(const string &name, const string &value);

  //! Get an entry
  template <typename T> 
  T get(const string &name);

  //! Set a value
  template <typename T> 
  void set(const string &name, const T& value);

  //! Parse file
  int parseFile(const string &filename);

};


template<typename T>
void ParamDatabase::set(const string& name, const T& val)
{
  list<PP1_record>::iterator it;
  for (it=database2.begin(); it!=database2.end(); it++) 
    {
      if (it->name.compare(name)==0) { // Found it!

	if (typeid(T)==typeid(bool)) {
	  if (it->type.compare("bool")) {
	    ostringstream msg;
	    msg << "ParamDatabase::set: expected <bool> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(int)) {
	  if (it->type.compare("int")) {
	    ostringstream msg;
	    msg << "ParamDatabase::set: expected <int> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(float)) {
	  if (it->type.compare("float")) {
	    ostringstream msg;
	    msg << "ParamDatabase::set: expected <float> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(double)) {
	  if (it->type.compare("double")) {
	    ostringstream msg;
	    msg << "ParamDatabase::set: expected <double> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(string)) {
	  if (it->type.compare("string")) {
	    ostringstream msg;
	    msg << "ParamDatabase::set: expected <string> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else if (typeid(T)==typeid(char)) {
	  if (it->type.compare("char")) {
	    ostringstream msg;
	    msg << "ParamDatabase::set: expected <char> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	}

	else {

	  ostringstream msg;
	  msg << "ParamDatabase::set: type mismatch on <"
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
  msg << "ParamDatabase::set: no database entry for variable <"
      << name << ">";
  throw msg.str().c_str();
}



template<typename T>
T ParamDatabase::get(const string &name)
{
  list<PP1_record>::iterator it;
  for (it=database2.begin(); it!=database2.end(); it++) 
    {
      if (it->name.compare(name)==0) { // Found it!

	if (typeid(T)==typeid(bool)) {
	  if (it->type.compare("bool")) {
	    ostringstream msg;
	    msg << "ParamDatabase::get: expected <bool> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(int)) {
	  if (it->type.compare("int")) {
	    ostringstream msg;
	    msg << "ParamDatabase::get: expected <int> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(float)) {
	  if (it->type.compare("float")) {
	    ostringstream msg;
	    msg << "ParamDatabase::get: expected <float> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(double)) {
	  if (it->type.compare("double")) {
	    ostringstream msg;
	    msg << "ParamDatabase::get: expected <double> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(string)) {
	  if (it->type.compare("string")) {
	    ostringstream msg;
	    msg << "ParamDatabase::get: expected <string> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

	else if (typeid(T)==typeid(char)) {
	  if (it->type.compare("char")) {
	    ostringstream msg;
	    msg << "ParamDatabase::get: expected <char> type and we have <"
		<< it->type << "> type!";
	    throw msg.str().c_str();
	  }
	  return *(static_cast<T *>(it->value));
	}

      }

    }

  ostringstream msg;
  msg << "ParamDatabase::get: no database entry for variable <"
      << name << ">";
  cerr << msg.str() << endl;
  throw msg.str().c_str();

  return (T)(0);
}

#endif
