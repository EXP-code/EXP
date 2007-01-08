#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <typeinfo>

#include <time.h>

#include <ParamDatabase.H>
#include <StringTok.H>

void PP1_record::parse(const PP1_record *p, const string& s)
{
  istringstream istr(s);

  if (type.compare("bool")==0) {
    bool *i = new bool;
    if (p) *i = *(static_cast<bool *> (p->value));
    else if (s.size()) {
      if (s[0]=='t' || s[0]=='T') *i = true;
      else if (s[0]=='f' || s[0]=='F') *i = false;
      else *i = atoi(s.c_str()) ? true : false;
    }
    else throw "Illegal parse <bool>";
    value = i;
  }
  else if (type.compare("int")==0) {
    int *i = new int;
    if (p) *i = *(static_cast <int *> (p->value));
    else if (s.size()) istr >>*i;
    else throw "Illegal parse <int>";
    value = i;
  }
  else if (type.compare("float")==0) {
    float *i = new float;
    if (p) *i = *(static_cast <float *> (p->value));
    else if (s.size()) istr >> *i;
    else throw "Illegal parse <float>";
    value = i;
  }
  else if (type.compare("double")==0) {
    double *i = new double;
    if (p) *i = *(static_cast <double *> (p->value));
    else if (s.size()) istr >> *i;
    else throw "Illegal parse <float>";
    value = i;
  }
  else if (type.compare("string")==0) {
    string *i = new string;
    if (p) *i = *(static_cast <string *> (p->value));
    else if (s.size()) *i = s;
    // no "throw": allow null strings
    value = i;
  }
  else if (type.compare("char")==0) {
    char *i;
    if (p) *i = *(static_cast <unsigned char *> (p->value));
    else if (s.size()) {
      i = new char [s.size()+1];
      strncpy(i, s.c_str(), s.size()+1);
    }
    else throw "Illegal parse <char>";
    value = i;
  }
}

ParamDatabase::ParamDatabase(database_record *init) : 
  ParamParse("=")
{
  if (init) {
    int indx = 0;
    while (init[indx].name[0] != '\0') {
      add_entry(init[indx].name, init[indx].type, init[indx].dflt);
      indx++;
    }
  }
}

void ParamDatabase::add_entry(const string &name, const string& type, 
			      const string &deflt)
{
				// Look for duplicate entry
  list<PP1_record>::iterator it;
  for (it=database2.begin(); it!=database2.end(); it++) 
    {
      if (it->name.compare(name)==0) 
	{
	  it->type = type;
	  it->parse(0, deflt);
	  return;
	}
    }
	
				// Need new record
  PP1_record pp2;
  pp2.name = name;
  pp2.type = type;
  try {
    pp2.parse(0, deflt);
  }
  catch (char *msg) {
    cerr << msg << endl;
  }

  database2.push_back(pp2);
}


void ParamDatabase::set_entry(const string &name, const string& value)
{
  list<PP1_record>::iterator it;
  for (it=database2.begin(); it!=database2.end(); it++) 
    {
      if (it->name.compare(name)==0) { // Found it!

	istringstream ins(value.c_str());

	if (it->type.compare("bool")==0) {
	  bool *i = new bool;
	  if (value[0]=='t' || value[0]=='T') *i = true;
	  else if (value[0]=='f' || value[0]=='F') *i = false;
	  else *i = atoi(value.c_str()) ? true : false;
	  it->value = (void *)i;
	}
	else if (it->type.compare("int")==0) {
	  int *i = new int;
	  ins >> *i;
	  it->value = (void *)i;
	}
	else if (it->type.compare("float")==0) {
	  float *i = new float;
	  ins >> *i;
	  it->value = (void *)i;
	}
	else if (it->type.compare("double")==0) {
	  double *i = new double;
	  ins >> *i;
	  it->value = (void *)i;
	}
	else if (it->type.compare("string")==0) {
	  string *i = new string(value);
	  it->value = (void *)i;
	}
	else if (it->type.compare("char")==0) {
	  char *i = new char [value.size()+1];
	  strncpy(i, value.c_str(), value.size()+1);
	  it->value = (void *)i;
	}
	else {
	  ostringstream msg;
	  msg << "ParamDatabase::set_entry: no type <"
	      << it->type << ">";
	  cerr << msg << endl;
	  throw msg.str().c_str();
	}

	return;
      }

    }

  ostringstream msg;
  msg << "ParamDatabase::set_entry: no database entry for variable <"
      << name << ">";
  throw msg.str().c_str();
}


int ParamDatabase::parseFile(const string &file)
{
  parse_file(file);

  spair ret;
  while (get_next(ret)) 
    set_entry(ret.first, ret.second);
  
  return 0;
}

