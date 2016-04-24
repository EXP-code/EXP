#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <typeinfo>
#include <utility>
#include <cstring>

using namespace std;

#include <unistd.h>
#include <getopt.h>
#include <time.h>

#include <ProgramParam.H>
#include <StringTok.H>

inline string ltoa(bool logic)
{
  if (logic) return "True";
  else return "False";
}

void PP2_record::parse(const PP2_record *p, const string& s)
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
    char *i=0;
    if (p) *i = *(static_cast <unsigned char *> (p->value));
    else if (s.size()) {
      i = new char [s.size()+1];
      strncpy(i, s.c_str(), s.size()+1);
    }
    else throw "Illegal parse <char>";
    value = i;
  }
}

ProgramParam::ProgramParam(const char *msg, program_option *init) : 
  ParamParse("="), desc(msg)
{
  if (init) {
    int indx = 0;
    while (init[indx].name[0] != '\0') {
      add_entry(init[indx].name, init[indx].type, 
		init[indx].dflt, init[indx].desc);
      indx++;
    }
  }
}

void ProgramParam::add_entry(const string &name, const string& type, 
			     const string &deflt, const string &help)
{
				// Look for duplicate entry
  for (auto &j : database2)
    {
      if (j.name.compare(name)==0) 
	{
	  j.type = type;
	  j.help = help;
	  try {
	    j.parse(0, deflt);
	  }
	  catch (char *msg) {
	    cerr << "On parsing <" << name << ">" << endl;
	    cerr << msg << endl;
	  }
	  return;
	}
    }
	
				// Need new record
  PP2_record pp2;
  pp2.name = name;
  pp2.help = help;
  pp2.type = type;
  try {
    pp2.parse(0, deflt);
  }
  catch (char *msg) {
    cerr << "On parsing <" << name << "> [2]" << endl;
    cerr << msg << endl;
  }

  database2.push_back(pp2);
}


void ProgramParam::set_entry_list(vector<sspair>& params)
{
  for (auto v : params)
    set_entry(v.first, v.second);
}


void ProgramParam::set_entry(const string &name, const string& value)
{
  for (auto &j : database2)
    {
      if (j.name.compare(name)==0) { // Found it!

	istringstream ins(value.c_str());

	if (j.type.compare("bool")==0) {
	  bool *i = new bool;
	  if (value[0]=='t' || value[0]=='T') *i = true;
	  else if (value[0]=='f' || value[0]=='F') *i = false;
	  else *i = atoi(value.c_str()) ? true : false;
	  j.value = (void *)i;
	}
	else if (j.type.compare("int")==0) {
	  int *i = new int;
	  ins >> *i;
	  j.value = (void *)i;
	}
	else if (j.type.compare("float")==0) {
	  float *i = new float;
	  ins >> *i;
	  j.value = (void *)i;
	}
	else if (j.type.compare("double")==0) {
	  double *i = new double;
	  ins >> *i;
	  j.value = (void *)i;
	}
	else if (j.type.compare("string")==0) {
	  string *i = new string(value);
	  j.value = (void *)i;
	}
	else if (j.type.compare("char")==0) {
	  char *i = new char [value.size()+1];
	  strncpy(i, value.c_str(), value.size()+1);
	  j.value = (void *)i;
	}
	else {
	  ostringstream msg;
	  msg << "ProgramParam::set_entry: no type <"
	      << j.type << ">";
	  std::cerr << msg.str() << std::endl;
	  throw msg.str().c_str();
	}

	return;
      }

    }

  ostringstream msg;
  msg << "ProgramParam::set_entry: no database entry for variable <"
      << name << ">";
  throw msg.str().c_str();
}


int ProgramParam::parse_args(int argc, char **argv)
{
  char *prog=argv[0];
  int iparmf=0,iret;
  string file;
  int c;

  while (1) {
    c = getopt(argc, argv, "f:dh");
    if (c == -1) break;

    switch (c) {
    case 'f': 
      iparmf=1;
      file = optarg;
      break;
    case 'd':
      print_default(prog);
      exit(0);
      break;
    case '?':
      usage(prog);
      exit(0);
      break;
    case 'h':
      usage(prog);
      help();
      exit(0);
    }
  }

    argc -= optind;
    if (iparmf)
      parse_file(file);
    else
      {
	iret = parse_argv(argc, &argv[optind]);
	argc -= iret;
      }
  
    if (argc != 0)
      usage(prog);

				// Set parameters

    spair ret;
    while (get_next(ret)) 
      set_entry(ret.first, ret.second);

    return 0;
}


void ProgramParam::help()
{
  cerr << endl;
  cerr << setw(12) << left << "Keyword" 
       << setw(8)  << left << "type"
       << setw(16) << left << "default value"
       << setw(40) << left << "description"
       << endl
       << setw(12) << left << "-------" 
       << setw(8)  << left << "----"
       << setw(16) << left << "-------------"
       << setw(40) << left << "-----------"
       << endl;


  for (auto  j : database2)
    {
      cerr << setw(12) << j.name
	   << setw(8)  << j.type;

      if (j.type.compare("bool")==0)
	cerr << setw(16) << ltoa(*(static_cast<bool*>(j.value)));
      else if (j.type.compare("int")==0) 
	cerr << setw(16) << *(static_cast<int*>(j.value));
      else if (j.type.compare("float")==0)
	cerr << setw(16) << *(static_cast<float*>(j.value));
      else if (j.type.compare("double")==0)
	cerr << setw(16) << *(static_cast<double*>(j.value));
      else if (j.type.compare("string")==0)
	cerr << setw(16) << *(static_cast<string*>(j.value));
      else if (j.type.compare("char")==0)
	cerr << setw(16) << *(static_cast<char*>(j.value));
      else
	cerr << setw(16) << "unknown";
  
      cerr << j.help << endl;
    }
}

void ProgramParam::print_default(char *prog)
{
  time_t tm = time(0);

  cerr << "# Default values for <" << prog << ">" << endl 
       << "# Local time: " << ctime(&tm) << endl;

  for (auto j : database2)
    {
      cerr << left << setw(20) << j.name << " = ";

      if (j.type.compare("bool")==0)
	cerr << setw(20) << ltoa(*(static_cast<bool*>(j.value)));
      else if (j.type.compare("int")==0) 
	cerr << setw(20) << *(static_cast<int*>(j.value));
      else if (j.type.compare("float")==0)
	cerr << setw(20) << *(static_cast<float*>(j.value));
      else if (j.type.compare("double")==0)
	cerr << setw(20) << *(static_cast<double*>(j.value));
      else if (j.type.compare("string")==0)
	cerr << setw(20) << *(static_cast<string*>(j.value));
      else if (j.type.compare("char")==0)
	cerr << setw(20) << *(static_cast<char*>(j.value));
      else
	cerr << setw(20) << "unknown";
  
      cerr << "# " << j.help << endl;
    }
}

void ProgramParam::print_params(char *prog, string& file)
{
  time_t tm = time(0);

  ofstream out(file.c_str());
  if (out) {

    out << "# Values for <" << prog << ">" << endl 
	<< "# Run time: " << ctime(&tm) << endl;

    for (auto j : database2)
      {
	out << left << setw(20) << j.name << " = ";

	if (j.type.compare("bool")==0)
	  out << setw(20) << *(static_cast<bool*>(j.value));
	else if (j.type.compare("int")==0) 
	  out << setw(20) << *(static_cast<int*>(j.value));
	else if (j.type.compare("float")==0)
	  out << setw(20) << *(static_cast<float*>(j.value));
	else if (j.type.compare("double")==0)
	  out << setw(20) << *(static_cast<double*>(j.value));
	else if (j.type.compare("string")==0)
	  out << setw(20) << *(static_cast<string*>(j.value));
	else if (j.type.compare("char")==0)
	  out << setw(20) << *(static_cast<char*>(j.value));
	else
	  out << setw(20) << "unknown";
	
	out << "# " << j.help << endl;
      }
  } else {
    ostringstream msg;
    msg << "ProgramParam::print_params: can not open file <"
	<< file << ">";
    throw msg.str().c_str();
  }
}

void ProgramParam::usage(char *prog)
{
  char usage_head[] = 
    "[-f file -d] [keyword=value [keyword=value] .. ]";

  cerr << "Usage: " << prog << " " << usage_head << endl;
  cerr << setw(18) << left << "     -f file"
       << "keyword/value parameter file" << endl;
  cerr << setw(18) << left << "     -d"
       << "print default parameters" << endl;
  cerr << setw(18) << left << "     -h"
       << "print default parameters with help" << endl;
  cerr << endl;
  cerr << desc << endl;
}
