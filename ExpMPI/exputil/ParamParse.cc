#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace std;

#include <ParamParse.H>
#include <StringTok.H>

Stanza::Stanza(const Stanza &p)
{
  name = p.name;
  elist = p.elist;
}

void Stanza::clear()
{
#if __GNUC__ && (__GNUC__ < 3)
  name.erase(name.begin(), name.end());
#else
  name.clear();
#endif
  elist.clear();
}

string trimLeft(const string value)
{
   // stripping only space (character #32)
   string::size_type where = value.find_first_not_of(' ');
  
   if (where == string::npos)
     // string has nothing but space
     return string();
   
   if (where == 0)
     // string has no leading space, don't copy its contents
     return value;
  
   return value.substr(where);
 }


string trimRight(const string value)
 {
   string::size_type where = value.find_last_not_of(" \t\n");
  
   if (where == string::npos)
     // string has nothing but space
     return string();
   
   if (where == (value.length() - 1))
     // string has no trailing space, don't copy its contents
     return value;
   
   return value.substr(0, where + 1);
 }


string trimComment(const string value)
 {
   string::size_type where = value.find_first_of("#!");
  
   if (where == string::npos)
     // no comments
     return value;
   
   return value.substr(0, max<int>(where - 1, 0));
 }

void ParamParse::parse_istream(istream* in)
{
  if (in==NULL) return;

  const int lbufsize = 512;
  char lbuf[lbufsize];

  spair cpair;
  string theEnd("\n");
				// Default name for first stanza
  if (database.size() == 0) {
    curstanza = new Stanza();
    database.push_back(curstanza);
    curstanza->name = "global";
    nstanza = 1;
  } else {
    if (!find_list("global")) 
      bomb("Couldn't find global stanza; this should never happen");
  }
  
  while (*in) {
				// Read in the line
    in->getline(lbuf, lbufsize);
    if (!*in) break;
    string::size_type beg;

    string line(lbuf);

				// Strip comments, 
				// leading and trailing whitespace
    line = trimComment(line);
    line = trimLeft(trimRight(line));


				// Skip empty line
    if (line.empty()) continue;

				// Look for beginning of new stanza
    if ( (beg=line.find("[")) != string::npos) {

				// Must be followed "]"
      string::size_type end = line.find("]");
      if (end == string::npos) {
	ostringstream sout;
	sout << "Syntax error: " << line;
	bomb(sout.str());
      }

				// Try to open the stanza
      string name = line.substr( beg+1, end-beg-1 );

				// Add new stanza needed
      if (!find_list(name)) {
	curstanza = new Stanza();
	database.push_back(curstanza);
	curstanza->name = name;
	nstanza++;
      }
      
				// Is it a data line
    } else if ( (beg=line.find(delim)) != string::npos) {

      StringTok<string> tokens(line);
      cpair.first = trimLeft(trimRight(tokens(delim)));
      cpair.second = trimLeft(trimRight(tokens(theEnd)));
      
      curstanza->elist.push_back(cpair);
    }
				// Forget about other lines
  }

				// Make the current stanza the first
  if (!find_list("global"))
    bomb ("Could not reset to global stanza");

}

ParamParse::ParamParse(istream* in, string Delim) : 
  curitem(NULL), enditem(NULL), curstanza(NULL), delim(Delim)
{
  parse_istream(in);
}

ParamParse::ParamParse(const string& filename, string Delim) : 
  curitem(NULL), enditem(NULL), curstanza(NULL), delim(Delim)
{
  parse_file(filename);
}

int ParamParse::parse_file(const string& filename)
{
  ifstream in(filename.c_str());
  parse_istream(&in);

  return curstanza->elist.size();
}

int ParamParse::parse_argv(int argc, char** argv)
{
  string data;
  for (int i=0; i<argc; i++) {
    if (i) data += "\n";
    data += argv[i];
  }
  istringstream in(data);
  parse_istream(&in);

  return curstanza->elist.size();
}


ParamParse::ParamParse(const string& Delim) : 
  curitem(NULL), enditem(NULL), curstanza(NULL)
{
  delim = Delim;
  curstanza = new Stanza();
  database.push_back(curstanza);
  curstanza->name = "global";
  nstanza = 1;
}

ParamParse::~ParamParse()
{
  list<Stanza*>::iterator it;
  for (it=database.begin(); it!=database.end(); it++) delete *it;
}


bool ParamParse::find_list(const string& stanza)
{
  list<Stanza*>::iterator it;
  for (it=database.begin(); it!=database.end(); it++) {
    if ((*it)->name.compare(stanza) == 0) {
      curstanza = *it;
      curitem = curstanza->elist.begin();
      enditem = curstanza->elist.end();
      return true;
    }
  }

  return false;
}


int ParamParse::find_item(const string& name, string& value)
{
  list<spair>::iterator p;
  for (p = curstanza->elist.begin(); p != curstanza->elist.end(); p++) {
    if (name.compare(p->first) == 0) {
      value = p->second;
      return 1;
    }
  }
  return 0;
}

void ParamParse::add_item(const string& name, const string& value)
{
  spair p;

  p.first = name;
  p.second = value;
  curstanza->elist.push_back(p);
}


int ParamParse::get_next(spair& ret)
{
  if (curitem != enditem) {
    ret = *(curitem++);
    return 1;
  }
  return 0;
}

void ParamParse::bomb(const string& msg)
{
  cerr << "ParamParse error: " << msg << "\n";
  exit(-1);
}

void ParamParse::print_database(ostream& out)
{
  int istanza = 0, iparam;
  list<Stanza*>::iterator it;
  list<spair>::iterator it1;

  out << setiosflags(ios::left);

  for (it=database.begin(); it!=database.end(); it++) {
    out << setw(2) << istanza << " [" << (*it)->name << "]" << endl;

    iparam = 0;
    for (it1=(*it)->elist.begin(); it1!=(*it)->elist.end(); it1++) {
      out << " " << setw(3) << iparam
	  << setw(15) << it1->first.c_str()  << " | "
	  << it1->second << endl;
      iparam++;
    }

    istanza++;
    out << endl;
    
  }

  out << resetiosflags(ios::left);
  
}




