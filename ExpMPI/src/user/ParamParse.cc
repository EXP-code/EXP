#include <iostream>
#include <iomanip>
#include <strstream>
#include <algorithm>

#include <ParamParse.H>
#include <StringTok.H>

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
   string::size_type where = value.find_last_not_of(' ');
  
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

ParamParse::ParamParse(istream* in, string Delim) : 
  curitem(NULL), enditem(NULL), curstanza(NULL), delim(Delim)
{
  if (in==NULL) return;

  const int lbufsize = 512;
  char lbuf[lbufsize];

  Stanza current;
  spair cpair;

				// Default name for first stanza
  current.name = "global";

  nstanza = 0;

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
	ostrstream sout;
	sout << "Syntax error: " << line;
	bomb(sout.str());
      }

				// Close the current stanza
      database.push_back(current);
      nstanza++;

				// Open the new stanza
      string name = line.substr( beg+1, end-beg-1 );
      current.name = name;
      current.elist.erase(current.elist.begin(), current.elist.end());

				// Is it a data line
    } else if ( (beg=line.find_first_of(delim)) != string::npos) {

      StringTok<string> tokens(line);
      cpair.first = trimLeft(trimRight(tokens(delim)));
      cpair.second = trimLeft(trimRight(tokens(delim)));
      
      current.elist.push_back(cpair);
    }
				// Forget about other lines
  }

				// Close the current stanza
  database.push_back(current);
  nstanza++;

}

void ParamParse::find_list(const string& stanza)
{
  list<Stanza>::iterator it;
  for (it=database.begin(); it!=database.end(); it++) {
    if (it->name.compare(stanza) == 0) {
      curstanza = it;
      curitem = it->elist.begin();
      enditem = it->elist.end();
      return;
    }
  }
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

int ParamParse::get_next(spair& ret)
{
  if (curitem != enditem) {
    ret = *(curitem++);
    return 1;
  }
  return 0;
}

void ParamParse::bomb(const char *msg)
{
  cerr << "ParamParse error: " << msg << "\n";
  exit(-1);
}

void ParamParse::print_database(ostream& out)
{
  int istanza = 0, iparam;
  list<Stanza>::iterator it;
  list<spair>::iterator it1;

  out << setiosflags(ios::left);

  for (it=database.begin(); it!=database.end(); it++) {
    out << setw(2) << istanza << " [" << it->name << "]" << endl;

    iparam = 0;
    for (it1=it->elist.begin(); it1!=it->elist.end(); it1++) {
      out << " " << setw(3) << iparam
	  << setw(15) << it1->first.c_str() << " | "
	  << setw(60) << it1->second.c_str()
	  << endl;
      iparam++;
    }

    istanza++;
    out << endl;
    
  }

  out << resetiosflags(ios::left);

}




