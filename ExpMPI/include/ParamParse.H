#ifndef _ParamParse_H
#define _ParamParse_H

#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <utility>

using namespace std;

typedef pair<string, string> spair;

class Stanza {
public:
  string name;
  list<spair> elist;

  Stanza() {};
  Stanza(const Stanza &);
};


class ParamParse
{
protected:
  list<Stanza> database;
  int nstanza;

  list<spair>::iterator curitem, enditem;
  list<Stanza>::iterator curstanza;
  string delim;

  void bomb(const string&);

public:

  ParamParse(istream* in, string Delim);
  
  void find_list(const string& stanza);
  int find_item(const string& name, string& value);
  int get_next(spair& ret);
  void print_database(ostream& out);
};

#endif
