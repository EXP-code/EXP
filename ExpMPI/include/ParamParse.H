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
  void clear();
};


class ParamParse
{
protected:
  list<Stanza*> database;
  int nstanza;

  list<spair>::iterator curitem, enditem;
  Stanza* curstanza;
  string delim;

  void parse_istream(istream* in);

  void bomb(const string&);

public:

  //! Construct database from keyvalue with delimiter Delim
  ParamParse(const string&);

  /** 
      Construct database from keyvalue with delimiter Delim 
      and read from stream
  */
  ParamParse(istream* in, string Delim);
  
  /** 
      Construct database from keyvalue with delimiter Delim 
      and read from file
  */
  ParamParse(const string& filename, string Delim);
  
  //! Destructor
  ~ParamParse();

  //! Add keyvalues to database from a file
  int parse_file(const string& filename);

  //! Add keyvalues to database from argc,argv
  int parse_argv(int argc, char** argv);

  //! Set current to stanza name (default "global")
  bool find_list(const string& stanza);

  //! Find an item in the current stanza
  int find_item(const string& name, string& value);

  //! Get an item from the current stanza
  int get_next(spair& ret);

  //! Print the database to ostream
  void print_database(ostream& out);

  //! Add an item to current stanza
  void add_item(const string& name, const string& value);
};

#endif
