#ifndef _ParamParse_H
#define _ParamParse_H

#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <pair.h>

typedef pair<string, string> spair;

//! Structure to keep a key/value list for a stanza in the parameter file
struct Stanza {
  //! Stanza name
  string name;
  //! List of key/value pairs
  list<spair> elist;
};


//! Caches key/value pairs
class ParamParse
{
protected:
  //! List of databases for each stanza
  list<Stanza> database;

  //! Number of stanzas in database
  int nstanza;

  //! Cache current and final items
  list<spair>::iterator curitem, enditem;

  //! Cache current stanza
  list<Stanza>::iterator curstanza;

  //! Delmiter string (any character in the string delimits fields)
  string delim;

  //! Error output
  void bomb(const char *msg);

public:

  //! Constructor
  ParamParse(istream* in, string Delim);
  
  //! Call with stanza name to set to beginning of stanza's parameter list
  void find_list(const string& stanza);

  //! Return value for given key name.  Returns 1 if found, 0 if not found.
  int find_item(const string& name, string& value);

  //! Return next key/value pair.  Returns 1 if new pair, 0 if end of list.
  int get_next(spair& ret);

  //! Print entire data base to output stream
  void print_database(ostream& out);
};

#endif
