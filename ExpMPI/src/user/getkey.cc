#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <StringTok.H>


string trimLeft(const string &value)
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


string trimRight(const string &value)
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


string trimComment(const string &value)
 {
   string::size_type where = value.find_first_of("#!");
  
   if (where == string::npos)
     // no comments
     return value;
   
   return value.substr(0, where - 1);
 }


int parse_string(string& line, vector<string>& word, vector<string>& valu)
{
    				// Look for trailing comment
  line = trimComment(line);
  line = trimLeft(trimRight(line));

				// Look for leading . or zero size
  if (line.substr(0,1) == "." || line.size() == 0) return 0;

				// Is there a token marker?
  if (line.find("=") == string::npos) return 0;

				// Parse for tokens
  StringTok<string> tokens(line);
  string first = trimLeft(trimRight(tokens("=")));
  string second = trimLeft(trimRight(tokens("=")));
    
  if (!first.empty()) {
    word.push_back(first);
    valu.push_back(second);
    return 1;
  } 

  return 0;
} 



int get_key_value_from_file(string file, 
			    vector<string>& word, vector<string>& valu)
{
  ifstream in(file.c_str());
  if (!in) {
    cerr << "get_key_value_from_file: error reading <" << file << ">\n";
    return 0;
  }
				// Start with empty elements
  word.erase(word.begin(), word.end());
  valu.erase(word.begin(), word.end());

  int num=0;
  const int lbufsize=1024;
  char lbuf[lbufsize];

  while (in) {
    in.getline(lbuf, lbufsize);
    if (!in) break;
    string line(lbuf);
    num += parse_string(line, word, valu);
  }

  return num;
}


int get_key_value(int argc, char **argv, 
		  vector<string>& word, vector<string>& valu)
{
				// Start with empty elements
  word.erase(word.begin(), word.end());
  valu.erase(word.begin(), word.end());

  int num=0;
  for (int i=0; i<argc; i++) {
    string line(argv[i]);
    num += parse_string(line, word, valu);
  }

  return num;
}


