#ifndef _SUTILS_H
#define _SUTILS_H

#include <algorithm> 
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <locale>

//! Split string on a character delimiter
std::vector<std::string> str_split(const std::string& s, char delimiter);

//! Return a lower case copy
std::string str_to_lower(const std::string& s);

//! Return an upper case copy
std::string str_to_upper(const std::string& s);

//! String trim white-space from start (in place version)
void ltrim(std::string &s);

//! String trim white-space string from end (in place version)
void rtrim(std::string &s);

//! String trim white-spacefrom both ends (in place)
void trim(std::string &s);

//! String trim white space from start (new copy)
std::string ltrim_copy(std::string s);

//! String trim white space from end (new copy)
std::string rtrim_copy(std::string s);

//! String trim white space from both ends (new copy)
std::string trim_copy(std::string s);


#endif
