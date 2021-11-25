#ifndef _SUTILS_H
#define _SUTILS_H

#include <algorithm> 
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <locale>

//! Split string on a character delimiter
std::vector<std::string> str_split(const std::string& s, char delimiter)
{
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

//! Return a lower case copy
std::string str_to_lower(const std::string& s)
{
  std::string d(s);
  std::for_each(d.begin(), d.end(), [](char & c){ c = std::tolower(c); });

  return d;
}

//! Return an upper case copy
std::string str_to_upper(const std::string& s)
{
  std::string d(s);
  std::for_each(d.begin(), d.end(), [](char & c){ c = std::toupper(c); });

  return d;
}

//! String trim white-space from start (in place version)
static inline void ltrim(std::string &s)
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
    return !std::isspace(ch);
  }));
}

//! String trim white-space string from end (in place version)
static inline void rtrim(std::string &s)
{
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
    return !std::isspace(ch);
  }).base(), s.end());
}

//! String trim white-spacefrom both ends (in place)
static inline void trim(std::string &s)
{
  ltrim(s);
  rtrim(s);
}

//! String trim white space from start (new copy)
static inline std::string ltrim_copy(std::string s)
{
  ltrim(s);
  return s;
}

//! String trim white space from end (new copy)
static inline std::string rtrim_copy(std::string s)
{
  rtrim(s);
  return s;
}

//! String trim white space from both ends (new copy)
static inline std::string trim_copy(std::string s)
{
  trim(s);
  return s;
}


#endif
