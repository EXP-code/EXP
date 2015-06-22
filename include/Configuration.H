#ifndef _Configuration_H
#define _Configuration_H

#include <iostream>
#include <string>
#include <fstream>
#include <locale>

#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/program_options.hpp>

using boost::property_tree::ptree;
using namespace std;

namespace po = boost::program_options;

class Configuration 
{
  ptree pt;

public:

  enum pType {XML, JSON};

  //! Null constructor
  Configuration() {}

  //! constructor
  Configuration(const string& filename, const std::string& type)   
  { 
    load(filename, type);
  }

  //! Read in a config
  void load(const string& filename, const std::string& type)   
  { 
    ifstream input(filename.c_str());
    switch (parse(type)) {
    case JSON:
      read_json(input, pt);
      break;
    case XML:
    default:
      read_xml(input, pt, boost::property_tree::xml_parser::trim_whitespace );
      break;
    }
  }

  
  //! Write current config to a file
  void save(const string& filename, const std::string& type) 
  {
    switch (parse(type)) {
    case JSON:
      write_json(filename+".json", pt);
      break;
    case XML:
      {
#if BOOST_VERSION >= 105600
	boost::property_tree::xml_writer_settings<std::string> w(' ', 2);
#else
	boost::property_tree::xml_writer_settings<char> w(' ', 2);
#endif
	write_xml(filename+".xml", pt, locale(), w); 
      }
      break;
    default:
      std::cerr << "No such type" << std::endl;
    }
  }
    
  //! import/export property tree function
  ptree& property_tree(void) { return pt; }
  
  void display(void) { display(0, pt); }

  //! Access/set a parameter value
  template<typename T>
  T entry(const std::string& query, const std::string& desc, const T& def)
  {
    if (pt.count(query)) {
      return pt.get<T>(query+".value");
    } else {
      pt.put(query+".desc", desc);
      pt.put(query+".value", def);
      return def;
    }
  }

private:

  // Parse the enum
  static pType parse(const std::string& stype) 
  {
    pType etype(JSON);
    if (stype.compare("JSON") == 0)
      etype = Configuration::JSON;
    else if (stype.compare("XML") == 0)
      etype = Configuration::XML;
    else {
      std::cerr << "No such type: " << stype << std::endl;
    }
    return etype;
  }

  void display(const int depth, const ptree& tree) 
  {
    for (auto const &v : tree.get_child("") ) {
      // Get the subtree
      ptree subtree  = v.second;

      // Get the node name
      string nodestr = tree.get<string>(v.first);
      
      // print current node
      //
      if ( nodestr.length() > 0 ) {
	std::cout << string("").assign(depth*4,' ') << "* ";
	std::cout << v.first;
	std::cout << "=\"" << tree.get<string>(v.first) << "\"";
	std::cout << std::endl;
      } else if (v.first.length()) {
	std::cout << string("").assign(depth*4,' ') << "* ";
	std::cout << v.first << std::endl;
      } else {
	std::cout << string("").assign(depth*4,' ') << "* ";
	std::cout << v.second.data() << std::endl;
      }

      // recursive go down the hierarchy
      //
      display(depth+1,subtree);
    }
  }
};

#endif
