#ifndef _EXPini_H
#define _EXPini_H

#include <algorithm>
#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>

#include <yaml-cpp/yaml.h>	// YAML support
#include <cxxopts.H>		// Option parser

/** Saves the current command-line cxxopts database in YAML format for
    user configuration.  For ease, the current convention is that all
    vector-values options will be stored as strings, and converted by
    the calling program to the correct value internally.
*/
void SaveConfig(const cxxopts::ParseResult& vm,  // Parsed options
		const cxxopts::Options& options, // Default options
		const std::string& config,       // Save file
		// Option groups to save, default by default
		const std::vector<std::string> groups={""},
		// Parameters to exclude
		const std::vector<std::string> exclude={"help", "template", "expert", "config"})
{    
  YAML::Emitter out;

  out << YAML::BeginMap;
  YAML::Node node;

  auto sep = std::string(20, '-');

  for (const auto kv: vm) {

    // Iterate through option groups
    for (auto g : options.groups()) {

      // Look for group in include list
      auto it = std::find(groups.begin(), groups.end(), g);

      // Found it: look for key in option list for this group
      if (it != groups.end()) {
	
	for (auto m : options.group_help(g).options) {
	  if (m.l == kv.key()) {
	    // Is this key in the excluded list?
	    auto jt = std::find(exclude.begin(), exclude.end(), kv.key());

	    // Not excluded: write to template file with desc as comments
	    if (jt == exclude.end()) {
	      // Write the key and then check the value
	      out << YAML::Key << kv.key();
	      // Are we a vector?
	      if (typeid(kv.value()) == typeid(std::vector<std::string>)) {
		out << YAML::BeginSeq;
		out << kv.value();
		out << YAML::EndSeq;
	      }
	      // We are a scalar
	      else {
		out << YAML::Value << kv.value();
	      }
	      out << YAML::Comment(m.desc);
	    }
	  }
	}
      }
    }
  }
  out << YAML::EndMap;
  
  std::ofstream temp(config);
  if (temp)
    temp << out.c_str();
  else
    std::cerr << "Could not save template file <" << config << ">"
	      << std::endl;
}


//! Read the YAML parameter config file and load the cxxopts database
//! with parameters
cxxopts::ParseResult LoadConfig(cxxopts::Options& options,
				const std::string& config)
{
  YAML::Node conf = YAML::LoadFile(config);

  const int count = conf.size()*2+1;
  std::vector<char*> data(count);
  int cnt = 1;

  data[0] = new char [11];
  strcpy(data[0], "LoadConfig"); // Emulate the caller name
  
  for (auto it=conf.begin(); it!=conf.end(); it++) {
    std::ostringstream s1, s2;
    s1 << "--" << it->first.as<std::string>();
    
    // Are we vector valued?
    if (it->second.IsSequence()) {
      std::ostringstream sout;
      for (auto jt=it->second.begin(); jt!=it->second.end(); jt++) {
	if (jt != it->second.begin()) s2 << ','; // Add the comma
	s2 << jt->as<std::string>();		 // Add the element
      }
    }
    // We are scalar valued
    else {
      s2 << it->second.as<std::string>();
    }

    // Check key for a yaml boolean type and exclude false values to
    // comply with cxxopts policy
    //
    std::string val = s2.str();
    std::transform(val.begin(), val.end(), val.begin(),
		   [](unsigned char c){ return std::tolower(c); });

    if (val.find("false") and val.find("off") and val.find("no")) {
      data[cnt+0] = new char [s1.str().size()+1];
      data[cnt+1] = new char [s2.str().size()+1];
    
      strcpy(data[cnt+0], s1.str().c_str());
      strcpy(data[cnt+1], s2.str().c_str());

      cnt += 2;
    }
  }
  
  auto vm = options.parse(cnt, &data[0]);

  for (auto & v : data) delete [] v;

  return vm;
}

#endif
