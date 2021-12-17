#ifndef _EXPini_H
#define _EXPini_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>

#include <yaml-cpp/yaml.h>	// YAML support
#include <cxxopts.H>		// Option parser

//! Saves the current command-line cxxopts database in YAML format for
//! user configuration
void SaveConfig(const cxxopts::ParseResult& vm,  // Parsed options
		const cxxopts::Options& options, // Default options
		const std::string& config,       // Save file
		// Option groups to save, default by default
		const std::vector<std::string> groups={""},
		// Parameters to exclude
		const std::vector<std::string> exclude={"template", "expert"})
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
	      out << YAML::Key   << kv.key();
	      out << YAML::Value << kv.value();
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

  int count = conf.size()*2, cnt = 0;
  char* data[count];
  
  for (auto it=conf.begin(); it!=conf.end(); it++, cnt+=2) {
    std::ostringstream s1, s2;
    s1 << "--" << it->first.as<std::string>();
    s2 << it->second.as<std::string>();
    
    data[cnt+0] = new char [s1.str().size()+1];
    data[cnt+1] = new char [s2.str().size()+1];
    
    strcpy(data[cnt+0], s1.str().c_str());
    strcpy(data[cnt+1], s2.str().c_str());
  }

  auto vm = options.parse(count, &data[0]);

  for (auto & v : data) delete [] v;

  return vm;
}

#endif
