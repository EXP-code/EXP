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
void SaveConfig(const cxxopts::ParseResult& vm,
		const std::string& config)
{    
  YAML::Emitter out;

  out << YAML::BeginMap;
  YAML::Node node;
  for (const auto &kv: vm) {
    out << YAML::Key   << kv.key();
    out << YAML::Value << kv.value();
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
