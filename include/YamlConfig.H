#ifndef _Config_H
#define _Config_H

#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>

#include <yaml-cpp/yaml.h>	// YAML support
#include "cxxopts.H"		// Option parser

/** Saves the current command-line cxxopts database in YAML format for
    user configuration.  For ease, the current convention is that all
    vector-values options will be stored as strings, and converted by
    the calling program to the correct value internally.
*/
void SaveConfig(		// Parsed options
		const cxxopts::ParseResult& vm,
				// Default options
		const cxxopts::Options& options,
				// Save file
		const std::string& config,
				// Option groups to save, default by default
		const std::vector<std::string> groups={""},
				// Parameters to exclude
		const std::vector<std::string> exclude={"help", "template", "expert", "config"});

//! Read the YAML parameter config file and load the cxxopts database
//! with parameters
cxxopts::ParseResult LoadConfig(cxxopts::Options& options,
				const std::string& config);

#endif
