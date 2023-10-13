#include <iostream>
#include <iomanip>
#include <string>

#include <yaml-cpp/yaml.h>	// YAML support
#include <cxxopts.H>		// Command-line parsing


static bool verbose = false;

void parse(const YAML::Node& cur, int level=0, bool lf=true)
{
  if (lf and verbose) {
    for (int i=0; i<level; i++) std::cout << "  ";
  }

  if (cur.IsMap()) {
    if (lf and verbose and level>0) std::cout << std::endl;
    for (YAML::const_iterator it=cur.begin(); it!=cur.end(); it++) {
      if (lf and verbose) {
	for (int i=0; i<=level; i++) std::cout << "  ";
      }
      lf = true;
      std::string lab = it->first.as<std::string>() + ":";
      if (verbose) std::cout << std::setw(16) << std::left << lab;
      parse(it->second, level+1);
    }
  } else if (cur.IsSequence()) {
    if (verbose) std::cout << std::endl;
    for (std::size_t i=0; i< cur.size(); i++) {
      if (verbose) {
	for (int i=0; i<=level; i++) std::cout << "  ";
	std::cout << "- ";
      }
      parse(cur[i], level+1, false);
    }
  } else if (cur.IsScalar()) {
    if (verbose) {
      std::string lab = cur.as<std::string>();
      std::cout << lab << std::endl;
    }
  } else if (cur.IsNull()) {
    if (verbose) std::cout << std::endl;
  }

  level += 1;
}

int main(int argc, char **argv)
{
  std::string file;

  cxxopts::Options options(argv[0], "Check EXP YAML config file");

  options.add_options()
    ("h,help",    "Produce help message")
    ("v,verbose", "Print parsed values to stdout")
    ("f,file",    "Name of config file", cxxopts::value<std::string>(file))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    return 1;
  } catch (std::exception e){
    std::cerr << "Exception thrown parsing config file:" 
	      << std::endl << e.what() << std::endl;
    return 2;
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }

  if (vm.count("verbose")) {
    verbose = true;
    std::cout << std::string(70, '=') << std::endl;
  }

  try {
    parse(YAML::LoadFile(file));
  }
  catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << "***" << argv[0] << ": parsing failure" << std::endl;
    exit(-1);
  }

  if (vm.count("verbose")) {
    std::cout << std::string(70, '=') << std::endl;
  }

  std::cout << file << " parsed successfully" << std::endl;
}
