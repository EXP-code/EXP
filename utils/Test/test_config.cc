#include <iostream>
#include <iomanip>
#include <string>

#include <yaml-cpp/yaml.h>	// YAML support
#include <cxxopts.H>		// Command-line parsing


static bool verbose = false;

void stanza_check(const YAML::Node& root)
{
  const std::vector<std::string> stanzas
    {
      "Global", "Components", "Output", "External", "Interaction"
    };

  std::vector<std::string> unexpected, remain(stanzas), duplicates;

  for (YAML::const_iterator it=root.begin(); it!=root.end(); it++) {
    auto lab = it->first.as<std::string>();
    auto loc = std::find(stanzas.begin(), stanzas.end(), lab);
    if (loc == stanzas.end()) unexpected.push_back(lab);
    else {
      auto loc = std::find(remain.begin(), remain.end(), lab);
      if (loc == remain.end()) duplicates.push_back(lab);
      else remain.erase(loc);
    }
  }
    
  std::ostringstream sout;

  if (remain.size()) {
    sout << "The following stanzas were not found:";
    for (auto s : remain) sout << " " << s;
    throw std::runtime_error(sout.str());
  }

  if (unexpected.size()) {
    if (verbose) std::cout << std::string(70, '=') << std::endl;
    std::cout << "The following stanzas are not used by EXP:";
    for (auto s : unexpected) std::cout << " " << s;
    std::cout << std::endl;
  }

  if (duplicates.size()) {
    sout << "The following stanzas are duplicated:";
    for (auto s : duplicates) sout << " " << s;
    throw std::runtime_error(sout.str());
  }

}

void parse(const YAML::Node& cur, int level=0, bool seq=false)
{
  // Key spacer
  auto spacer = [&level]() {
    if (verbose) {
      for (int i=0; i<level; i++) std::cout << "  ";
    }
  };

  // Initial key spacer at this level unless we're a sequence
  // if (not seq) spacer();

  // Check for node type
  if (cur.IsMap()) {
				// End the previous level
    if (not seq and verbose and level>0)
      std::cout << std::endl;

    for (YAML::const_iterator it=cur.begin(); it!=cur.end(); it++) {
      if (not seq) spacer();
      seq = false;
      if (verbose) std::cout << std::left << it->first.as<std::string>() + ": ";
      parse(it->second, level+1);
    }
  } else if (cur.IsSequence()) {
    if (not seq and verbose) std::cout << std::endl;
    for (std::size_t i=0; i< cur.size(); i++) {
      if (not seq) spacer();
      if (verbose) std::cout << "- ";
      parse(cur[i], level+1, true);
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
    ("n,noEXP",   "Check general YAML without EXP check")
    ("f,file",    "Name of config file", cxxopts::value<std::string>(file))
    ;

  options.parse_positional({"file"});

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

  YAML::Node root = YAML::LoadFile(file);

  try {
    parse(root);
    if (vm.count("noEXP")==0) stanza_check(root);
  }
  catch (const std::runtime_error& err) {
    if (verbose) std::cout << std::string(70, '=') << std::endl;
    std::cout << err.what() << std::endl;
    std::cout << "***" << argv[0] << ": parsing failure" << std::endl;
    exit(-1);
  }

  if (vm.count("verbose"))
    std::cout << std::string(70, '=') << std::endl;

  std::cout << file << " parsed successfully" << std::endl;
}
