#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <yaml-cpp/yaml.h>	// YAML support
#include <cxxopts.H>		// Command-line parsing


class EXPparser
{
private:
  std::ostream& ostr;
  YAML::Node root;
  bool verbose;

  void exp_check(const YAML::Node& root);
  void exp_parse(const YAML::Node& cur, int level=0, bool seq=false);

  // Level spacer
  void spacer(int lev, bool seq=false)
  {
    for (int i=0; i<lev; i++) ostr << "  ";
    if (seq) ostr << "- ";
    else     ostr << "  ";
  }

public:
  //! Constructor
  EXPparser(std::string file, std::ostream& ostr) : ostr(ostr)
  {
    root = YAML::LoadFile(file);
  }

  //! Parse a YAML file
  void parse() { exp_parse(root); }

  //! Check for EXP stanzas
  void check() { exp_check(root); }
};

void EXPparser::exp_check(const YAML::Node& root)
{
  // All EXP stanzas
  const std::vector<std::string> stanzas
    {
      "Global", "Components", "Output", "External", "Interaction"
    };

  // Optional EXP stanzas
  const std::vector<std::string> optional
    {
      "External", "Interaction"
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

  // Check EXP stanzas not found
  if (remain.size()) {

    // Find and remove optional stanzas
    std::vector<std::string> opt;
    for (auto s : optional) {
      auto it = std::find(remain.begin(), remain.end(), s);
      if (it != remain.end()) {
	opt.push_back(s);
	remain.erase(it);
      }
    }
    
    // Print info message for optional stanzas
    if (opt.size()) {
      if (ostr.rdbuf()  == std::cout.rdbuf())
	std::cout << std::string(70, '=') << std::endl;

      std::cout << "The following optional stanzas were not found:";
      for (auto s : opt) std::cout << " " << s;
      std::cout << std::endl;
    }

    if (remain.size()) {
      sout << "The following required stanzas were not found:";
      for (auto s : remain) sout << " " << s;
      throw std::runtime_error(sout.str());
    }
  }

  if (unexpected.size()) {
    if (ostr.rdbuf()  == std::cout.rdbuf())
      std::cout << std::string(70, '=') << std::endl;
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

void EXPparser::exp_parse(const YAML::Node& cur, int level, bool seq)
{
  // Check for node type
  //
  if (cur.IsMap()) {
    for (YAML::const_iterator it=cur.begin(); it!=cur.end(); it++) {
      if (seq) seq = false;	// Clear the seq flag
      else spacer(level);	// Otherwise, print a level spacer
      // Print the key
      ostr << std::left << it->first.as<std::string>() + ": ";
      if (it->second.IsMap() or it->second.IsSequence()) ostr << std::endl;
      exp_parse(it->second, level+1);
    }
  } else if (cur.IsSequence()) {
    for (std::size_t i=0; i< cur.size(); i++) {
      spacer(level+1, true);
      exp_parse(cur[i], level+1, true);
    }
  } else if (cur.IsScalar()) {
    std::string lab = cur.as<std::string>();
    ostr << lab << std::endl;
  } else if (cur.IsNull()) {
    ostr << std::endl;
  }

  level += 1;
}

int main(int argc, char **argv)
{
  std::string file, outyaml;

  cxxopts::Options options(argv[0],
			   "Check EXP YAML config file for problems.  This does not check the validity\nof the parameters keys or values but only the YAML syntax.  A successful\ncompletion implies that your config file will be parsed correctly by EXP.\n");

  options.add_options()
    ("h,help",    "Produce help message")
    ("v,verbose", "Print parsed values to stdout")
    ("n,noEXP",   "Check general YAML without EXP check")
    ("f,file",    "Name of YAML config file; alternative to final positional argument",
     cxxopts::value<std::string>(file))
    ("o,outyaml", "Store parsed YAML into this file",
     cxxopts::value<std::string>(outyaml))
    ;

  options.parse_positional({"file"});
  options.positional_help("[<YAML file>]");
  options.show_positional_help();

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

  std::ofstream ofs;

  if (vm.count("outyaml")) {
    ofs.open(outyaml);
    if (not ofs) {
      std::cout << argv[0] << ": could not open <" << outyaml
		<< "> for YAML output" << std::endl;
      exit(-1);
    }
  } else {
    ofs.setstate(std::ios_base::badbit);
  }

  bool verbose = false;
  if (vm.count("verbose")) {
    verbose = true;
    std::cout << std::string(70, '=') << std::endl;

    ofs.copyfmt(std::cout);
    ofs.clear(std::cout.rdstate());
    ofs.basic_ios<char>::rdbuf(std::cout.rdbuf());
  }

  // Finally, do the work ...
  //
  try {
    EXPparser p(file, ofs);
    p.parse();
    if (vm.count("noEXP")==0) p.check();
  }
  catch (const std::runtime_error& err) {
    if (verbose) std::cout << std::string(70, '=') << std::endl;
    std::cout << err.what() << std::endl;
    std::cout << "***" << argv[0] << ": parsing failure" << std::endl;
    exit(-1);
  }

  if (verbose) std::cout << std::string(70, '=') << std::endl;

  std::cout << file << " parsed successfully" << std::endl;
}
