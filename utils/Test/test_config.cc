#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <yaml-cpp/yaml.h>	// YAML support
#include <cxxopts.H>		// Command-line parsing


class EXPparser
{
private:
  std::ostream& ostr;		// Reference to the 'console' output stream
  YAML::Node root;		// Root node read from YAMLloadFile
  int wid;			// Pretty-printing width (default: 0)

  //! Check for EXP stanzas
  void exp_check(const YAML::Node& root);

  //! Parse a YAML file
  void exp_parse(const YAML::Node& cur, int level=0, bool seq=false);

  //! Level spacer for console and output files
  void spacer(int lev, bool seq=false)
  {
    for (int i=0; i<lev; i++) ostr << "  ";
    if (seq) ostr << "- ";
    else     ostr << "  ";
  }

public:

  /** Constructor.  Instatiate with a file and output stream.  The
      stream may be a file or std::cout.  If the file is labeled
      'bad', it's a good /dev/null */
  EXPparser(std::string file, std::ostream& ostr) : ostr(ostr), wid(0)
  {
    root = YAML::LoadFile(file);
  }

  //! Parse a YAML file
  void parse() { exp_parse(root); }

  //! Check for EXP stanzas
  void check() { exp_check(root); }

  //! Set pretty-printing width
  int setPretty(int width=24) {
    int last = wid;
    wid = width;
    return last;
  }
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

  // Iterate over all top-level nodes
  //
  for (YAML::const_iterator it=root.begin(); it!=root.end(); it++) {
    // The key
    auto lab = it->first.as<std::string>();

    // Look for the key in the EXP stanza list
    auto loc = std::find(stanzas.begin(), stanzas.end(), lab);

    // Save it to the unexpected list, if not found
    if (loc == stanzas.end()) unexpected.push_back(lab);

    // Remove the found key from the working list (remain)
    else {
      auto loc = std::find(remain.begin(), remain.end(), lab);
      if (loc == remain.end()) duplicates.push_back(lab);
      else remain.erase(loc);
    }
  }
    
  std::ostringstream sout;

  // If keys remain, they may be foreign to EXP or optional
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

    // Throw an error; this will break EXP
    if (remain.size()) {
      sout << "The following required stanzas were not found:";
      for (auto s : remain) sout << " " << s;
      throw std::runtime_error(sout.str());
    }
  }

  // Print info message about foreign stanzas.  These might be comment
  // or provenance metadata, e.g.
  if (unexpected.size()) {
    if (ostr.rdbuf()  == std::cout.rdbuf())
      std::cout << std::string(70, '=') << std::endl;
    std::cout << "The following stanzas are not used by EXP:";
    for (auto s : unexpected) std::cout << " " << s;
    std::cout << std::endl;
  }

  // Throw an error; this will most likely break EXP
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
      ostr << std::setw(std::max<int>(0, wid-2*level)) << std::left
	   << it->first.as<std::string>() + ": ";
      if (it->second.IsMap() or it->second.IsSequence()) ostr << std::endl;
      exp_parse(it->second, level+1);
    }
  } else if (cur.IsSequence()) {
    for (std::size_t i=0; i< cur.size(); i++) {
      spacer(level+1, true);	// A bit awkward but shorter
      exp_parse(cur[i], level+1, true);
    }
  } else if (cur.IsScalar()) {	// Print the scalar value and done
    std::string lab = cur.as<std::string>();
    if (seq)
      ostr << std::setw(std::max<int>(0, wid-2*level)) << std::left
	   << lab << std::endl;
    else
      ostr << lab << std::endl;
  } else if (cur.IsNull()) {	// Print LF and done
    ostr << std::endl;
  }

  level += 1;
}

int main(int argc, char **argv)
{
  std::string file, outyaml;
  int width;

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
    ("w,width",   "Pretty-printing field width (default: 0, try 24)",
     cxxopts::value<int>(width))
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
    // This makes 'ofs' a proxy for /dev/null
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
    if (vm.count("width")) p.setPretty(width);
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
