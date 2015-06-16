/**
   Compile string:
   g++ -g -o testConfig test_config.cc -lboost_program_options
*/

#include <iostream>
#include <string>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/program_options.hpp>
#include <locale>

using boost::property_tree::ptree;
using namespace std;

namespace po = boost::program_options;

class configuration 
{
  ptree pt;

public:

  enum pType {XML, JSON};

  //! constructor
  configuration(const string& filename, const std::string& type)   
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

  
  // Write current config to a file
  void save(const string& filename, const std::string& type) 
  {
    switch (parse(type)) {
    case JSON:
      write_json(filename+".json", pt);
      break;
    case XML:
      {
	boost::property_tree::xml_writer_settings<char> w(' ',2);
	write_xml(filename+".xml", pt, locale(), w); 
      }
      break;
    default:
      std::cerr << "No such type" << std::endl;
    }
  }
    
  // import export property tree function
  ptree property_tree(void) { return pt; }
  
  void display(void) { display(0, pt); }

private:

  // Parse the enum
  static pType parse(const std::string& stype) 
  {
    pType etype;
    if (stype.compare("JSON") == 0)
      etype = configuration::JSON;
    else if (stype.compare("XML") == 0)
      etype = configuration::XML;
    else {
      std::cerr << "No such type: " << stype << std::endl;
    }
    return etype;
  }

  void display(const int depth, const ptree& tree) 
  {
    BOOST_FOREACH( ptree::value_type const&v, tree.get_child("") ) {
      ptree subtree = v.second;
      string nodestr = tree.get<string>(v.first);
      
      // print current node
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
      display(depth+1,subtree);
    }
  }
};

#define DEFAULT_LABEL "not found"

int main(int argc, char **argv)
{
  std::string label, in_type, out_type, dkey, vkey;
  std::vector<std::string> input, used;
  bool testGet = false;
				// Default values
  used.push_back("input");
  used.push_back("output");

  //--------------------------------------------------
  // Declare the supported options.
  //--------------------------------------------------

  po::variables_map vm;
  po::positional_options_description vm_pos;
  po::options_description vm_cmdline;

  po::options_description desc("Available options");
  desc.add_options()
    ("help,h",										"Produce help message")
    ("test,t",										"Use keys to access and print entries")
    ("key,k", 		
     po::value<std::string>(&dkey)->default_value("parameters.output.outchkpt.nint.desc"),
     "Label to look up in the database")
    ("KEY,K", 		
     po::value<std::string>(&vkey)->default_value("parameters.output.outchkpt.nint.value"),
     "Integer value to look up in the database")
    ("inType,i",
     po::value<std::string>(&in_type)->default_value("JSON"),
     "Input configuration type")
    ("outType,o",
     po::value<std::string>(&out_type)->default_value("JSON"),
     "Output configuration type")
    ("input", po::value<std::vector<std::string> >(&input)->composing(), 
     "");

  vm_pos.add("input", 2);

  po::parsed_options parsed = po::command_line_parser(argc, argv)
    .options(desc)
    .positional(vm_pos)
    .allow_unregistered()
    .run();

  try {
    po::store(parsed, vm);
    po::notify(vm);    
  } catch(boost::program_options::error& e){
    std::cerr << "Invalid_option_value exception thrown parsing config file:"
	      << std::endl << e.what() << std::endl;
    return 2;
  } catch(std::exception e){
    std::cerr <<"Exception thrown parsing config file:" 
	      << std::endl << e.what() << std::endl;
    return 2;
  }

  if (vm.count("help")) {
    std::cout << std::setfill('-') << setw(76) << '-' 
	      << std::endl << std::setfill(' ')
	      << "EXP option test parser" << std::endl
	      << std::endl
	      << std::setfill('-') << setw(76) << '-' 
	      << std::endl << std::setfill(' ')
	      << "Usage: " << argv[0] << " [options] file" << std::endl
	      << std::setfill('-') << setw(76) << '-' 
	      << std::endl << std::setfill(' ')
	      << desc << std::endl;
    return 1;
  }

  if (vm.count("test")) testGet = true;

  for (size_t i=0; i<std::min<size_t>(input.size(), 2); i++) 
    used[i] = input[i];

  configuration cfg(used[0], in_type);
  cfg.display();
  cfg.save(used[1], out_type);

  // example: how to grab data from property tree
  //
  if (testGet) {
    std::string lab = 
      cfg.property_tree().get<std::string>(dkey, DEFAULT_LABEL);
    int val = 
      cfg.property_tree().get<int>(vkey, -1);

    std::cout << std::endl;
    std::cout << "Label = " << lab << std::endl;
    std::cout << "Value = " << val << std::endl;
    std::cout << std::endl;
  }

  return 0;
}
