#include <boost/random/mersenne_twister.hpp>
#include "Configuration.H"

#define DEFAULT_LABEL "not found"

string outdir, runtag;
boost::mt19937 random_gen;

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

  Configuration cfg(used[0], in_type);

  ptree tree = cfg.property_tree();

  std::cout << "Top level nodes" << std::endl 
	    << "---------------" << std::endl;

  for (auto &v : tree.get_child("") ) {
    // Get the node name
    string nodestr = tree.get<string>(v.first);
      
    // print current node
    //
    if ( nodestr.length() > 0 ) {
      std::cout << string("").assign(4,' ') << "* ";
      std::cout << v.first;
      std::cout << "=\"" << tree.get<string>(v.first) << "\"";
      std::cout << std::endl;
    } else if (v.first.length()) {
      std::cout << string("").assign(4,' ') << "* ";
      std::cout << v.first << std::endl;
    } else {
      std::cout << string("").assign(4,' ') << "* ";
      std::cout << v.second.data() << std::endl;
    }
  }

  std::cout << std::endl
	    << "Example get" << std::endl 
	    << "-----------" << std::endl << std::left
	    << std::setw(20) << "Name" 
	    << std::setw(20) << "Value" << "Comment field"  << std::endl
	    << std::setw(20) << "-----" 
	    << std::setw(20) << "-----" << "-------------"  << std::endl
	    << std::setw(20) << "param1" 
	    << std::setw(20) << tree.get<int>("param1.value", -1) 
	    << tree.get<string>("param1.desc", "") << std::endl
	    << std::setw(20) << "param2"
	    << std::setw(20) << tree.get<float>("param2.value", 0.0)
	    << tree.get<string>("param2.desc", "") << std::endl
	    << std::setw(20) << "param3"
	    << std::setw(20) << tree.get<string>("param3.value", "No string") 
	    << tree.get<string>("param3.desc", "[None]") << std::endl
	    << std::setw(20) << "param4"
	    << std::setw(20) << tree.get<int>("param4.value", -1)
	    << tree.get<string>("param4.desc", "") << std::endl
	    << std::endl;

  return 0;
}
