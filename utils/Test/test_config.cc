#include <EXPini.H>

int main(int argc, char **argv)
{
  double x, y, z;
  std::string config;

  cxxopts::Options options(argv[0], "Test config routine");

  options.add_options()
    ("h,help", "This help message")
    ("v,verbose", "Print verbose debugging output")
    ("t,template", "Make a YAML template config file")
    ("c,config", "Config file", cxxopts::value<std::string>(config)->default_value("config.yaml"))
    ("x,xval", "x value", cxxopts::value<double>(x)->default_value("1.1"))
    ("y,yval", "y value", cxxopts::value<double>(y)->default_value("2.1"))
    ("z,zval", "z value", cxxopts::value<double>(z)->default_value("3.1"))
    ;

  auto vm = options.parse(argc, argv);

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 1;
  }

  if (vm.count("template")) {
    SaveConfig(vm, "template.yaml");
    return 0;
  }

  if (vm.count("config")) {
    vm = LoadConfig(options, config);
  }
  
  std::cout << "Current config:" << std::endl;
  for (const auto &kv: vm) {
    std::cout << kv.key() << " = " << kv.value() << std::endl;
  }

}
