#include <EXPini.H>

int main(int argc, char **argv)
{
  double x, y, z, u;
  std::string config;
  bool b;

  cxxopts::Options options(argv[0], "Test config routine");

  options.add_options()
    ("h,help", "This help message")
    ("expert", "Print out the expert options, too")
    ("v,verbose", "Print verbose debugging output")
    ("t,template", "Make a YAML template config file")
    ("c,config", "Config file", cxxopts::value<std::string>(config)->default_value("config.yaml"))
    ("x,xval", "x value", cxxopts::value<double>(x)->default_value("1.1"))
    ("y,yval", "y value", cxxopts::value<double>(y)->default_value("2.1"))
    ("z,zval", "z value", cxxopts::value<double>(z)->default_value("3.1"))
    ;

  options.add_options("expert")
    ("C,crazy", "Do something totally crazy")
    ("u,uval", "u value", cxxopts::value<double>(u)->default_value("-13.1"))
    ("b,bool", "logic", cxxopts::value<bool>(b)->default_value("false"))
    ;

  auto vm = options.parse(argc, argv);

  if (vm.count("template")) {
    if (vm.count("expert"))
      SaveConfig(vm, options, "template.yaml", {"", "expert"});
    else
      SaveConfig(vm, options, "template.yaml");
    return 0;
  }

  if (vm.count("expert")) {
    std::cout << options.help({"", "expert"}) << std::endl;
    return 1;
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 1;
  }

  if (vm.count("config")) {
    vm = LoadConfig(options, config);
  }
  
  std::cout << "Current config:" << std::endl;
  for (const auto &kv: vm) {
    std::cout << kv.key() << " = " << kv.value() << std::endl;
  }

}
