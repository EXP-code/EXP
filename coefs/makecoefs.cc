#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <cxxopts.H>
#include <libvars.H>

#include <CoefContainer.H>
#include <BasisFactory.H>
#include <FieldGenerator.H>
#include <ParticleReader.H>

int main(int argc, char **argv)
{
  std::string runtag, outdir, config, type, delim, files, comp;

  // Parse Command line
  //
  const std::string overview = "A quick test of coefficient making and field generation\nfrom one or more snapshots\n";

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v, verbose", "provide details of running pipeline")
    ("e, extend", "add coefficients to an existing HDF5 file")
    ("s, surface", "make VTK field files")
    ("f,files", "list of snapshot files",
     cxxopts::value<std::string>(files)->default_value("file.list"))
    ("t,type", "snapshot type",
     cxxopts::value<std::string>(type)->default_value("PSPspl"))
    ("o,outdir", "output directory",
     cxxopts::value<std::string>(outdir)->default_value("."))
    ("r,runtag", "runtag for output data files",
     cxxopts::value<std::string>(runtag)->default_value("testrun"))
    ("b,basis", "YAML config file for basis",
     cxxopts::value<std::string>(config)->default_value("basis.yaml"))
    ("c,comp", "The component name",
     cxxopts::value<std::string>(comp)->default_value("dark"))
     ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 1;
  }

  bool verbose = vm.count("verbose");
  bool surface = vm.count("surface");

  // Check to see if HDF5 file exists
  //
  bool extend = vm.count("extend");
  bool hdf5_exists = std::filesystem::exists(runtag + ".h5");

  if (hdf5_exists and not extend) {
    std::cout << "HDF5 file <" << runtag + ".h5" << "> exists.  Use --extend flag if you would like to append the new coefficients" << std::endl;

    return(1);
  }

  try {

    // Read the YAML from a file (could use a YAML emitter to do this on
    // the fly inside of the code)
    //
    YAML::Node yaml = YAML::LoadFile(config);
    
    // Create the basis object
    //
    auto basis = Basis::Basis::factory(yaml);
    
    // Will contain the coefficients
    //
    Coefs::CoefsPtr coefs;
    
    // Loop through list of snapshots
    //
    for (auto batch : PR::ParticleReader::parseFileList(files, delim)) {
      
      // Create a particle reader object
      //
      auto reader = PR::ParticleReader::createReader(type, batch);
      reader->SelectType(comp);
      
      // Make the coefficients from this reader
      //
      auto coef = basis->createFromReader(reader);
      
      // Add the coefficients to the container (makes the instance if
      // necessary)
      //
      coefs = Coefs::Coefs::addcoef(coefs, coef); 

      // Some running commentary
      //
      if (vm.count("verbose")) {
	std::cout << "---- Finished the batch:";
	for (auto f : batch) std::cout << " " << f;
	std::cout << std::endl;
      }
    }

    // Write an H5 file
    //
    if (hdf5_exists and extend)
      coefs->ExtendH5Coefs(runtag);
    else
      coefs->WriteH5Coefs(runtag);

    // Make some surface files
    //
    if (surface) {

      std::vector<double> time = coefs->Times();
      std::vector<double> pmin = {-1.0, -1.0, 0.0};
      std::vector<double> pmax = { 1.0,  1.0, 0.0};
      std::vector<int>    grid = { 40,   40,  0};

      Field::FieldGenerator fields(time, pmin, pmax, grid);

      if (verbose) {
	std::cout << "---- Making surfaces for times:";
	for (auto t : time) std::cout << " " << t;
	std::cout << std::endl;
      }

      fields.file_slices(basis, coefs, runtag);
    }
  }
  catch (const std::runtime_error& error) {
    std::cout << "makecoefs: found a problem" << std::endl
	      << error.what() << std::endl;
    return(1);
  }

  // Done
  //
  return(0);
}
