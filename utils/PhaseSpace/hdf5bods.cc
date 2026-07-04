/*
Draft EXP ascii particle phase-space data HDF5 converter

This is a possible template for implementing HDF5 input files for
EXP. The code reads an ASCII file with particle data and writes it to
an HDF5 file using the HighFive library. It also provides a function
to read the HDF5 file back into ASCII format, allowing for round-trip
conversion.

Compilation notes
-----------------

To compile this code, you need to have HighFive and HDF5 development
libraries installed.  It also uses the 'cxxopts.H' and 'ScopeTimer.H'
which need to be local.

Example compilation command for Ubuntu/Debian systems:

g++ -std=c++17 -I/usr/include/hdf5/serial exp2hdf5.cc -o exp2hdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl -lhdf5 -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/serial

Key features
------------

1. Reads ASCII header with 3 integers for particle count, aux int
   fields, and aux float fields
2. Parses flexible number of auxiliary integer and floating-point fields
3. Uses HighFive for a simple C++ HDF5 interface
4. Applies GZIP compression (level 4) for reduced file size.
5. Preserves numerical precision with double-precision floats if desired.
   This is less relevant for simulation initial conditions.
6. Stores metadata as HDF5 attributes for round-trip integrity
7. Comprehensive std error handling with informative messages
8. Supports both float32 and float64 precision with automatic detection
9. Round-trip conversion maintains numerical accuracy when using float64

Performance considerations
--------------------------

- Uses chunked storage with automatic chunk sizing (10th of dataset or 1024,
  whichever is larger).
- Uses shuffle filter to improve compression efficiency.
- GZIP compression level 9 provides good compression ratio. Level 5 seems
  good enough for most cases.
- For very large datasets, consider tuning chunk size in DataSetCreateProps
- Memory usage is O(N) where N = total number of values across all fields

File structure in HDF5
----------------------
/
├── Attributes
│ ├── num_particles (int)
│ ├── num_aux_ints (int)
│ └── num_aux_floats (int)
└── particles/ (Group)
├── m (Dataset)
├── x, y, z (Datasets)
├── u, v, w (Datasets)
├── aux_int_0, aux_int_1, ... (Datasets)
└── aux_float_0, aux_float_1, ... (Datasets)

This schemea could (should be?) changed to use the Gadget style.  On
the other hand, this design is clean and specific to exp.  Also, users
can easily wrie and read this with (e.g.) h5py.  I'm incliened to keep
it.
*/

// C++ std
#include <iostream>
#include <variant>
#include <chrono>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <filesystem>

// HighFive
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

// Command-line parsing
#include "cxxopts.H"

// For command timing
#include "ScopeTimer.H"

using namespace HighFive;

// Specify floating point precision
enum class FloatPrecision { FLOAT32, FLOAT64 };

// Type variant for flexible float handling
template<typename T>
struct ParticleDataTemplate
{
  std::vector<T> m;             // mass/species index
  std::vector<T> x, y, z;       // position
  std::vector<T> u, v, w;       // velocity
  std::vector<std::vector<int>> aux_ints;      // auxiliary integer fields
  std::vector<std::vector<T>> aux_floats;      // auxiliary float fields

  size_t num_particles = 0;
  size_t num_aux_ints = 0;
  size_t num_aux_floats = 0;
};

// Read ASCII and write HDF5 with specified precision
template<typename T>
void ascii_to_hdf5_impl(const std::string& ascii_file, 
                        const std::string& hdf5_file,
                        FloatPrecision precision)
{
  std::ifstream infile(ascii_file);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open ASCII file: " + ascii_file);
  }

  // Read header
  int num_particles, num_aux_ints, num_aux_floats;
  if (!(infile >> num_particles >> num_aux_ints >> num_aux_floats)) {
    throw std::runtime_error("Failed to read header from ASCII file");
  }

  // Validate header
  if (num_particles <= 0 || num_aux_ints < 0 || num_aux_floats < 0) {
    throw std::runtime_error("Invalid header values");
  }

  // Initialize particle data structure
  ParticleDataTemplate<T> data;

  data.num_particles  = num_particles;
  data.num_aux_ints   = num_aux_ints;
  data.num_aux_floats = num_aux_floats;

  data.m.resize(num_particles);
  data.x.resize(num_particles);
  data.y.resize(num_particles);
  data.z.resize(num_particles);
  data.u.resize(num_particles);
  data.v.resize(num_particles);
  data.w.resize(num_particles);

  data.aux_ints.resize(num_aux_ints);
  for (auto& vec : data.aux_ints) {
    vec.resize(num_particles);
  }

  data.aux_floats.resize(num_aux_floats);
  for (auto& vec : data.aux_floats) {
    vec.resize(num_particles);
  }

  // Read particle data
  for (int i = 0; i < num_particles; ++i) {
    if (!(infile >> data.m[i] >> data.x[i] >> data.y[i] >> data.z[i]
              >> data.u[i] >> data.v[i] >> data.w[i])) {
      throw std::runtime_error("Failed to read particle data at row " + 
                               std::to_string(i));
    }

    for (int j = 0; j < num_aux_ints; ++j) {
      if (!(infile >> data.aux_ints[j][i])) {
        throw std::runtime_error("Failed to read aux int field at row " + 
                                 std::to_string(i));
      }
    }

    for (int j = 0; j < num_aux_floats; ++j) {
      if (!(infile >> data.aux_floats[j][i])) {
        throw std::runtime_error("Failed to read aux float field at row " + 
                                 std::to_string(i));
      }
    }
  }
  infile.close();

  // Create HDF5 file with compression enabled
  File file(hdf5_file, File::ReadWrite | File::Create | File::Truncate);

  // Store header information and precision metadata
  file.createAttribute<int>("num_particles", DataSpace::From(num_particles))
    .write(num_particles);
  file.createAttribute<int>("num_aux_ints", DataSpace::From(num_aux_ints))
    .write(num_aux_ints);
  file.createAttribute<int>("num_aux_floats", DataSpace::From(num_aux_floats))
    .write(num_aux_floats);

  // Store precision type: 0 = float32, 1 = float64
  int precision_flag = (precision == FloatPrecision::FLOAT32) ? 0 : 1;
  file.createAttribute<int>("float_precision", DataSpace::From(precision_flag))
    .write(precision_flag);

  // Create a group for particle data
  Group particles_group = file.createGroup("particles");

  // Define compression filters
  DataSetCreateProps props;
  props.add(Chunking(std::vector<hsize_t>{(hsize_t)std::max(num_particles / 10, 1024)}));
  props.add(Shuffle());
  props.add(Deflate(4));	// This is a good compromise between
				// speed and compression ratio

  // Write core phase-space fields
  particles_group.createDataSet("m", data.m, props);
  particles_group.createDataSet("x", data.x, props);
  particles_group.createDataSet("y", data.y, props);
  particles_group.createDataSet("z", data.z, props);
  particles_group.createDataSet("u", data.u, props);
  particles_group.createDataSet("v", data.v, props);
  particles_group.createDataSet("w", data.w, props);

  // Write auxiliary integer fields
  for (size_t j = 0; j < data.aux_ints.size(); ++j) {
    std::string dset_name = "aux_int_" + std::to_string(j);
    particles_group.createDataSet(dset_name, data.aux_ints[j], props);
  }

  // Write auxiliary float fields
  for (size_t j = 0; j < data.aux_floats.size(); ++j) {
    std::string dset_name = "aux_float_" + std::to_string(j);
    particles_group.createDataSet(dset_name, data.aux_floats[j], props);
  }

  std::string precision_str = (precision == FloatPrecision::FLOAT32) ? "float32" : "float64";
  std::cout << "Successfully wrote " << num_particles << " particles to " 
            << hdf5_file << " (" << precision_str << ")" << std::endl;
}

// Dispatcher function - user specifies precision
void ascii_to_hdf5(const std::string& ascii_file, 
                   const std::string& hdf5_file,
                   FloatPrecision precision = FloatPrecision::FLOAT64)
{
  if (precision == FloatPrecision::FLOAT32) {
    ascii_to_hdf5_impl<float>(ascii_file, hdf5_file, precision);
  } else {
    ascii_to_hdf5_impl<double>(ascii_file, hdf5_file, precision);
  }
}

// Variant for flexible reading
using FloatData = std::variant<std::vector<float>, std::vector<double>>;

struct ParticleDataVariant
{
  FloatData m, x, y, z, u, v, w;
  std::vector<std::vector<int>> aux_ints;
  std::vector<FloatData> aux_floats;

  size_t num_particles = 0;
  size_t num_aux_ints = 0;
  size_t num_aux_floats = 0;
  FloatPrecision precision;
};

// Read HDF5 with automatic precision detection
ParticleDataVariant read_hdf5_data(const std::string& hdf5_file)
{
  File file(hdf5_file, File::ReadOnly);

  // Read metadata
  int num_particles = file.getAttribute("num_particles").read<int>();
  int num_aux_ints = file.getAttribute("num_aux_ints").read<int>();
  int num_aux_floats = file.getAttribute("num_aux_floats").read<int>();

  // Auto-detect precision
  int precision_flag = 1;  // Default to float64
  try {
    precision_flag = file.getAttribute("float_precision").read<int>();
  } catch (...) {
    // If attribute doesn't exist, assume float64 (backward compatibility)
  }
  FloatPrecision precision = (precision_flag == 0) ? FloatPrecision::FLOAT32 
                                                    : FloatPrecision::FLOAT64;

  // Open particles group
  Group particles_group = file.getGroup("particles");

  ParticleDataVariant data;
  data.num_particles  = num_particles;
  data.num_aux_ints   = num_aux_ints;
  data.num_aux_floats = num_aux_floats;
  data.precision      = precision;

  // Lambda to read float or double based on precision
  auto read_dataset = [&particles_group, precision](const std::string& name) -> FloatData {
    if (precision == FloatPrecision::FLOAT32) {
      return particles_group.getDataSet(name).read<std::vector<float>>();
    } else {
      return particles_group.getDataSet(name).read<std::vector<double>>();
    }
  };

  // Read core PSP fields
  data.m = read_dataset("m");
  data.x = read_dataset("x");
  data.y = read_dataset("y");
  data.z = read_dataset("z");
  data.u = read_dataset("u");
  data.v = read_dataset("v");
  data.w = read_dataset("w");

  // Read auxiliary integer fields
  data.aux_ints.resize(num_aux_ints);
  for (int j = 0; j < num_aux_ints; ++j) {
    std::string dset_name = "aux_int_" + std::to_string(j);
    data.aux_ints[j] = particles_group.getDataSet(dset_name)
                                      .read<std::vector<int>>();
  }

  // Read auxiliary float fields
  data.aux_floats.resize(num_aux_floats);
  for (int j = 0; j < num_aux_floats; ++j) {
    std::string dset_name = "aux_float_" + std::to_string(j);
    data.aux_floats[j] = read_dataset(dset_name);
  }

  return data;
}

// Write ASCII from variant data using std::visit
//
void hdf5_to_ascii(const std::string& hdf5_file, const std::string& ascii_file)
{
  ParticleDataVariant data = read_hdf5_data(hdf5_file);

  std::ofstream outfile(ascii_file);
  if (!outfile.is_open()) {
    throw std::runtime_error("Could not open ASCII file for writing: " + ascii_file);
  }

  // Write header
  outfile << data.num_particles << " " << data.num_aux_ints 
          << " " << data.num_aux_floats << "\n";

  outfile.precision(16);  // High precision for floats

  // Lambda visitor to handle both float and double
  auto write_value = [](const FloatData& variant_data, size_t index) {
    return std::visit([index](const auto& vec) { return static_cast<double>(vec[index]); },
                      variant_data);
  };

  // Write particle data
  int num_particles = data.num_particles;
  for (int i = 0; i < num_particles; ++i) {
    // Write core phase-space fields using visitor
    outfile << write_value(data.m, i) << " "
            << write_value(data.x, i) << " " 
            << write_value(data.y, i) << " " 
            << write_value(data.z, i) << " "
            << write_value(data.u, i) << " " 
            << write_value(data.v, i) << " " 
            << write_value(data.w, i);

    // Write auxiliary integer fields
    for (int j = 0; j < data.num_aux_ints; ++j) {
      outfile << " " << data.aux_ints[j][i];
    }

    // Write auxiliary float fields using visitor
    for (int j = 0; j < data.num_aux_floats; ++j) {
      outfile << " " << write_value(data.aux_floats[j], i);
    }
    outfile << "\n";
  }
  outfile.close();

  std::string precision_str = (data.precision == FloatPrecision::FLOAT32) 
                              ? "float32" : "float64";
  std::cout << "Successfully wrote " << num_particles << " particles to " 
            << ascii_file << " (" << precision_str << ")" << std::endl;
}

int main(int argc, char* argv[])
{
  std::string prefix = "particles", output;
  std::string suffix = "bods";
  std::string ascii_restored = "particles_restored.txt";

  // Parse command-line arguments for input/output files and mode
  //
  cxxopts::Options options("exp2hdf5", "ASCII to HDF5 particle converter with round-trip support and precision handling");
  
  options.add_options()
    ("i,input", "Input prefix", cxxopts::value<std::string>(prefix)->default_value("particles"))
    ("o,output", "Output prefix (optional, otherwise input prefix is used)", cxxopts::value<std::string>(prefix))
    ("a,suffix", "Input suffix", cxxopts::value<std::string>(suffix)->default_value("bods"))
    ("roundtrip", "Perform round-trip conversion (ASCII -> HDF5 -> ASCII)")
    ("verify", "Verify that the restored ASCII file matches the original in 'roundtrip' mode")
    ("to_hdf5", "Convert ASCII to HDF5")
    ("double", "Use double precision for HDF5 output (float is default)")
    ("to_ascii", "Convert HDF5 to ASCII");
    ("h,help", "Print usage");

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;

    // Append custom examples
    std::cout << R"(
Examples:
  Convert a standard EXP input body file nameed 'mybods.bods' to hdf5 format
    $ hdf5bods --to_hdf5 -i mybods
  The resulting HDF5 file using float32 internally will be called 'mybods.h5'

  Do the same conversion but use full double precision
    $ hdf5bods --to_hdf5 --double -i mybods
  The resulting HDF5 file will use float64 internally

  Convert the HDF5 file back to the standard EXP ascii format:
    $ hdf5bods --to_ascii -i input

  The suffix on input can be customized.  For example, the following converts
  a standard EXP input body file nameed 'mybods.asc' to hdf5 format
    $ hdf5bods --to_hdf5 --suffix=asc -i mybods
  The resulting HDF5 file using float32 format will be called 'mybods.h5'

)" << std::endl;


    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }

  // This is for testing and timing
  //
  if (vm.count("roundtrip")) {

    // Time the whole test
    ScopeTimer timer("full test");

    // Filenames
    std::string ascii = prefix + "." + suffix;
    std::string h5_32 = prefix + "_f32.h5";
    std::string h5_64 = prefix + "_f64.h5";
    std::string rest  = prefix + "_restored." + suffix;

    try {
      // Convert to HDF5 with float32 precision
      std::cout << "=== Converting ASCII to HDF5 (float32) ===" << std::endl;
      {
	ScopeTimer timer("float 32 conversion");
	ascii_to_hdf5(ascii, h5_32, FloatPrecision::FLOAT32);
      }

      // Convert to HDF5 with float64 precision
      std::cout << "\n=== Converting ASCII to HDF5 (float64) ===" << std::endl;
      {
	ScopeTimer timer("float 64 conversion");
	ascii_to_hdf5(ascii, h5_64, FloatPrecision::FLOAT64);
      }
      
      // Round-trip with float64 version
      std::cout << "\n=== Converting HDF5 (float64) back to ASCII ===" << std::endl;
      {
	ScopeTimer timer("hdf5 to ascii");
	hdf5_to_ascii(h5_64, rest);
      }
      
      std::cout << "\n=== Round-trip conversion complete! ===" << std::endl;
      
      if (vm.count("verify")) {
	std::cout << "\n=== Verifying restored ASCII file against original ===" << std::endl;
	std::ifstream original(ascii);
	std::ifstream restored(rest );

	if (!original.is_open() || !restored.is_open()) {
	  throw std::runtime_error("Could not open files for verification");
	}

	std::string line_orig, line_rest;
	size_t line_num = 0;
	bool mismatch_found = false;
	
	std::vector<double> max_diff(7, 0.0); // For m, x, y, z, u, v, w
	std::vector<double> vec_orig(7), vec_rest(7);

	while (std::getline(original, line_orig) && std::getline(restored, line_rest)) {
	  ++line_num;
	  
	  std::istringstream ss_orig(line_orig), ss_rest(line_rest);
	  
	  if (line_num == 1) {
	    int total_particles_orig, total_aux_ints_orig, total_aux_floats_orig;
	    int total_particles_rest, total_aux_ints_rest, total_aux_floats_rest;
	    
	    ss_orig >> total_particles_orig >> total_aux_ints_orig >> total_aux_floats_orig;
	    ss_rest >> total_particles_rest >> total_aux_ints_rest >> total_aux_floats_rest;

	    if (total_particles_orig != total_particles_rest ||
		total_aux_ints_orig != total_aux_ints_rest ||
		total_aux_floats_orig != total_aux_floats_rest) {
	      std::cerr << "Header mismatch at line 1:\n"
			<< "Original: " << line_orig << "\n"
			<< "Restored: " << line_rest << "\n";
	      mismatch_found = true;
	      break;
	    }
	    
	    continue;
	  }

	  for (int i = 0; i < 7; ++i) {
	    ss_orig >> vec_orig[i];
	    ss_rest >> vec_rest[i];
	  }

	  for (int i = 0; i < 7; ++i) {
	    if (std::abs(vec_orig[i] - vec_rest[i]) > max_diff[i]) {
	      max_diff[i] = std::abs(vec_orig[i] - vec_rest[i]);
	    }
	    
	  }
	}
	
	std::cout << "Maximum absolute differences for core fields"
		  << std::endl << "(m, x, y, z, u, v, w): "
		  << max_diff[0] << ", "
		  << max_diff[1] << ", "
		  << max_diff[2] << ", "
		  << max_diff[3] << ", "
		  << max_diff[4] << ", "
		  << max_diff[5] << ", "
		  << max_diff[6] << std::endl;
      }

    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;
    }


    std::cout << "\n=== Compression ratios ===" << std::endl;
    
    namespace fs = std::filesystem;

    try {
      fs::path ascPath = ascii;
      fs::path f32Path = h5_32;
      fs::path f64Path = h5_64;

      // Returns size in bytes (std::uintmax_t)
      std::uintmax_t ascSize = fs::file_size(ascPath); 
      std::uintmax_t f32Size = fs::file_size(f32Path); 
      std::uintmax_t f64Size = fs::file_size(f64Path); 

      // Print compression ratios
      std::cout << "Ascii/hdf5(float)  = "
		<< std::setprecision(3)
		<< static_cast<double>(ascSize)/f32Size
		<< std::endl;

      std::cout << "Ascii/hdf5(double) = "
		<< std::setprecision(3)
		<< static_cast<double>(ascSize)/f64Size
		<< std::endl;
      std::cout << std::endl;
    }
    catch (const fs::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }

    return 0;
  }

  // Production mode: convert ASCII to HDF5 or vice versa based on
  // mode selection
  
if (vm.count("to_hdf5")) {
    ScopeTimer timer("conversion to hdf5");
    std::string ascii = prefix + "." + suffix;
    std::string ofile = prefix + ".h5";
    if (vm.count("output")) ofile = output + ".h5";
    if (vm.count("double"))
      ascii_to_hdf5(ascii, ofile, FloatPrecision::FLOAT64);
    else
      ascii_to_hdf5(ascii, ofile, FloatPrecision::FLOAT32);
  } else if (vm.count("to_ascii")) {
    ScopeTimer timer("conversion to ascii");
    std::string ifile = prefix + ".h5";
    std::string ofile = prefix + "." + suffix;
    hdf5_to_ascii(ifile, ofile);
  } else {
    std::cerr << "No conversion mode specified. Use --ascii_to_hdf5 or --hdf5_to_ascii." << std::endl;
    return 1;
  }

  return 0;
}

