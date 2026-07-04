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

Example standalone compilation command for Ubuntu/Debian systems:

g++ -std=c++17 -I/usr/include/hdf5/serial hdf5bods.cc -o hdf5bods -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl -lhdf5 -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/serial

It is currently build in the utils/PhaseSpace directory of the EXP source tree.

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

This schema is clean, to the point, and specific to EXP while vaguely
echoing the Gadget-style.  Also, users can easily write and read this with
h5py and we can easily update our EXP-specific IC generators to use this
style by porting over the ascii_to_hdf5 routine().

*/

// C++ std
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <stdexcept>
#include <variant>
#include <memory>
#include <typeinfo>

// For OpenMP parallelization
#include <omp.h>

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

// Helper: Parse a single line of particle data
template<typename T>
void parse_particle_line(const std::string& line,
                         int particle_idx,
                         int num_aux_ints,
                         int num_aux_floats,
                         ParticleDataTemplate<T>& data)
{
  std::istringstream iss(line);
  T val;

  // Parse core PSP fields
  if (!(iss >> data.m[particle_idx] 
            >> data.x[particle_idx] >> data.y[particle_idx] >> data.z[particle_idx]
            >> data.u[particle_idx] >> data.v[particle_idx] >> data.w[particle_idx])) {
    throw std::runtime_error("Failed to parse core fields at particle " + 
                             std::to_string(particle_idx));
  }

  // Parse auxiliary integer fields
  for (int j = 0; j < num_aux_ints; ++j) {
    if (!(iss >> data.aux_ints[j][particle_idx])) {
      throw std::runtime_error("Failed to parse aux int at particle " + 
                               std::to_string(particle_idx));
    }
  }

  // Parse auxiliary float fields
  for (int j = 0; j < num_aux_floats; ++j) {
    if (!(iss >> data.aux_floats[j][particle_idx])) {
      throw std::runtime_error("Failed to parse aux float at particle " + 
                               std::to_string(particle_idx));
    }
  }
}

// Read entire ASCII file into memory
std::vector<std::string> read_ascii_lines(const std::string& ascii_file,
                                          int& num_particles,
                                          int& num_aux_ints,
                                          int& num_aux_floats)
{
  std::ifstream infile(ascii_file);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open ASCII file: " + ascii_file);
  }

  // Read header
  if (!(infile >> num_particles >> num_aux_ints >> num_aux_floats)) {
    throw std::runtime_error("Failed to read header from ASCII file");
  }

  // Validate header
  if (num_particles <= 0 || num_aux_ints < 0 || num_aux_floats < 0) {
    throw std::runtime_error("Invalid header values");
  }

  // Read all lines
  std::vector<std::string> lines;
  lines.reserve(num_particles);
  std::string line;
  while (std::getline(infile, line)) {
    if (!line.empty() && line[0] != '#') {  // Skip empty/comment lines
      lines.push_back(line);
    }
  }
  infile.close();

  if ((int)lines.size() != num_particles) {
    throw std::runtime_error("Mismatch between header count and actual particle count");
  }

  return lines;
}

// Read ASCII and write HDF5 with specified precision
template<typename T>
void ascii_to_hdf5_impl(const std::string& ascii_file, 
                        const std::string& hdf5_file,
                        FloatPrecision precision,
			bool verbose)
{
  // Header values
  int num_particles, num_aux_ints, num_aux_floats;

  // Read all ASCII lines into memory
  if (verbose) std::cout << "Reading ASCII file..." << std::endl;
  std::vector<std::string> lines = read_ascii_lines(ascii_file, num_particles, 
                                                     num_aux_ints, num_aux_floats);

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

  // Parse particle data in parallel
  if (verbose) std::cout << "Parsing particles (parallel)..." << std::endl;

  int parse_errors = 0;
  #pragma omp parallel for schedule(dynamic, 256) reduction(+:parse_errors) \
      num_threads(omp_get_max_threads())
  for (int i = 0; i < num_particles; ++i) {
    try {
      parse_particle_line<T>(lines[i], i, num_aux_ints, num_aux_floats, data);
    } catch (const std::exception& e) {
      parse_errors += 1;
      #pragma omp critical
      {
        std::cerr << "Thread error at particle " << i << ": " << e.what() << std::endl;
      }
    }
  }

  if (parse_errors > 0) {
    throw std::runtime_error("Encountered " + std::to_string(parse_errors) +
                             " particle parse errors; aborting conversion.");
  }
  // Create HDF5 file with compression enabled (SEQUENTIAL - thread-safe)
  if (verbose) std::cout << "Writing HDF5 file..." << std::endl;
  File file(hdf5_file, File::ReadWrite | File::Create | File::Truncate);

  // Store header information and precision metadata
  file.createAttribute<int>("num_particles", DataSpace::From(num_particles))
    .write(num_particles);
  file.createAttribute<int>("num_aux_ints", DataSpace::From(num_aux_ints))
    .write(num_aux_ints);
  file.createAttribute<int>("num_aux_floats", DataSpace::From(num_aux_floats))
    .write(num_aux_floats);

  // Create a group for particle data
  Group particles_group = file.createGroup("particles");

  // Define compression filters
  DataSetCreateProps props;
  hsize_t chunk = static_cast<hsize_t>(std::max(num_particles / 10, 1024));
  if (chunk > static_cast<hsize_t>(num_particles)) chunk = static_cast<hsize_t>(num_particles);
  props.add(Chunking(std::vector<hsize_t>{chunk}));
  props.add(Shuffle());
  props.add(Deflate(4));	// This is a good compromise between
				// speed and compression ratio

  // Write datasets sequentially (HDF5 is not fully thread-safe for writing)
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

  if (verbose) {
    std::string precision_str = (precision == FloatPrecision::FLOAT32) ? "float32" : "float64";
    std::cout << "Successfully wrote " << num_particles << " particles to " 
	      << hdf5_file << " (" << precision_str << ")" << std::endl;
  }
}

// Dispatcher function - user specifies precision
void ascii_to_hdf5(const std::string& ascii_file, 
                   const std::string& hdf5_file,
                   FloatPrecision precision = FloatPrecision::FLOAT64,
		   bool verbose = true)
{
  if (precision == FloatPrecision::FLOAT32) {
    ascii_to_hdf5_impl<float>(ascii_file, hdf5_file, precision, verbose);
  } else {
    ascii_to_hdf5_impl<double>(ascii_file, hdf5_file, precision, verbose);
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

// Detect float precision by inspecting dataset type
FloatPrecision detect_precision(const DataSet& dataset)
{
  const auto datatype = dataset.getDataType();
  if (datatype.getClass() != DataTypeClass::Float) {
    throw std::runtime_error("Unexpected dataset type for precision detection (expected float)");
  }

  // Get the size in bytes
  const size_t type_size = datatype.getSize();

  // Float32 is 4 bytes, Float64 is 8 bytes
  if (type_size == 4) {
    return FloatPrecision::FLOAT32;
  } else if (type_size == 8) {
    return FloatPrecision::FLOAT64;
  }

  throw std::runtime_error("Unexpected float type size: " + std::to_string(type_size));
}


// Read HDF5 with automatic precision detection
ParticleDataVariant read_hdf5_data(const std::string& hdf5_file)
{
  File file(hdf5_file, File::ReadOnly);

  // Read metadata
  int num_particles  = file.getAttribute("num_particles").read<int>();
  int num_aux_ints   = file.getAttribute("num_aux_ints").read<int>();
  int num_aux_floats = file.getAttribute("num_aux_floats").read<int>();

  // Open particles group
  Group particles_group = file.getGroup("particles");

  // Detect precision by inspecting first dataset "m"
  DataSet m_dataset = particles_group.getDataSet("m");
  FloatPrecision precision = detect_precision(m_dataset);

  ParticleDataVariant data;
  data.num_particles  = num_particles;
  data.num_aux_ints   = num_aux_ints;
  data.num_aux_floats = num_aux_floats;
  data.precision      = precision;

  // Precision is available in data.precision; caller decides whether to print it.

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
void hdf5_to_ascii(const std::string& hdf5_file, const std::string& ascii_file,
		   bool verbose = true)
{
  ParticleDataVariant data = read_hdf5_data(hdf5_file);

  std::ofstream outfile(ascii_file);
  if (!outfile.is_open()) {
    throw std::runtime_error("Could not open ASCII file for writing: " + ascii_file);
  }

  // Write header
  //
  outfile << data.num_particles << " " << data.num_aux_ints 
          << " " << data.num_aux_floats << "\n";
  outfile.close();  // Close header write

  outfile.precision(16);  // High precision for floats

  // Lambda visitor to handle both float and double
  //
  auto write_value = [](const FloatData& variant_data, size_t index)
  {
    return std::visit([index](const auto& vec)
    { return static_cast<double>(vec[index]); }, variant_data);
  };

  // Strategy: Pre-build all output lines in parallel, then write sequentially
  std::vector<std::string> output_lines(data.num_particles);

  #pragma omp parallel for schedule(static, 256) \
      num_threads(omp_get_max_threads())
  for (int i = 0; i < (int)data.num_particles; ++i) {
    std::ostringstream oss;
    oss.precision(16);

    // Write core phase-space fields
    oss << write_value(data.m, i) << " "
        << write_value(data.x, i) << " " 
        << write_value(data.y, i) << " " 
        << write_value(data.z, i) << " "
        << write_value(data.u, i) << " " 
        << write_value(data.v, i) << " " 
        << write_value(data.w, i);

    // Write auxiliary integer fields
    for (int j = 0; j < data.num_aux_ints; ++j) {
      oss << " " << data.aux_ints[j][i];
    }

    // Write auxiliary float fields
    for (int j = 0; j < data.num_aux_floats; ++j) {
      oss << " " << write_value(data.aux_floats[j], i);
    }

    output_lines[i] = oss.str();
  }

  // Write all lines sequentially to avoid I/O contention
  outfile.open(ascii_file, std::ios::app);
  if (!outfile.is_open()) {
    throw std::runtime_error("Could not open ASCII file for appending: " + ascii_file);
  }
  for (const auto& line : output_lines) {
    outfile << line << "\n";
  }
  outfile.close();

  std::string precision_str = (data.precision == FloatPrecision::FLOAT32) 
                              ? "float32" : "float64";
  if (verbose)
    std::cout << "Successfully wrote " << data.num_particles << " particles to " 
	      << ascii_file << " (" << precision_str << ")" << std::endl;
}

// Main function with command-line parsing
//
int main(int argc, char* argv[])
{
  std::string prefix = "particles", output;
  std::string suffix = "bods";
  std::string ascii_restored = "particles_restored.txt";
  int num_threads = omp_get_max_threads();
  bool quiet = false;

  // Parse command-line arguments for input/output files and mode
  //
  cxxopts::Options options(argv[0], "ASCII \u2192 HDF5 and HDF5 \u2192 ASCII particle converter for EXP body files\n"
			            "with built-in round-trip testing and float-size selection\n");
  
  options.add_options()
    ("i,input",   "Input prefix", cxxopts::value<std::string>(prefix)->default_value("particles"))
    ("o,output",  "Output prefix (optional, otherwise input prefix is used)", cxxopts::value<std::string>(output))
    ("a,suffix",  "Input suffix", cxxopts::value<std::string>(suffix)->default_value("bods"))
    ("t,threads", "Number of OpenMP threads (default: max available)", cxxopts::value<int>(num_threads)->default_value(std::to_string(num_threads)))
    ("roundtrip", "Perform round-trip conversion (ASCII -> HDF5 -> ASCII)")
    ("verify",    "Verify that the restored ASCII file matches the original in 'roundtrip' mode")
    ("to_hdf5",   "Convert ASCII to HDF5")
    ("double",    "Use double precision for HDF5 output (float is default)")
    ("to_ascii",  "Convert HDF5 to ASCII")
    ("q,quiet",   "Suppress verbose output")
    ("h,help",    "Print usage");

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;

    // Append note and custom examples
    std::cout << R"(
The best performance is achieved using the default 'float' rather than 'double'
precision for the HDF5 file, and that should be sufficient for initial data in
practice.  However, the round-trip conversion will be exact only if double
precision is used.

Examples:

  Convert a standard EXP input body file named 'mybods.bods' to HDF5 format
    $ hdf5bods --to_hdf5 -i mybods
  The resulting HDF5 file using float32 internally will be called 'mybods.h5'

  Do the same conversion but use full double precision without being chatty
    $ hdf5bods --quiet --to_hdf5 --double -i mybods
  The resulting HDF5 file will use float64 internally

  Convert the HDF5 file back to the standard EXP ascii format:
    $ hdf5bods --to_ascii -i input

  The suffix on input can be customized.  For example, the following converts
  a standard EXP input body file named 'mybods.asc' to HDF5 format
    $ hdf5bods --to_hdf5 --suffix=asc -i mybods
  The resulting HDF5 file using float32 format will be called 'mybods.h5'

  Test round-trip conversion and verify that the restored ASCII file matches the
  original:
    $ hdf5bods --roundtrip --verify -i mybods

)" << std::endl;

    return 0;
  }

  if (vm.count("quiet")) {
    quiet = true;
  }

  if (vm.count("threads")) {
    omp_set_num_threads(num_threads);
  }
  
  if (!quiet)
    std::cout << "=== Using " << num_threads << " OpenMP threads ==="
	      << std::endl;

  // This is for testing and timing
  //
  if (vm.count("roundtrip")) {

    // Time the whole test
    std::unique_ptr<ScopeTimer> time_ptr;
    if (!quiet) time_ptr = std::make_unique<ScopeTimer>("full test");

    // Filenames
    std::string ascii = prefix + "." + suffix;
    std::string h5_32 = prefix + "_f32.h5";
    std::string h5_64 = prefix + "_f64.h5";
    std::string rest  = prefix + "_restored." + suffix;

    try {
      // Convert to HDF5 with float32 precision
      std::cout << "=== Converting ASCII to HDF5 (float32) ===" << std::endl;
      {
	std::unique_ptr<ScopeTimer> time_ptr;
	if (!quiet) time_ptr = std::make_unique<ScopeTimer>("float 32 conversion");
	ascii_to_hdf5(ascii, h5_32, FloatPrecision::FLOAT32, !quiet);
      }

      // Convert to HDF5 with float64 precision
      std::cout << "\n=== Converting ASCII to HDF5 (float64) ===" << std::endl;
      {
	std::unique_ptr<ScopeTimer> time_ptr;
	if (!quiet) time_ptr = std::make_unique<ScopeTimer>("float 64 conversion");
	ascii_to_hdf5(ascii, h5_64, FloatPrecision::FLOAT64, !quiet);
      }
      
      // Round-trip with float64 version
      std::cout << "\n=== Converting HDF5 (float64) back to ASCII ===" << std::endl;
      {
	std::unique_ptr<ScopeTimer> time_ptr;
	if (!quiet) time_ptr = std::make_unique<ScopeTimer>("hdf5 to ascii");
	hdf5_to_ascii(h5_64, rest, !quiet);
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
  //
  if (vm.count("to_hdf5")) {
    std::unique_ptr<ScopeTimer> time_ptr;
    if (!quiet) time_ptr = std::make_unique<ScopeTimer>("conversion to hdf5");
    std::string ascii = prefix + "." + suffix;
    std::string ofile = prefix + ".h5";
    if (vm.count("output")) ofile = output + ".h5";
    if (vm.count("double"))
      ascii_to_hdf5(ascii, ofile, FloatPrecision::FLOAT64, !quiet);
    else
      ascii_to_hdf5(ascii, ofile, FloatPrecision::FLOAT32, !quiet);
  } else if (vm.count("to_ascii")) {
    std::unique_ptr<ScopeTimer> time_ptr;
    if (!quiet) time_ptr = std::make_unique<ScopeTimer>("conversion to ascii");
    std::string ifile = prefix + ".h5";
    std::string ofile = prefix + "." + suffix;
    hdf5_to_ascii(ifile, ofile);
  } else {
    std::cerr << "No conversion mode specified. Use --to_hdf5 to convert to HDF5, --to_ascii to convert from HDF5 to ascii, or the --roundtrip test." << std::endl;
    return 1;
  }

  return 0;
}

