/*
EXP ascii particle phase-space data HDF5 converter

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
├── index (optional Dataset)
├── aux_int_0, aux_int_1, ... (Datasets)
└── aux_float_0, aux_float_1, ... (Datasets)

This schema is clean, to the point, and specific to EXP while vaguely
echoing the Gadget-style.  Also, users can easily write and read this with
h5py and we can easily update our EXP-specific IC generators to use this
style by porting over the ascii_to_hdf5 routine().

*/

// C++ std
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <map>
#include <stdexcept>
#include <variant>
#include <memory>
#include <optional>
#include <typeinfo>

// For OpenMP parallelization
#include <omp.h>

// HighFive
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

// Required for H5Zfilter_avail
#include <H5Zpublic.h>

// Command-line parsing
#include "cxxopts.H"

// For command timing
#include "ScopeTimer.H"

using namespace HighFive;


// Factory designed to iterate across filter list available on Ubuntu
// 26.04 LTS and other systems with HDF5 1.14.x and HighFive 2.6.x
HighFive::DataSetCreateProps
createFilterProps(const std::vector<hsize_t>& chunk_dims, 
		  unsigned int filter_id,
		  // Scale from 1 (fast) to 9 (max compression)
		  int compression_level = 5)
{
  HighFive::DataSetCreateProps props;
  props.add(HighFive::Chunking(chunk_dims));

  // For most compressors, applying shuffle improves compression for floating-point data.
  // (Blosc can do its own shuffling internally; shuffle-only/checksum-only are handled below.)
  if (filter_id != 3 && filter_id != 4 && filter_id != 32001) {
    props.add(HighFive::Shuffle());
  }
  hid_t plist_id = props.getId();
  std::vector<unsigned int> cd_values;

  switch (filter_id) {
  case 1: // Deflate / GZIP (Built-in)
    // Expects 1 parameter: compression level (0-9)
    cd_values = { static_cast<unsigned int>(compression_level) };
    break;
    
  case 2: // SZIP (Built-in)
    // SZIP is complex and highly variant; HighFive includes native wrappers.
    // For raw testing, standard pixels-per-block option is typically 16.
    cd_values = { 141, 16 }; 
    break;
    
  case 3: // Shuffle Filter alone (No trailing compression)
    props.add(HighFive::Shuffle());
    return props;
    
  case 4: // Fletcher32 Checksum alone (Data validation, not compression)
    // Fletcher32 uses built-in ID 4 and accepts 0 configuration arguments.
    // Call the native HDF5 library directly using HighFive's internal handle.
    H5Pset_filter(plist_id, 4, H5Z_FLAG_OPTIONAL, 0, NULL);
    return props;
    
  case 307: // BZip2
            // Expects 1 parameter: Block size in 100KB steps (1-9)
    cd_values = { static_cast<unsigned int>(compression_level) };
    break;
    
  case 32004: // LZ4
    // Expects 1 parameter: internal chunk/block size. 
    // 0 falls back to the default (64KB), optimized for CPU caches.
    cd_values = { 0 }; 
    break;
    
  case 32001: // Blosc (v1 Meta-Compressor)
    // Expects 7 parameters. Slots 0-3 are reserved.
    // [4]=level(1-9), [5]=shuffle type(1=byte, 2=bit), [6]=codec code
    // Codecs: 0=blosclz, 1=lz4, 2=lz4hc, 3=snappy, 4=zlib, 5=zstd
    cd_values = { 0, 0, 0, 0, 
		  static_cast<unsigned int>(compression_level), 
		  1,  // 1 = Byte Shuffle (Highly recommended for float/double)
		  1 }; // 1 = Redirect internal processing to LZ4
    break;
    
  default:
    throw std::invalid_argument("Unknown or unsupported testing filter ID requested.");
    }
  
  // Pass the custom properties down to the HDF5 backend pipeline
  if (!cd_values.empty()) {
    H5Pset_filter(plist_id, filter_id, H5Z_FLAG_OPTIONAL, cd_values.size(), cd_values.data());
  }
  
  return props;
}

// Specify floating point precision
enum class FloatPrecision { FLOAT32, FLOAT64 };

// Type variant for flexible float handling
template<typename T>
struct ParticleDataTemplate
{
  std::optional<std::vector<unsigned long>> index; // optional particle index
  std::vector<T> m;             // mass
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
                         bool has_index,
                         int num_aux_ints,
                         int num_aux_floats,
                         ParticleDataTemplate<T>& data)
{
  std::istringstream iss(line);
  // Parse optional index field
  if (has_index) {
    unsigned long idx = 0;
    if (!(iss >> idx)) {
      throw std::runtime_error("Failed to parse index at particle " +
                               std::to_string(particle_idx));
    }
    (*data.index)[particle_idx] = idx;
  }

  // Parse core PSP fields
  if (!(iss
	>> data.m[particle_idx] 
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
  // Open ASCII file
  std::ifstream infile(ascii_file);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open ASCII file: " + ascii_file);
  }

  // Read header
  std::string header_line;
  if (std::getline(infile, header_line)) {
    std::istringstream header_stream(header_line);
    if (!(header_stream >> num_particles >> num_aux_ints >> num_aux_floats)) {
      throw std::runtime_error("Failed to parse header line: " + header_line);
    }
  }
  else {
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
			bool has_index,
			unsigned filter_id,
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
  if (has_index) data.index = std::vector<unsigned long>(num_particles);

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
      parse_particle_line<T>(lines[i], i, has_index, num_aux_ints, num_aux_floats, data);
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
  //
  Group particles_group = file.createGroup("particles");

  // Calculate optimized, capped chunk size
  //
  hsize_t chunk_size = 262144; // Cachesafe cap (~2MB for double)
  if (num_particles < chunk_size) {
    chunk_size = std::max<hsize_t>(1024, num_particles);
  }
  if (chunk_size > num_particles) {
    chunk_size = num_particles; // Absolute safety fallback
  }

  // Define compression filter
  //
  auto props = createFilterProps({chunk_size}, filter_id, 5);

  // Write datasets sequentially (HDF5 is not fully thread-safe for writing)
  particles_group.createDataSet("m", data.m, props);
  particles_group.createDataSet("x", data.x, props);
  particles_group.createDataSet("y", data.y, props);
  particles_group.createDataSet("z", data.z, props);
  particles_group.createDataSet("u", data.u, props);
  particles_group.createDataSet("v", data.v, props);
  particles_group.createDataSet("w", data.w, props);
  if (data.index.has_value()) {
    particles_group.createDataSet("index", *data.index, props);
  }

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
		   bool has_index = false,
		   unsigned filter_id = 1,
		   bool verbose = true)
{
  if (precision == FloatPrecision::FLOAT32) {
    ascii_to_hdf5_impl<float>(ascii_file, hdf5_file, precision, has_index, filter_id, verbose);
  } else {
    ascii_to_hdf5_impl<double>(ascii_file, hdf5_file, precision, has_index, filter_id, verbose);
  }
}

// Variant for flexible reading
using FloatData = std::variant<std::vector<float>, std::vector<double>>;

struct ParticleDataVariant
{
  std::optional<std::vector<unsigned long>> index;
  FloatData m, x, y, z, u, v, w;
  std::vector<std::vector<int>> aux_ints;
  std::vector<FloatData> aux_floats;

  size_t num_particles = 0;
  size_t num_aux_ints = 0;
  size_t num_aux_floats = 0;
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

  ParticleDataVariant data;
  data.num_particles  = num_particles;
  data.num_aux_ints   = num_aux_ints;
  data.num_aux_floats = num_aux_floats;

  // Lambda to read float or double based on precision
  auto read_dataset = [&particles_group](const std::string& name) -> FloatData {
    // Detect precision by inspecting dataset
    DataSet dataset = particles_group.getDataSet(name);
    FloatPrecision precision = detect_precision(dataset);

    if (precision == FloatPrecision::FLOAT32) {
      return particles_group.getDataSet(name).read<std::vector<float>>();
    } else {
      return particles_group.getDataSet(name).read<std::vector<double>>();
    }
  };

  // Read core PSP fields
  if (particles_group.exist("index")) {
    data.index = particles_group.getDataSet("index").read<std::vector<unsigned long>>();
  }
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

    // Write optional index and core phase-space fields
    if (data.index.has_value()) oss << (*data.index)[i] << " ";
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

  if (verbose)
    std::cout << "Successfully wrote " << data.num_particles
	      << " particles to " << ascii_file << std::endl;
}

// Main function with command-line parsing
//
int main(int argc, char* argv[])
{
  std::string prefix = "particles", output;
  std::string suffix = "bods";
  int num_threads = omp_get_max_threads();
  unsigned filter_id = 1;
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
    ("f,filter",  "HDF5 filter ID to use (default: 1 = GZIP)", cxxopts::value<unsigned>(filter_id)->default_value("1"))
    ("roundtrip", "Perform round-trip conversion (ASCII -> HDF5 -> ASCII)")
    ("verify",    "Verify that the restored ASCII file matches the original in 'roundtrip' mode")
    ("to_hdf5",   "Convert ASCII to HDF5")
    ("double",    "Use double precision for HDF5 output (float is default)")
    ("filterlist","List available HDF5 filters and exit")
    ("aindex",    "Input ASCII has leading particle index column; write the particle index dataset to HDF5")
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

  if (vm.count("filterlist")) {

    // Map of commonly registered official HDF5 Filter IDs
    std::map<unsigned int, std::string> known_filters = {
      {1, "Deflate / GZIP (Built-in)"},
      {2, "SZIP (Built-in)"},
      {3, "Shuffle (Built-in)"},
      {4, "Fletcher32 Checksum (Built-in)"},
      {32001, "Blosc"},
      {32004, "LZ4"},
      {32008, "Bitshuffle"},
      {32013, "Zfp"},
      {32015, "Zstd"},
      {307, "BZip2"}
    };
    
    std::cout << "Checking available HDF5 filters on your system:\n";
    std::cout << "-----------------------------------------------\n";

    for (const auto& [id, name] : known_filters) {
      // H5Zfilter_avail returns a positive number if available, 0 if not
      htri_t available = H5Zfilter_avail(id);
        
      if (available > 0) {
	std::cout << "[AVAILABLE] ID: " << id << " -> " << name << "\n";
      } else {
	std::cout << "[NOT FOUND] ID: " << id << " -> " << name << "\n";
      }
    }
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
	ascii_to_hdf5(ascii, h5_32, FloatPrecision::FLOAT32, vm.count("aindex"), filter_id, !quiet);
      }

      // Convert to HDF5 with float64 precision
      std::cout << "\n=== Converting ASCII to HDF5 (float64) ===" << std::endl;
      {
	std::unique_ptr<ScopeTimer> time_ptr;
	if (!quiet) time_ptr = std::make_unique<ScopeTimer>("float 64 conversion");
	ascii_to_hdf5(ascii, h5_64, FloatPrecision::FLOAT64, vm.count("aindex"), filter_id, !quiet);
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

	std::vector<double> max_diff(7, 0.0); // For m, x, y, z, u, v, w
	std::vector<double> vec_orig(7), vec_rest(7);
	std::vector<double> max_aux_int_diff;
	std::vector<double> max_aux_float_diff;
	std::vector<int> aux_int_orig, aux_int_rest;
	std::vector<double> aux_float_orig, aux_float_rest;
	int num_aux_ints = -1, num_aux_floats = -1;

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
	      break;
	    }

	    num_aux_ints = total_aux_ints_orig;
	    num_aux_floats = total_aux_floats_orig;
	    max_aux_int_diff.assign(num_aux_ints, 0.0);
	    max_aux_float_diff.assign(num_aux_floats, 0.0);
	    aux_int_orig.assign(num_aux_ints, 0);
	    aux_int_rest.assign(num_aux_ints, 0);
	    aux_float_orig.assign(num_aux_floats, 0.0);
	    aux_float_rest.assign(num_aux_floats, 0.0);
	    
	    continue;
	  }

	  if (vm.count("aindex")) {
	    unsigned long idx_orig = 0, idx_rest = 0;
	    if (!(ss_orig >> idx_orig) || !(ss_rest >> idx_rest)) {
	      throw std::runtime_error("Failed to parse index field at line " + std::to_string(line_num));
	    }
	  }

	  for (int i = 0; i < 7; ++i) {
	    if (!(ss_orig >> vec_orig[i]) || !(ss_rest >> vec_rest[i])) {
	      throw std::runtime_error("Failed to parse core field at line " + std::to_string(line_num));
	    }
	  }

	  for (int i = 0; i < 7; ++i) {
	    if (std::abs(vec_orig[i] - vec_rest[i]) > max_diff[i]) {
	      max_diff[i] = std::abs(vec_orig[i] - vec_rest[i]);
	    }
	    
	  }

	  for (int i = 0; i < num_aux_ints; ++i) {
	    if (!(ss_orig >> aux_int_orig[i]) || !(ss_rest >> aux_int_rest[i])) {
	      throw std::runtime_error("Failed to parse auxiliary integer field at line " + std::to_string(line_num));
	    }
	  }

	  for (int i = 0; i < num_aux_ints; ++i) {
	    auto diff_i = std::llabs(static_cast<long long>(aux_int_orig[i]) -
				     static_cast<long long>(aux_int_rest[i]));
	    double diff = static_cast<double>(diff_i);
	    if (diff > max_aux_int_diff[i]) max_aux_int_diff[i] = diff;
	  }

	  for (int i = 0; i < num_aux_floats; ++i) {
	    if (!(ss_orig >> aux_float_orig[i]) || !(ss_rest >> aux_float_rest[i])) {
	      throw std::runtime_error("Failed to parse auxiliary float field at line " + std::to_string(line_num));
	    }
	  }

	  for (int i = 0; i < num_aux_floats; ++i) {
	    double diff = std::abs(aux_float_orig[i] - aux_float_rest[i]);
	    if (diff > max_aux_float_diff[i]) max_aux_float_diff[i] = diff;
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

	std::cout << "Maximum absolute differences for auxiliary integer fields: ";
	if (max_aux_int_diff.empty()) {
	  std::cout << "(none)";
	} else {
	  for (size_t i = 0; i < max_aux_int_diff.size(); ++i) {
	    if (i) std::cout << ", ";
	    std::cout << max_aux_int_diff[i];
	  }
	}
	std::cout << std::endl;

	std::cout << "Maximum absolute differences for auxiliary float fields: ";
	if (max_aux_float_diff.empty()) {
	  std::cout << "(none)";
	} else {
	  for (size_t i = 0; i < max_aux_float_diff.size(); ++i) {
	    if (i) std::cout << ", ";
	    std::cout << max_aux_float_diff[i];
	  }
	}
	std::cout << std::endl;
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
      ascii_to_hdf5(ascii, ofile, FloatPrecision::FLOAT64, vm.count("aindex"), filter_id, !quiet);
    else
      ascii_to_hdf5(ascii, ofile, FloatPrecision::FLOAT32, vm.count("aindex"), filter_id, !quiet);
  } else if (vm.count("to_ascii")) {
    std::unique_ptr<ScopeTimer> time_ptr;
    if (!quiet) time_ptr = std::make_unique<ScopeTimer>("conversion to ascii");
    std::string ifile = prefix + ".h5";
    std::string ofile = prefix + "." + suffix;
    hdf5_to_ascii(ifile, ofile, !quiet);
  } else {
    std::cerr << "No conversion mode specified. Use --to_hdf5 to convert to HDF5, --to_ascii to convert from HDF5 to ascii, or the --roundtrip test." << std::endl;
    return 1;
  }

  return 0;
}
