// Generate slab initial conditions for varous slab models. The
// output can be in ASCII or HDF5 format, with optional compression
// filters for HDF5.

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <random>
#include <string>
#include <memory>
#include <cmath>

#include <Eigen/Dense>

// HighFive
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

// Required for H5Zfilter_avail
#include <H5Zpublic.h>

#include "massmodel1d.H"
#include "interp.H"
#include "cxxopts.H"

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

  // For most compressors, applying shuffle improves compression for
  // floating-point data.  (Blosc can do its own shuffling internally;
  // shuffle-only/checksum-only are handled below.)
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
    // SZIP is complex and highly variant; HighFive includes native
    // wrappers.  For raw testing, standard pixels-per-block option is
    // typically 16.
    cd_values = { 141, 16 }; 
    break;
    
  case 3: // Shuffle Filter alone (No trailing compression)
    props.add(HighFive::Shuffle());
    return props;
    
  case 4: // Fletcher32 Checksum alone (Data validation, not compression)
    // Fletcher32 uses built-in ID 4 and accepts 0 configuration
    // arguments.  Call the native HDF5 library directly using
    // HighFive's internal handle.
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

// Template to handle float and double precision for particle data
template<typename T>
struct ParticleData
{
  std::optional<std::vector<unsigned long>> index; // optional particle index
  std::vector<T> m;             // mass
  std::vector<T> x, y, z;       // position
  std::vector<T> u, v, w;       // velocity
  std::vector<std::vector<int>> aux_ints;   // auxiliary integer fields
  std::vector<std::vector<T>>   aux_floats; // auxiliary float fields

  size_t num_particles  = 0;
  size_t num_aux_ints   = 0;
  size_t num_aux_floats = 0;
};

// Write particle data to HDF5 file with specified filter and precision
template<typename T>
void write_hdf5_data(const std::string& outfile,
		     const ParticleData<T>& data,
		     unsigned filter_id,
		     bool verbose)
{
  // Create HDF5 file with compression enabled (SEQUENTIAL - thread-safe)
  //
  if (verbose) std::cout << "Writing HDF5 file..." << std::endl;

  std::string hdf5_file = outfile + ".h5";
  File file(hdf5_file, File::ReadWrite | File::Create | File::Truncate);

  // Store header information and precision metadata
  //
  file.createAttribute<int>("num_particles", DataSpace::From(data.num_particles))
    .write(data.num_particles);
  file.createAttribute<int>("num_aux_ints", DataSpace::From(data.num_aux_ints))
    .write(data.num_aux_ints);
  file.createAttribute<int>("num_aux_floats", DataSpace::From(data.num_aux_floats))
    .write(data.num_aux_floats);
  
  // Create a group for particle data
  //
  Group particles_group = file.createGroup("particles");
  
  // Calculate optimized, capped chunk size
  //
  hsize_t chunk_size = 262144; // Cachesafe cap (~2MB for double)
  if (data.num_particles < chunk_size) {
    chunk_size = std::max<hsize_t>(1024, data.num_particles);
  }
  if (chunk_size > data.num_particles) {
    chunk_size = data.num_particles; // Absolute safety fallback
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

  // Write optional index dataset if present
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
    std::string precision_str = "float32";
    if (sizeof(T) == 8) precision_str = "float64";
    std::cout << "Successfully wrote " << data.num_particles << " particles to " 
	      << hdf5_file << " (" << precision_str << ")" << std::endl;
  }
}

int
main(int argc, char **argv)
{
  unsigned int seed;
  int Ntable, Num_particles, Num_aux_ints, Num_aux_floats;
  double Dratio, Hratio, R, Hmax, DispX, DispZ, fJ, Lx, Ly;
  std::string outfile, config, modfile, modelType;
  unsigned filter_id = 1;
  bool Mu, HDF5 = true;

  // Parse command line
  //
  std::string message = "Generate slab initial conditions for varous slab models. The\n"
    "output can be in ASCII or HDF5 format, with optional compression\n"
    "filters for HDF5.\n";

  cxxopts::Options options(argv[0], message);

  options.add_options()
    ("h,help",     "Print this help message")
    ("5,hdf5",     "Write HDF5 output (default)")
    ("A,ascii",    "Write old-style ASCII output (default is HDF5)")
    ("v,verbose",  "Verbose output")
    ("8,double",   "Use double precision for HDF5 output (default is float)")
    ("f,filter",   "HDF5 filter ID to use (default: 1 = GZIP)", cxxopts::value<unsigned>(filter_id)->default_value("1"))
    ("N,number",   "Number of bodies",
     cxxopts::value<int>(Num_particles)->default_value("10000"))
    ("nauxint",    "Number of auxiliary integer fields",
     cxxopts::value<int>(Num_aux_ints)->default_value("0"))
    ("nauxfloat",  "Number of auxiliary float fields",
     cxxopts::value<int>(Num_aux_floats)->default_value("0"))
    ("n,ntable",   "Number of points in model table",
     cxxopts::value<int>(Ntable)->default_value("400"))
    ("m,model",    "Model type (LowIso, Sech2, Sech2mu, Sech2Halo)",
     cxxopts::value<std::string>(modelType)->default_value("Sech2"))
    ("d,dratio",    "Ratio of disk to halo density",
     cxxopts::value<double>(Dratio)->default_value("3.0"))
    ("r,hratio",    "Ratio of halo to disk scale height",
     cxxopts::value<double>(Hratio)->default_value("10.0"))
    ("R,vratio",    "Ratio of vertical to horizontal velocity dispersion",
     cxxopts::value<double>(R)->default_value("1.0"))
    ("H,hmax",      "Maximum vertical size in scale heights",
     cxxopts::value<double>(Hmax)->default_value("10.0"))
    ("x,Lx",        "Slab length in the x-direction",
     cxxopts::value<double>(Lx)->default_value("1.0"))
    ("y,Ly",        "Slab length in the y-direction",
     cxxopts::value<double>(Ly)->default_value("1.0"))
    ("X,DispX",     "In-plane velocity variance",
     cxxopts::value<double>(DispX)->default_value("1.0"))
    ("F,fJ",        "Ratio of Jeans length to box scale",
     cxxopts::value<double>(fJ)->default_value("1.0"))
    ("M,Mu",        "Surface density norm for Sech2Halo",
     cxxopts::value<bool>(Mu)->default_value("true"))
    ("c,config",    "Config file",
     cxxopts::value<std::string>(config))
    ("s,seed",      "Random # seed",
     cxxopts::value<unsigned int>(seed))
    ("i,modfile",   "Slab model file for LowIso",
     cxxopts::value<std::string>(modfile)->default_value("slab.model"))
    ("o,outfile",   "output file prefix",
     cxxopts::value<std::string>(outfile)->default_value("slab.bods"))
     ;
  
  cxxopts::ParseResult vm;

  // Parse the command line
  //
  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    return -1;
  }
  
  // Get help
  //
  if (vm.count("help")) {
    std::cout << std::endl << options.help() << std::endl;
    return 0;
  }

  // Select output format
  //
  if (vm.count("hdf5"))  HDF5 = true;
  if (vm.count("ascii")) HDF5 = false;

  // Define model
  //
  double h = 1.0;
  DispZ = DispX*R*R;
  std::shared_ptr<OneDModel> model;

  if (modelType.compare("LowIso") == 0) {
    model = std::make_shared<LowIso>(modfile);
  }
  else if (modelType.compare("Sech2mu") == 0) {
    auto tmp = std::make_shared<Sech2mu>(DispZ, DispX);
    h = tmp->get_scale_height();
    if (Hmax>0) tmp->set_hmax(Hmax);
    model = tmp;
  }
  else if (modelType.compare("Sech2") == 0) {
    auto tmp = std::make_shared<Sech2>(DispZ);
    h = tmp->get_scale_height();
    if (Hmax>0) tmp->set_hmax(Hmax);
    model = tmp;
  }
  else if (modelType.compare("Sech2Halo") == 0) {
    Sech2Halo::MU = false;
    if (Mu) Sech2Halo::MU = true;
    auto tmp = std::make_shared<Sech2Halo>(DispZ, Dratio, Hratio);
    h = tmp->get_scale_height();
    if (Hmax>0) tmp->set_hmax(Hmax);
    model = tmp;
  }
  else {
    std::cerr << "non-existent model: " << modelType << std::endl;
    exit(-1);
  }

  // Jeans' length for selecting scale height
  // 
  double maxZ = model->get_max_radius();
  double mu   = model->get_mass(maxZ);
  double KJ   = 2.0*M_PI*mu/DispX;
  double LJ   = 2.0*M_PI/KJ;

  std::cout.setf(ios::left);
  char prev = cout.fill('.');

  if (vm.count("verbose")) {
    std::cout << std::endl << "Slab model parameters:" << std::endl;
    std::cout << std::setw(40) << "Model type" << modelType << std::endl;
    std::cout << std::setw(40) << "Surface mass density" << mu << std::endl;
    std::cout << std::setw(40) << "Jeans' length" << LJ << std::endl;
    std::cout << std::setw(40) << "Scale height" << h << std::endl;
    std::cout << std::setw(40) << "Maximum thickness" << maxZ << std::endl;
  }
  
  cout.fill(prev);

  // Make mass table
  //
  std::vector<double> Z(Ntable);
  std::vector<double> M(Ntable);
  double z, dz = 2.0*maxZ/(Ntable-1.0);
  for (int i=0; i<Ntable; i++) {
    Z[i] = -maxZ + dz*i;
    M[i] = model->get_mass(Z[i]);
  }

  std::random_device rd;
  if (vm.count("seed")==0) seed = rd();
  std::mt19937 gen(seed);

  std::uniform_real_distribution<> Unit(0.0, 1.0);
  std::normal_distribution Vv{0.0, sqrt(DispZ)};
  std::normal_distribution Vh{0.0, sqrt(DispX)};

  double KE = 0.0;
  double VC = 0.0;
  double mass = mu/Num_particles;

  if (vm.count("verbose")) {
    std::cout << std::endl << "Generating " << Num_particles << " particles..." << std::endl;
  }

  if (HDF5) {

    // A lambda function to handle both float and double precision for
    // particle generation and HDF5 output.  The lambda captures the
    // necessary variables from the surrounding scope and the type of
    // the value passed to it determines the precision of the particle
    // data.  This would be easier in C++20 with concepts, but we can
    // use a simple lambda to avoid code duplication.
    //
    auto create_and_write = [&](const auto& value) {
      
      using T = std::decay_t<decltype(value)>;
        
      ParticleData<T> data;

      data.m.resize(Num_particles);
      data.x.resize(Num_particles);
      data.y.resize(Num_particles);
      data.z.resize(Num_particles);
      data.u.resize(Num_particles);
      data.v.resize(Num_particles);
      data.w.resize(Num_particles);

      data.num_particles  = Num_particles;
      data.num_aux_ints   = Num_aux_ints;
      data.num_aux_floats = Num_aux_floats;
      
      data.aux_ints.resize(Num_aux_ints);
      for (auto& vec : data.aux_ints) {
	vec.resize(Num_particles);
	vec.assign(Num_particles, 0); // Initialize auxiliary integer
				      // fields to zero
      }
      
      data.aux_floats.resize(Num_aux_floats);
      for (auto& vec : data.aux_floats) {
	vec.resize(Num_particles);
	vec.assign(Num_particles, 0.0f); // Initialize auxiliary float
					 // fields to zero
      }
    
      // This is fast enough that multithreading is not strictly
      // necessary, but we can use OpenMP to parallelize the particle
      // generation loop.  The reduction clause is used to accumulate
      // KE and VC across threads.
#pragma omp parallel for schedule(dynamic, 256) reduction(+:KE, VC) num_threads(omp_get_max_threads())
      for (int n=0; n<Num_particles; n++) {

	data.m[n] = mass;
	data.x[n] = Lx*Unit(gen);
	data.y[n] = Ly*Unit(gen);
	data.z[n] = odd2(Unit(gen)*M.back(), M, Z);
	data.u[n] = Vh(gen);
	data.v[n] = Vh(gen);
	data.w[n] = Vv(gen);
      
	KE += data.w[n] * data.w[n];
	VC += data.z[n] * model->get_dpot(data.z[n]);
      }

      write_hdf5_data<T>(outfile, data, filter_id, vm.count("verbose") > 0);
    };
      
    // Write HDF5 with float32 precision.  This induces the compiler
    // to build the template for float and double precision, allowing
    // the user to specifies --double, we instantiate with double type
    // instead.
    if (vm.count("double")) {
      double precision = 1.0; // Placeholder to indicate double precision
      create_and_write(precision);
    } else {
      float precision = 1.0f; // Placeholder to indicate float precision
      create_and_write(precision);
    }

  }
  // Old-style ascii output for backwards compatibility.  This is not
  // recommended at this point (it storage inefficient and slower to
  // read), but is retained for users who may have scripts that use
  // original ascii format.
  else {
    
    if (vm.count("verbose")) {
      std::cout << "Writing ASCII output to " << outfile << std::endl;
    }

    std::ofstream out(outfile);
    if (!out) {
      std::cerr << "Can't open <" << outfile << ">" << std::endl;
      exit(-1);
    }
    out.precision(6);
    out.setf(ios::scientific);
    
    // Header line
    out << std::setw(10) << Num_particles << std::setw(15) << 0 << std::setw(15) << 0 << std::endl;

    Eigen::MatrixXd posvel(Num_particles, 6);

    // Generate particles
#pragma omp parallel for schedule(dynamic, 256) reduction(+:KE, VC) num_threads(omp_get_max_threads())
    for (int n=0; n<Num_particles; n++) {
      posvel.row(n) <<
	Lx*Unit(gen), Ly*Unit(gen), odd2(Unit(gen)*M.back(), M, Z),
	Vh(gen), Vh(gen), Vv(gen);

      KE += posvel(n, 5)*posvel(n, 5);
      VC += posvel(n, 2)*model->get_dpot(posvel(n, 2));
    }
  
    std::cout << "Done generating particles" << std::endl;

    for (int n=0; n<Num_particles; n++) {
      out << std::setw(15) << mass;
      for (int j=0; j<6; j++) out << std::setw(15) << posvel(n, j);
      out << std::endl;
    }
  }

  // Both ASCII and HDF5 output will print virial parameters to stdout
  std::cout << std::endl
	    << "Virial parameters: KE=" << 0.5*mass*KE
	    << " VC=" << mass*VC
	    << " 2T/VC=" << KE/VC << std::endl;

  return 0;
}

