#include <type_traits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

// For unlink
#include <unistd.h>

// H5 API
#include <H5Cpp.h>

// HighFive API for HDF5
#include <highfive/highfive.hpp>
#include <highfive/eigen.hpp>

#include "expand.H"
#include <global.H>

#include <OutHDF5.H>

const std::set<std::string>
OutHDF5::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "nbeg",
  "noids",
  "real4",
  "real8",
  "timer",
  "gadget",
  "checkpt",
  "preserve",
  "H5compress",
  "H5chunk"
};

OutHDF5::OutHDF5(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutHDF5::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
				// Get file name
    if (Output::conf["filename"])
      filename = Output::conf["filename"].as<std::string>();
    else{
      filename.erase();
      filename = "H5_" + runtag;
    }
    
    if (Output::conf["nint"])
      nint = Output::conf["nint"].as<int>();
    else
      nint = 100;

    if (Output::conf["nintsub"]) {
#ifdef ALLOW_NINTSUB
      nintsub = Output::conf["nintsub"].as<int>();
      if (nintsub <= 0) nintsub = 1;
#else
      nintsub_warning("OutHDF5");
      nintsub = std::numeric_limits<int>::max();
#endif
    } else
      nintsub = std::numeric_limits<int>::max();

    if (Output::conf["nbeg"])
      nbeg = Output::conf["nbeg"].as<int>();
    else
      nbeg = 0;

    if (Output::conf["real4"]) {
      real4 = Output::conf["real4"].as<bool>();
      real8 = not real4;
    } else {
      real4 = true;
      real8 = false;
    }

    if (Output::conf["real8"]) {
      real8 = Output::conf["real8"].as<bool>();
      real4 = not real8;
    }

    // This will override user setting of real4...
    if (Output::conf["checkpt"]) {
      chkpt = true;
      real8 = true;
      real4 = false;
    }

    write_flags = H5F_ACC_TRUNC;
    if (Output::conf["preserve"]) {
      bool preserve = Output::conf["preserve"].as<bool>();
      write_flags = preserve ? H5F_ACC_EXCL : H5F_ACC_TRUNC;
    }

    if (Output::conf["noids"]) {
      ids = not Output::conf["noids"].as<bool>();
    }

    if (Output::conf["timer"])
      timer = Output::conf["timer"].as<bool>();
    else
      timer = false;

    if (Output::conf["gadget"])
      gadget4 = Output::conf["gadget"].as<bool>();
    else
      gadget4 = false;

    // Default HDF5 compression is no compression.  By default,
    // shuffle is on unless turned off manually.
    //
    int H5compress = 0, H5chunk = 0;
    bool H5shuffle= true;

    if (Output::conf["H5compress"])
      H5compress  = Output::conf["H5compress"].as<int>();

    if (Output::conf["H5shuffle"])
      H5shuffle  = Output::conf["H5shuffle"].as<bool>();

    if (Output::conf["H5chunk"])
      H5chunk  = Output::conf["H5chunk"].as<int>();

    // Register compression parameters win the Component instance
    if (H5chunk>0) 
      for (auto c : comp->components) c->setH5(H5compress, H5shuffle, H5chunk);
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutHDF5: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("OutHDF5::initialize: error parsing YAML");
  }

  // Determine last file
  // 
  if (restart && nbeg==0) {

    // Only root node looks for files
    //
    if (myid==0) {

      for (nbeg=0; nbeg<100000; nbeg++) {
	// Output name
	//
	ostringstream fname;
	fname << outdir
	      <<  filename << "_" << setw(5) << setfill('0') << nbeg
	      << ".1";

	// See if we can open file
	//
	ifstream in(fname.str().c_str());

	if (!in) {
	  cout << "OutHDF5: will begin with nbeg=" << nbeg << endl;
	  break;
	}
      }
    }

    // Communicate starting file to all nodes
    //
    MPI_Bcast(&nbeg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
}


void OutHDF5::Run(int n, int mstep, bool last)
{
  if (!dump_signal and !last) {
    if (n % nint) return;
    if (restart && n==0) return;
    if (multistep>1 && mstep % nintsub !=0) return;
  }

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();
  
  std::ofstream out;
  std::ostringstream fname;

  // Checkpoint mode
  if (chkpt) {
    // Create checkpoint filename
    fname << "checkpoint_" << runtag;
    if (numprocs>1) fname << "." << myid+1;

    // Create backup filename
    std::string currfile = outdir + fname.str();
    std::string backfile = currfile + ".bak";

    // Remove old backup file
    if (unlink(backfile.c_str())) {
      if (VERBOSE>5) perror("OutHDF5::Run()");
      std::cout << "OutHDF5::Run(): error unlinking old backup file <" 
		<< backfile << ">, it may not exist" << std::endl;
    } else {
      if (VERBOSE>5) {
	std::cout << "OutHDF5::Run(): successfully unlinked <"
		  << backfile << ">" << std::endl;
      }
    }

    // Rename current file to backup file
    if (rename(currfile.c_str(), backfile.c_str())) {
      if (VERBOSE>5) perror("OutHDF5::Run()");
      std::cout << "OutHDF5: renaming backup file <" 
		<< backfile << ">, it may not exist" << std::endl;
    } else {
      if (VERBOSE>5) {
	std::cout << "OutHDF5::Run(): successfully renamed <"
		  << currfile << "> to <" << backfile << ">" << std::endl;
      }
    }
  }
  // Standard snapshot mode
  else {
    fname << filename << "_" << setw(5) << setfill('0') << nbeg++;
    if (numprocs>1) fname << "." << myid+1;
  }

  // Full path
  std::string path = outdir + fname.str();

  if (gadget4) RunGadget4(path);
  else         RunPSP(path);

  chktimer.mark();

  dump_signal = 0;

  if (timer) {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> intvl = end - beg;
    if (myid==0)
      std::cout << "OutHDF5 [T=" << tnow << "] timing=" << intvl.count()
		<< std::endl;
  }
}

void OutHDF5::RunGadget4(const std::string& path)
{
  int nOK = 0;

  std::unique_ptr<HighFive::File> file;

  // Begin HDF5 file writing
  //
  try {

    // Silence the HDF5 error stack
    //
    HighFive::SilenceHDF5 quiet;
      
    // Try opening the file as HDF5
    //
    file = std::make_unique<HighFive::File>(path,
					    HighFive::File::ReadWrite |
					    HighFive::File::Create);
    try {

      // Create a new group for Header
      //
      HighFive::Group header = file->createGroup("Header");
      int dp = 1;
      header.createAttribute("Flag_DoublePrecision", dp);

      double hubble = 1, zero = 0;
      header.createAttribute("HubbleParam", hubble);

      header.createAttribute("Omega0", zero);

      header.createAttribute("OmegaBaryon", zero);

      header.createAttribute("OmegaLambda", zero);

      header.createAttribute("Redshift", zero);

      if (masses.size()==0) checkParticleMasses();
      header.createAttribute("MassTable", masses);

      header.createAttribute("NumFilesPerSnapshot", numprocs);
      
      std::vector<unsigned long> nums(masses.size());
      {
	int n=0;
	for (auto c : comp->components) nums[n++] = c->Number();
      }
      header.createAttribute("NumPart_ThisFile", nums);
      {
	int n=0;
	for (auto c : comp->components) nums[n++] = c->CurTotal();
      }
      header.createAttribute("NumPart_Total", nums);

      header.createAttribute("Time", tnow);
      
      // Create a new group for Config
      //
      HighFive::Group config = file->createGroup("Config");

      int style = 0;
      config.createAttribute("PSPstyle", style);

      int ntypes = comp->components.size();
      config.createAttribute("NTYPES", ntypes);

      config.createAttribute("DOUBLEPRECISION", dp);

      std::vector<int> Niattrib, Ndattrib;
      for (auto & c : comp->components) {
	Niattrib.push_back(c->niattrib);
	Ndattrib.push_back(c->ndattrib);
      }

      config.createAttribute("Niattrib", Niattrib);

      config.createAttribute("Ndattrib", Ndattrib);

      // Create a new group for Parameters
      //
      HighFive::Group params = file->createGroup("Parameters");
      
      std::string gcommit(GIT_COMMIT), gbranch(GIT_BRANCH), gdate(COMPILE_TIME);
      params.createAttribute<std::string>("Git_commit", gcommit);
      params.createAttribute<std::string>("Git_branch", gbranch);
      params.createAttribute<std::string>("Compile_date", gdate);

      std::vector<std::string> names, forces, configs;
      for (auto c : comp->components) {
	names.push_back(c->name);
	forces.push_back(c->id);
	YAML::Emitter out;
	out << c->fconf; // where node is your YAML::Node
	configs.push_back(out.c_str());
      }

      params.createAttribute("ComponentNames", names);
      params.createAttribute("ForceMethods", forces);
      params.createAttribute("ForceConfigurations", configs);
      
    } catch (HighFive::Exception& err) {
      std::string msg("OutHDF5: error writing HDF5 file, ");
      throw std::runtime_error(msg + err.what());
    }
	
  } catch (HighFive::Exception& err) {
    std::cerr << "OutHDF5 [" << myid << "]: can't open file <" << path
	      << "> . . . quitting" << std::endl;
    nOK = 1;
  }

  MPI_Allreduce(MPI_IN_PLACE, &nOK, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Exit on file open failure
  //
  if (nOK) {
    MPI_Finalize();
    exit(33);
  }
  
  int count = 0;
  for (auto c : comp->components) {
    
#ifdef HAVE_LIBCUDA
    if (use_cuda) {
      if (c->force->cudaAware() and not comp->fetched[c]) {
	comp->fetched[c] = true;
	c->CudaToParticles();
      }
    }
#endif

    std::ostringstream sout;
    sout << "PartType" << count++;

    auto pgroup = file->createGroup(sout.str());

    if (real4)
      c->write_HDF5<float >(pgroup, multim[count], ids);
    else
      c->write_HDF5<double>(pgroup, multim[count], ids);
  }
}


// Template to write a scalar attribute for types used here
template <typename T>
void writeScalar(H5::Group& group, const std::string& name, T value)
  {
    // Create a dataspace for the attribute (scalar in this case)
    H5::DataSpace space(H5S_SCALAR);
	
    H5::DataType type;

    // Create a int datatype
    if constexpr (std::is_same_v<T, int>) {
      type = H5::PredType::NATIVE_INT;
    }
    else if constexpr (std::is_same_v<T, unsigned long>) {
      type = H5::PredType::NATIVE_ULONG;
    }
    else if constexpr (std::is_same_v<T, float>) {
      type = H5::PredType::NATIVE_FLOAT;
    }
    else if constexpr (std::is_same_v<T, double>) {
      type = H5::PredType::NATIVE_DOUBLE;
    }
    else if constexpr (std::is_same_v<T, std::string>) {
      type = H5::StrType(0, H5T_VARIABLE);
    }
    else {
      assert(0);
    }

    // Create the attribute
    H5::Attribute attribute = group.createAttribute(name, type, space);

    // Treat string separately
    if constexpr (std::is_same_v<T, std::string>) {
      const char* strPtr = value.c_str();
      attribute.write(type, &strPtr);
    } else
      attribute.write(type, &value);
  };


// Template to write a std::vector attribute for types used here
template <typename T>
void writeVector(H5::Group& group, const std::string& name, std::vector<T>& value)
  {
    // Create a dataspace for the attribute
    hsize_t attrDims[1] = {value.size()};
    H5::DataSpace attribute_space(1, attrDims);

    H5::DataType type;

    // Create a int datatype
    if constexpr (std::is_same_v<T, int>) {
      type = H5::PredType::NATIVE_INT;
    }
    else if constexpr (std::is_same_v<T, unsigned long>) {
      type = H5::PredType::NATIVE_ULONG;
    }
    else if constexpr (std::is_same_v<T, float>) {
      type = H5::PredType::NATIVE_FLOAT;
    }
    else if constexpr (std::is_same_v<T, double>) {
      type = H5::PredType::NATIVE_DOUBLE;
    }
    else if constexpr (std::is_same_v<T, std::string>) {
      type = H5::StrType(0, H5T_VARIABLE);
    }
    else {
      assert(0);
    }

    // Create the attribute
    H5::Attribute attribute = group.createAttribute(name, type, attribute_space);

    // Treat string separately
    if constexpr (std::is_same_v<T, std::string>) {
      std::vector<const char*> strVec;
      for (auto & v : value) strVec.push_back(v.c_str());
      attribute.write(type, strVec.data());
    } else
      attribute.write(type, value.data());
  };


void OutHDF5::RunPSP(const std::string& path)
{
  int nOK = 0;

  // Begin HDF5 file writing
  //
  try {
    // Create a file
    H5::H5File file(path.c_str(), write_flags);

    try {

      // Create a new group for Header
      //
      H5::Group header = file.createGroup("Header");

      int dp = real4 ? 0 : 1;
      writeScalar(header, "Flag_DoublePrecision", dp);

      double hubble = 1, zero = 0;

      writeScalar(header, "HubbleParam", hubble);
      writeScalar(header, "Omega0",      zero  );
      writeScalar(header, "OmegaBaryon", zero  );
      writeScalar(header, "OmegaLambda", zero  );
      writeScalar(header, "Redshift",    zero  );

      if (masses.size()==0) checkParticleMasses();
      writeVector(header, "MassTable", masses);

      writeScalar(header, "NumFilesPerSnapshot", numprocs);
      
      std::vector<unsigned long> nums(masses.size());
      {
	int n=0;
	for (auto c : comp->components) nums[n++] = c->Number();
      }
      writeVector(header, "NumPart_ThisFile", nums);
      {
	int n=0;
	for (auto c : comp->components) nums[n++] = c->CurTotal();
      }
      writeVector(header, "NumPart_Total", nums);

      writeScalar(header, "Time", tnow);
      
      // Create a new group for Config
      //
      H5::Group config = file.createGroup("Config");

      int style = 1;
      writeScalar(config, "PSPstyle", style);

      int ntypes = comp->components.size();
      writeScalar(config, "NTYPES", ntypes);

      writeScalar(config, "DOUBLEPRECISION", dp);

      std::vector<int> Niattrib, Ndattrib;
      for (auto & c : comp->components) {
	Niattrib.push_back(c->niattrib);
	Ndattrib.push_back(c->ndattrib);
      }
      writeVector(config, "Niattrib", Niattrib);
      writeVector(config, "Ndattrib", Ndattrib);

      // Create a new group for Parameters
      //
      H5::Group params = file.createGroup("Parameters");
      
      std::string gcommit(GIT_COMMIT), gbranch(GIT_BRANCH), gdate(COMPILE_TIME);
      writeScalar(params, "Git_commit",   gcommit);
      writeScalar(params, "Git_branch",   gbranch);
      writeScalar(params, "Compile_data", gdate  );

      std::vector<std::string> names, forces, configs;
      for (auto c : comp->components) {
	names.push_back(c->name);
	forces.push_back(c->id);
	YAML::Emitter out;
	out << c->fconf; // where node is your YAML::Node
	configs.push_back(out.c_str());
      }

      writeVector(params, "ComponentNames", names);
      writeVector(params, "ForceMethods",  forces);
      writeVector(params, "ForceConfigurations", configs);
      
    } catch (H5::Exception& error) {
      throw std::runtime_error(std::string("OutHDF5: error writing HDF5 file ") + error.getDetailMsg());
    }

    int count = 0;
    for (auto c : comp->components) {

#ifdef HAVE_LIBCUDA
      if (use_cuda) {
	if (c->force->cudaAware() and not comp->fetched[c]) {
	  comp->fetched[c] = true;
	  c->CudaToParticles();
	}
      }
#endif

      std::ostringstream sout;
      sout << "PartType" << count++;

      H5::Group group = file.createGroup(sout.str());

      if (real4)
	c->write_H5<float>(group);
      else
	c->write_H5<double>(group);
    }
  } catch (H5::Exception& error) {
    throw std::runtime_error(std::string("OutHDF5: error writing HDF5 file ") + error.getDetailMsg());
  }
}

void OutHDF5::checkParticleMasses()
{
  for (auto c : comp->components) {
    // First bunch of particles
    int number = -1;
    PartPtr *p = c->get_particles(&number);
    
    double minMass = std::numeric_limits<double>::max();
    double maxMass = 0.0;

    // Keep going...
    while (p) {
      for (int k=0; k<number; k++) {
	auto &P = *p;
	minMass = std::min(P->mass, minMass);
	maxMass = std::max(P->mass, maxMass);
      }
      // Next bunch of particles
      p = c->get_particles(&number);
    }

    MPI_Bcast(&minMass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxMass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if ( (maxMass - minMass)/maxMass < 1.0e-12) {
      masses.push_back(maxMass);
      multim.push_back(true);
    } else {
      masses.push_back(0.0);
      multim.push_back(false);
    }
  }
}
