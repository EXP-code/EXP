#include <filesystem>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <expand.H>

#include <OutVel.H>

const std::set<std::string>
OutVel::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "name",
  "dof",
  "scale",
  "rmin",
  "rmax",
  "ascl",
  "delta",
  "lmax",
  "nmax",
  "model"
};

OutVel::OutVel(const YAML::Node& conf) : Output(conf)
{
  nint    = 10;
  nintsub = std::numeric_limits<int>::max();
  tcomp   = NULL;
  lmax    = 4;
  nmax    = 10;
  rmin    = 1.0e-4;
  rmax    = 2.0;
  ascl    = 0.01;
  delta   = 0.005;
  scale   = 0.05;
  model   = "file";

  initialize();

  // Target output file
  //
  outfile = "velcoef." + tcomp->name;

  // Check for valid model type
  //
  std::vector<std::string> modelTypes { "file", "expon" };
  auto it = std::find(modelTypes.begin(), modelTypes.end(), model);
  if (it == modelTypes.end()) {
    std::ostringstream sout;
    sout << "OutVel: found type <" << model << ">.  Must be one of ";
    for (auto v : modelTypes) sout << " " << v;
    throw std::runtime_error(sout.str());
  }

  // Check dof value
  //
  if (dof!=2 or dof!=3) {
    std::ostringstream sout;
    sout << "OutVel: found " << dof << " for dof.  Must be 2 or 3.";
    throw std::runtime_error(sout.str());
  }

  // Create the basis
  //
  basis = std::make_shared<OrthoBasisClasses::VelocityBasis>(conf);
}

void OutVel::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["nint"])         nint     = conf["nint" ].as<int>();
    if (conf["nintsub"]) {
      nintsub  = conf["nintsub"].as<int>();
      if (nintsub <= 0) nintsub = 1;
    }
    if (conf["name"])
      {				// Search for desired component
	std::string tmp = conf["name"].as<std::string>();
	for (auto c : comp->components) {
	  if (!(c->name.compare(tmp))) tcomp  = c;
	}
      }

    if (!tcomp) {
      std::string message = "OutVel: no component to trace. Please specify "
	"the component name using the 'name' parameter.";
      throw std::runtime_error(message);
    }

    if (conf["filename"])
      {
	filename = conf["filename"].as<std::string>();
      }
    else
      {
	filename = outdir + "velcoef." + tcomp->name + "." + runtag;
      }

  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutVel: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}

void OutVel::Run(int n, int mstep, bool last)
{
  // Skip this master step
  //
  if (n % nint != 0 && !last) return;

  // Skip this sub step
  //
  if (mstep < std::numeric_limits<int>::max() and mstep % nintsub != 0) return;

  // Zero coefficients
  //
  basis->reset_coefs();
    
  // Compute coefficients
  //
  PartMapItr ibeg = tcomp->Particles().begin();
  PartMapItr iend = tcomp->Particles().end();

  constexpr std::complex<double> I(0.0, 1.0);

#pragma omp parallel
  {
    for (PartMapItr it=ibeg; it!=iend; ++it) {
      
#pragma omp single nowait
      {
	double M = it->second->mass;

	double x = it->second->pos[0];
	double y = it->second->pos[1];
	double z = it->second->pos[2];
	
	double u = it->second->vel[0];
	double v = it->second->vel[1];
	double w = it->second->vel[2];
	
	basis->accumulate(M, x, y, z, u, v, w);
      }
    }
  }
  
  // Make coefficients
  //
  basis->make_coefs();

  // Check if file exists
  //
  if (std::filesystem::exists(outfile)) {
    ExtendH5Coefs();
  }
  // Otherwise, extend the existing HDF5 file
  //
  else {
    WriteH5Coefs();
  }
}

void OutVel::WriteH5Coefs()
{
  try {
    // Create a new hdf5 file
    //
    HighFive::File file(outfile,
			HighFive::File::ReadWrite |
			HighFive::File::Create);
      
    // Write the specific parameters
    //
    WriteH5Params(file);
      
    // Group count variable
    //
    unsigned count = 0;
    HighFive::DataSet dataset = file.createDataSet("count", count);
      
    // Create a new group for coefficient snapshots
    //
    HighFive::Group group = file.createGroup("snapshots");
    
    // Write the coefficients
    //
    count = WriteH5Times(group, count);
      
    // Update the count
    //
    dataset.write(count);
    
  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
  
}
  

void OutVel::ExtendH5Coefs()
{
  try {
    // Open an hdf5 file
    //
    HighFive::File file(outfile, HighFive::File::ReadWrite);
      
    // Get the dataset
    HighFive::DataSet dataset = file.getDataSet("count");
      
    unsigned count;
    dataset.read(count);
    
    HighFive::Group group = file.getGroup("snapshots");
    
    // Write the coefficients
    //
    count = WriteH5Times(group, count);
    
    // Update the count
    //
    dataset.write(count);
    
  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
  
}
  
void OutVel::WriteH5Params(HighFive::File& file)
{
  // Identify myself
  //
  std::string whoami("Velocity orthgonal function coefficients");
  file.createAttribute<std::string>("whoami", HighFive::DataSpace::From(whoami)).write(whoami);

  // Write the model type
  //
  file.createAttribute<std::string>("model", HighFive::DataSpace::From(model)).write(model);
  
  // We write the coefficient mnemonic
  //
  file.createAttribute<std::string>("name", HighFive::DataSpace::From(tcomp->name)).write(tcomp->name);
    
  // Stash the config string
  //
  std::ostringstream yaml; yaml << conf;
  file.createAttribute<std::string>("config", HighFive::DataSpace::From(yaml.str())).write(yaml.str());
  

  // Write the remaining parameters
  //
  file.createAttribute<int>   ("lmax",  HighFive::DataSpace::From(lmax)  ).write(lmax);
  file.createAttribute<int>   ("nmax",  HighFive::DataSpace::From(nmax)  ).write(nmax);
  file.createAttribute<int>   ("dof",   HighFive::DataSpace::From(dof)   ).write(dof);
  file.createAttribute<double>("scale", HighFive::DataSpace::From(scale) ).write(scale);
  file.createAttribute<double>("rmin",  HighFive::DataSpace::From(rmin)  ).write(rmin);
  file.createAttribute<double>("rmax",  HighFive::DataSpace::From(rmax)  ).write(rmax);
  file.createAttribute<double>("ascl",  HighFive::DataSpace::From(ascl)  ).write(ascl);
  file.createAttribute<double>("delta", HighFive::DataSpace::From(delta) ).write(delta);
}
  
unsigned OutVel::WriteH5Times(HighFive::Group& snaps, unsigned count)
{
  std::ostringstream stim;
  stim << std::setw(8) << std::setfill('0') << std::right << count++;
  HighFive::Group stanza = snaps.createGroup(stim.str());
      
  // Add time attribute
  //
  stanza.createAttribute<double>("Time", HighFive::DataSpace::From(tnow)).write(tnow);
      
  // Coefficient size (allow Eigen::Tensor to be easily recontructed from metadata)
  //
  const auto& d = coefs[0]->dimensions();
  std::array<long int, 3> shape {d[0], d[1], d[2]};
  stanza.createAttribute<double>("shape", HighFive::DataSpace::From(shape)).write(shape);

  // Add coefficient data from flattened tensor
  //
  HighFive::DataSet dataset = stanza.createDataSet("coefficients", store);
  
  return count;
}
