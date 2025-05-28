#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

#include "expand.H"
#include <global.H>

#include <OutHDF5.H>

const std::set<std::string>
OutHDF5::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "nbeg",
  "real4",
  "real8",
  "timer",
  "threads"
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

    if (Output::conf["timer"])
      timer = Output::conf["timer"].as<bool>();
    else
      timer = false;

    if (Output::conf["threads"])
      threads = Output::conf["threads"].as<int>();
    else
      threads = 0;
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

  // Output name prefix
  fname << filename << "." << setw(5) << setfill('0') << nbeg++
	<< "." << myid+1;

  // Master file name
  std::string path = outdir + fname.str();

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

      // Create a new group for Config
      //
      HighFive::Group config = file->createGroup("Config");
      
      int ncomp = comp->components.size();
      config.createAttribute<int>("NTYPES",
				  HighFive::DataSpace::From(ncomp)).write(ncomp);

      std::string gcommit(GIT_COMMIT), gbranch(GIT_BRANCH), gdate(COMPILE_TIME);
      config.createAttribute<std::string>("Git_commit",
					  HighFive::DataSpace::From(gcommit)).write(gcommit);

      config.createAttribute<std::string>("Git_branch",
					  HighFive::DataSpace::From(gbranch)).write(gbranch);

      config.createAttribute<std::string>("Compile_date",
					  HighFive::DataSpace::From(gdate)).write(gdate);

      // Create a new group for Header
      //
      HighFive::Group header = file->createGroup("Header");
      int dp = 1;
      header.createAttribute<int>("Flag_DoublePrecision",
				  HighFive::DataSpace::From(dp)).write(dp);

      double hubble = 1, zero = 0;
      header.createAttribute<double>("HubbleParam",
				     HighFive::DataSpace::From(hubble)).write(hubble);

      header.createAttribute<double>("Omega0",
				     HighFive::DataSpace::From(zero)).write(zero);

      header.createAttribute<double>("OmegaBaryon",
				     HighFive::DataSpace::From(zero)).write(zero);

      header.createAttribute<double>("OmegaLambda",
				     HighFive::DataSpace::From(zero)).write(zero);

      header.createAttribute<double>("Redshift",
				     HighFive::DataSpace::From(zero)).write(zero);

      std::vector<double> masses(ncomp, 0.0);
      header.createAttribute<std::vector<double>>("MassTable",
						  HighFive::DataSpace::From(masses)).write(masses);

      header.createAttribute<int>("NumFilesPerSnapshot",
				  HighFive::DataSpace::From(numprocs)).write(numprocs);
      
      std::vector<unsigned> nums(ncomp);
      {
	int n=0;
	for (auto c : comp->components) nums[n++] = c->Number();
      }
      header.createAttribute<std::vector<unsigned>>("NumPart_ThisFile",
						    HighFive::DataSpace::From(nums)).write(nums);
      {
	int n=0;
	for (auto c : comp->components) nums[n++] = c->CurTotal();
      }
      header.createAttribute<std::vector<unsigned>>("NumPart_Total",
						    HighFive::DataSpace::From(nums)).write(nums);

      header.createAttribute<double>("Time",
				     HighFive::DataSpace::From(tnow)).write(tnow);
      
      // Create a new group for Parameters
      //
      HighFive::Group params = file->createGroup("Parameters");
      
    } catch (HighFive::Exception& err) {
      std::string msg("Coefs::factory: error reading HDF5 file, ");
      throw std::runtime_error(msg + err.what());
    }
	
  } catch (HighFive::Exception& err) {
    std::cerr << "OutHDF5 [" << myid << "]: can't open file <" << path
	      << "> . . . quitting" << std::endl;
    nOK = 1;
  }

  MPI_Allreduce(0, &nOK, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

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
      c->write_HDF5<float>(pgroup, true, true);
    else
      c->write_HDF5<double>(pgroup, true, true);
  }

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

