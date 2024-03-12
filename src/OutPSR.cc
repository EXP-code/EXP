#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include "expand.H"
#include <global.H>
#include <OutputContainer.H>

#include <AxisymmetricBasis.H>
#include <OutPSR.H>

const std::set<std::string>
OutPSR::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "nbeg",
  "real4",
  "timer",
  "threads"
};

OutPSR::OutPSR(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutPSR::initialize()
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
      filename = "SPL." + runtag;
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
      nintsub_warning("OutPSR");
      nintsub = std::numeric_limits<int>::max();
#endif
    } else
      nintsub = std::numeric_limits<int>::max();

    if (Output::conf["nbeg"])
      nbeg = Output::conf["nbeg"].as<int>();
    else
      nbeg = 0;

    if (Output::conf["real4"])
      real4 = Output::conf["real4"].as<bool>();
    else
      real4 = true;

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
    if (myid==0) std::cout << "Error parsing parameters in OutPSR: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("OutPSR::initialize: error parsing YAML");
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
	fname << outdir << filename << "." << setw(5) << setfill('0') << nbeg;

	// See if we can open file
	//
	ifstream in(fname.str().c_str());

	if (!in) {
	  cout << "OutPSR: will begin with nbeg=" << nbeg << endl;
	  break;
	}
      }
    }

    // Communicate starting file to all nodes
    //
    MPI_Bcast(&nbeg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
}


void OutPSR::write_thread(void)
{
  int count = 0;

  for (auto c : comp->components) {
				// Write component header
    std::ostringstream cname;
    cname << fname << "_" << count++ << "-" << myid;
    
				// Open particle file and write
    std::string blobfile = outdir + cname.str();
    std::ofstream pout(blobfile);

    if (pout.fail()) {
      std::cerr << "[" << myid << "] OutPSR: can't open file <" << cname.str() 
		<< "> . . . quitting" << std::endl;
    } else {

      c->write_binary_particles(&pout, real4);
      
      if (pout.fail()) {
	std::cout << "OutPSR: error writing binary particles to <"
		  << blobfile << std::endl;
      }
    }
  }
}

void OutPSR::Run(int n, int mstep, bool last)
{
  if (!dump_signal and !last) {
    if (n % nint) return;
    if (restart && n==0) return;
    if (multistep>1 && mstep % nintsub !=0 ) return;
  }

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();
  
  std::ofstream out;
  std::ostringstream sfname;

  // Output name prefix
  sfname << filename << "." << setw(5) << setfill('0') << nbeg++;

  fname = sfname.str();

  // Master file name
  std::string master = outdir + fname;

  int nOK = 0;

  if (myid==0) {
				// Open file and write master header
    out.open(master);

    if (out.fail()) {
      std::cerr << "OutPSR: can't open file <" << master
		<< "> . . . quitting" << std::endl;
      nOK = 1;
    }
				// Used by OutCHKPT to not duplicate a dump
    if (not real4) lastPSR = fname;
				// Open file and write master header
    if (nOK==0) {
      struct MasterHeader header;
      header.time  = tnow;
      header.ntot  = comp->ntot;
      header.ncomp = comp->ncomp;

      out.write((char *)&header, sizeof(MasterHeader));
    }
  }
  
  MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nOK) {
    MPI_Finalize();
    exit(33);
  }

#ifdef HAVE_LIBCUDA
  for (auto c : comp->components) {
    if (use_cuda) {
      if (c->force->cudaAware() and not comp->fetched[c]) {
	comp->fetched[c] = true;
	c->CudaToParticles();
      }
    }
  }
#endif
  
  auto thrd = std::make_shared<std::thread>(&OutPSR::write_thread, this);

  output->cproc.push_back(thrd);

  chktimer.mark();

  dump_signal = 0;

  if (timer) {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> intvl = end - beg;
    if (myid==0)
      std::cout << "OutPSR [T=" << tnow << "] timing=" << intvl.count()
		<< std::endl;
  }
}

