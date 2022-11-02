#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "expand.H"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutPS.H>

const std::set<std::string> OutPS::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "timer"
};

OutPS::OutPS(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutPS::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("OutPS", "parameter", unmatched, __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
				// Get file name
    if (Output::conf["filename"]) {
      filename = Output::conf["filename"].as<std::string>();
    } else {
      filename.erase();
      filename = outdir + "OUT." + runtag;
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
      nintsub_warning("OutPS");
      nintsub = std::numeric_limits<int>::max();
#endif
    } else
      nintsub = std::numeric_limits<int>::max();

    if (Output::conf["timer"])
      timer = Output::conf["timer"].as<bool>();
    else
      timer = false;
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutPS: "
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


void OutPS::Run(int n, int mstep, bool last)
{
  if (!dump_signal and !last) {
    if (n % nint) return;
    if (restart && n==0) return;
    if (multistep && mstep % nintsub !=0) return;
  }

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();

  std::ofstream out;

  psdump = n;

  int nOK = 0;

  if (myid==0) {
				// Open file and write master header
    out.open(filename, ios::out | ios::app);

    if (out.fail()) {
      std::cout << "OutPS: can't open file <" << filename
		<< "> . . . quitting" << std::endl;
      nOK = 1;
    }

    // lastPS = filename;
				// Open file and write master header
    
    if (nOK==0) {
      struct MasterHeader header;
      header.time = tnow;
      header.ntot = comp->ntot;
      header.ncomp = comp->ncomp;
      
      out.write((char *)&header, sizeof(MasterHeader));
    }
  }
  
  MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nOK) {
    MPI_Finalize();
    exit(33);
  }

  for (auto c : comp->components) {
#ifdef HAVE_LIBCUDA
    if (use_cuda) {
      if (c->force->cudaAware() and not comp->fetched[c]) {
	comp->fetched[c] = true;
	c->CudaToParticles();
      }
    }
#endif
    c->write_binary(&out, true);	// Write floats rather than doubles
  }

  if (myid==0) {
    try {
      out.close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "OutPS: exception closing file <" << filename
		<< ": " << e.what() << std::endl;
    }
  }

  chktimer.mark();

  dump_signal = 0;

  if (timer) {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> intvl = end - beg;
    if (myid==0)
      std::cout << "OutPS [T=" << tnow << "] timing=" << intvl.count()
		<< std::endl;
  }
}

