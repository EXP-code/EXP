#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutPS.H>

OutPS::OutPS(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutPS::initialize()
{
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


void OutPS::Run(int n, bool last)
{
  if (n % nint && !last && !dump_signal) return;
  if (restart  && n==0  && !dump_signal) return;

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();

  ofstream *out;

  psdump = n;

  int nOK = 0;

  if (myid==0) {
				// Open file and write master header
    out = new ofstream(filename.c_str(), ios::out | ios::app);

    if (!*out) {
      cerr << "OutPS: can't open file <" << filename.c_str() 
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
      
      out->write((char *)&header, sizeof(MasterHeader));
    }
  }
  
  MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nOK) {
    MPI_Finalize();
    exit(33);
  }

  for (auto c : comp->components) {
    c->write_binary(out, true);	// Write floats rather than doubles
  }

  if (myid==0) {
    out->close();
    delete out;
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

