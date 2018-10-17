#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutPS.H>

OutPS::OutPS(string& line) : Output(line)
{
  initialize();
}

void OutPS::initialize()
{
  string tmp;
				// Get file name
  if (!Output::get_value(string("filename"), filename)) {
    filename.erase();
    filename = outdir + "OUT" + runtag;
  }

  if (Output::get_value(string("nint"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 100;

  if (Output::get_value(string("timer"), tmp))
    timer = atoi(tmp.c_str()) ? true : false;
  else
    timer = false;
}


void OutPS::Run(int n, bool last)
{
  if (n % nint && !last && !dump_signal) return;
  if (restart  && n==0  && !dump_signal) return;

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();

  ofstream *out;

  psdump = n;

  if (myid==0) {
				// Open file and write master header
    out = new ofstream(filename.c_str(), ios::out | ios::app);

    if (!*out) {
      cerr << "OutPS: can't open file <" << filename.c_str() 
	   << "> . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 33);
    }

    // lastPS = filename;
				// Open file and write master header
    
    struct MasterHeader header;
    header.time = tnow;
    header.ntot = comp->ntot;
    header.ncomp = comp->ncomp;

    out->write((char *)&header, sizeof(MasterHeader));
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

