#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutPSP.H>

OutPSP::OutPSP(string& line) : Output(line)
{
  initialize();
}

void OutPSP::initialize()
{
  string tmp;
				// Get file name
  if (!Output::get_value(string("filename"), filename)) {
    filename.erase();
    filename = outdir + "OUTS." + runtag;
  }

  if (Output::get_value(string("nint"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 100;

  if (Output::get_value(string("nbeg"), tmp))
    nbeg = atoi(tmp.c_str());
  else
    nbeg = 0;

  if (Output::get_value(string("real4"), tmp))
    real4 = atoi(tmp.c_str()) ? true : false;
  else
    real4 = false;

  if (Output::get_value(string("real8"), tmp))
    real4 = atoi(tmp.c_str()) ? false : true;
  else
    real4 = false;

  if (Output::get_value(string("timer"), tmp))
    timer = atoi(tmp.c_str()) ? true : false;
  else
    timer = false;
				// Determine last file

  if (restart && nbeg==0 && myid==0) {

    for (nbeg=0; nbeg<100000; nbeg++) {

				// Output name
      ostringstream fname;
      fname << filename << "." << setw(5) << setfill('0') << nbeg;

				// See if we can open file
      ifstream in(fname.str().c_str());

      if (!in) {
	cout << "OutPSP: will begin with nbeg=" << nbeg << endl;
	break;
      }
    }
  }
}


void OutPSP::Run(int n, bool last)
{
  if (n % nint && !last && !dump_signal) return;
  if (restart  && n==0  && !dump_signal) return;

  std::chrono::high_resolution_clock::time_point beg, end;
  if (timer) beg = std::chrono::high_resolution_clock::now();
  
  MPI_Offset offset;
  MPI_Status status;
  MPI_File   file;

  psdump = n;

  // Output name
  //
  ostringstream fname;
  fname << filename << "." << setw(5) << setfill('0') << nbeg++;

  // Open shared file and write master header
  //
  int ret = MPI_File_open(MPI_COMM_WORLD, fname.str().c_str(),
			  MPI_MODE_CREATE | MPI_MODE_WRONLY,
			  MPI_INFO_NULL, &file);

  if (ret != MPI_SUCCESS) {
    cerr << "OutPSP: can't open file <" << fname.str() 
	 << "> . . . quitting\n";
    MPI_Abort(MPI_COMM_WORLD, 33);
  }
				// Used by OutCHKPT to not duplicate a dump
  lastPS = fname.str();
				// Open file and write master header
  
  if (myid==0) {
    struct MasterHeader header;
    header.time  = tnow;
    header.ntot  = comp->ntot;
    header.ncomp = comp->ncomp;
    
    MPI_File_write_at(file, offset, (char *)&header, sizeof(MasterHeader),
		      MPI_CHAR, &status);
  }
  
  offset += sizeof(MasterHeader);

  for (auto c : comp->components) {
    c->write_binary_mpi(file, offset, real4); 
  }

  MPI_File_close(&file);

  chktimer.mark();

  dump_signal = 0;

  if (timer) {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> intvl = end - beg;
    if (myid==0)
      std::cout << "OutPSP [T=" << tnow << "] timing=" << intvl.count()
	      << std::endl;
  }
}
