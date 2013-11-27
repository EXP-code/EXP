#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutPSN.H>

OutPSN::OutPSN(string& line) : Output(line)
{
  initialize();
}

void OutPSN::initialize()
{
  string tmp;
				// Get file name
  if (!Output::get_value(string("filename"), filename)) {
    filename.erase();
    filename = outdir + "OUT." + runtag;
  }

  if (Output::get_value(string("nint"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 100;

  if (Output::get_value(string("nbeg"), tmp))
    nbeg = atoi(tmp.c_str());
  else
    nbeg = 0;

				// Determine last file

  if (restart && nbeg==0 && myid==0) {

    for (nbeg=0; nbeg<100000; nbeg++) {

				// Output name
      ostringstream fname;
      fname << filename << "." << setw(5) << setfill('0') << nbeg;

				// See if we can open file
      ifstream in(fname.str().c_str());

      if (!in) {
	cout << "OutPSN: will begin with nbeg=" << nbeg << endl;
	break;
      }
    }
  }
}


void OutPSN::Run(int n, bool last)
{
  if (n % nint && !last && !dump_signal) return;
  if (restart  && n==0  && !dump_signal) return;

  ofstream *out;
  list<Component*>::iterator cc;
  Component* c;

  psdump = n;

  if (myid==0) {
				// Output name
    ostringstream fname;
    fname << filename << "." << setw(5) << setfill('0') << nbeg++;

				// Open file and write master header
    out = new ofstream(fname.str().c_str());

    if (!*out) {
      cerr << "OutPSN: can't open file <" << fname.str() 
	   << "> . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 33);
    }
				// Used by OutCHKPT to not duplicate a dump
    lastPS = fname.str();
				// Open file and write master header
    
    struct MasterHeader header;
    header.time = tnow;
    header.ntot = comp.ntot;
    header.ncomp = comp.ncomp;

    out->write((char *)&header, sizeof(MasterHeader));
  }
  
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    c->write_binary(out, true);	// Write floats rather than doubles
  }

  if (myid==0) {
    out->close();
    delete out;
  }

  chktimer.mark();

  dump_signal = 0;
}

