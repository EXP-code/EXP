#include <iostream>
#include <fstream>
#include <iomanip>
#include <strstream>

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
    filename = "OUT\0";
  }

  if (Output::get_value(string("nint"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 100;

  if (Output::get_value(string("nbeg"), tmp))
    nbeg = atoi(tmp.c_str());
  else
    nbeg = 0;

}


void OutPSN::Run(int n, bool last)
{
  if (n % nint && !last) return;

  synchronize_velocity(1);

  ofstream *out;
  list<Component*>::iterator cc;
  Component* c;

  if (myid==0) {
				// Output name
    ostrstream fname;
    fname << filename << "." << setw(5) << setfill('0') << nbeg++ << '\0';

				// Open file and write master header
    out = new ofstream(fname.str(), ios::out | ios::noreplace);

    if (!*out) {
      cerr << "OutPSN: can't open file <" << fname.str() 
	   << "> . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 33);
    }
				// Open file and write master header
    
    struct MasterHeader header;
    header.time = tpos;
    header.ntot = comp.ntot;
    header.ncomp = comp.ncomp;

    out->write(&header, sizeof(MasterHeader));
  }
  
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    c->write_binary(out);
  }

  if (myid==0) {
    out->close();
    delete out;
  }

  synchronize_velocity(0);

}

