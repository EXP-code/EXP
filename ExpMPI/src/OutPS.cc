#include <iostream>
#include <fstream>
#include <iomanip>
#include <strstream>

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
    filename = "OUT\0";
  }

  if (Output::get_value(string("nint"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 100;

}


void OutPS::Run(int n, bool last)
{
  if (n % nint && !last) return;

  synchronize_velocity(1);

  ofstream *out;
  list<Component*>::iterator cc;
  Component* c;

  if (myid==0) {
				// Open file and write master header
    out = new ofstream(filename.c_str(), ios::out | ios::app);

    if (!*out) {
      cerr << "OutPS: can't open file <" << filename.c_str() 
	   << "> . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 33);
    }
				// Open file and write master header
    
    struct MasterHeader header;
    header.time = tnow;
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

