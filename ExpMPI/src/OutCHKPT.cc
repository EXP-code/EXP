#include <iostream>
#include <fstream>
#include <iomanip>
#include <strstream>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutCHKPT.H>

OutCHKPT::OutCHKPT(string& line) : Output(line)
{
  initialize();
}

void OutCHKPT::initialize()
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


void OutCHKPT::Run(int n, bool last)
{
  if (n % nint && !last) return;

  ofstream *out;
  list<Component*>::iterator cc;
  Component* c;

  if (myid==0) {
    
				// Attempt to move file to backup
    ostrstream msg;
    msg << "mv " << filename << " " << filename << ".bak" << '\0';
    system(msg.str());

				// Open file and write master header
    out = new ofstream(filename.c_str());

    if (!*out) {
      cerr << "OutCHKPT: can't open file <" << filename.c_str() 
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

}

