#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

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
    filename = outdir + "OUT." + runtag + ".chkpt";
  }

  if (Output::get_value(string("nint"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 100;

}


void OutCHKPT::Run(int n, bool last)
{
  if (n % nint && !last) return;

  synchronize_velocity(1);

  ofstream *out;
  list<Component*>::iterator cc;
  Component* c;

  if (myid==0) {
    
				// Attempt to move file to backup
    ostringstream msg;
    msg << "mv " << filename << " " << filename << ".bak";
    system(msg.str().c_str());

				// Open file and write master header
    out = new ofstream(filename.c_str());

    if (!*out) {
      cerr << "OutCHKPT: can't open file <" << filename.c_str() 
	   << "> . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 33);
    }
				// Open file and write master header
    
    struct MasterHeader header;
    header.time = tpos;
    header.ntot = comp.ntot;
    header.ncomp = comp.ncomp;

    out->write((char *)&header, sizeof(MasterHeader));
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

