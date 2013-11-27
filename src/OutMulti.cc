#include <iostream>
#include <fstream>
#include <iomanip>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutMulti.H>

OutMulti::OutMulti(string& line) : Output(line)
{
  initialize();
}

void OutMulti::initialize()
{
  string tmp;
				// Get file name
  if (!Output::get_value(string("filename"), filename)) {
    filename.erase();
    filename = outdir + "OUT.multi." + runtag;
  }

  if (Output::get_value(string("nint"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 100;

}


void OutMulti::Run(int n, bool last)
{
  if (n % nint && !last && !dump_signal) return;
  if (restart  && n==0  && !dump_signal) return;

  ofstream *out;
  list<Component*>::iterator cc;
  Component* c;

  if (myid==0) {
  }
  
  vector<unsigned> counts(multistep+1, 0), histo(multistep+1, 0);

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    PartMapItr it = c->Particles().begin();
    for (int q=0; q<c->Number(); q++) counts[c->Part((it++)->first)->level]++;
  }

  MPI_Reduce(&counts[0], &histo[0], multistep+1, MPI_UNSIGNED, MPI_SUM,
	     0, MPI_COMM_WORLD);

  if (myid==0) {
				// Open file and write master header
    ofstream out(filename.c_str(), ios::out | ios::app);

    if (!out) {
      cerr << "OutMulti: can't open file <" << filename.c_str() 
	   << "> . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 33);
    }

				// Dump the histogram
    for (int n=0; n<=multistep; n++)
      out << setw(18) << tnow << setw(5) << n 
	  << setw(8) << histo[n] << endl;
    out << endl;
  }

  dump_signal = 0;
}

