#include <iostream>
#include <fstream>
#include <iomanip>
#include <strstream>

#include <expand.h>

#include <OrbTrace.H>

OrbTrace::OrbTrace(string& line) : Output(line)
{
  norb = 5;
  nskip = 0;
  filename = "ORBTRACE";

  initialize();

  ostrstream name;
  name << filename << "." << myid << '\0';
  filename = string(name.str());

				// Try to open the first time . . .
  ofstream out(filename.c_str(), ios::out | ios::app);
  if (!out) {
    ostrstream outs;
    outs << "Process " << myid << ": can't open file <" << filename << ">\n";
    bomb(outs.str());
  }
  
  list<Component*>::iterator cc;
  Component *c;
  
  int npos = 1;
  int nbodies, nmax, ncur;

  out << "# " << setw(4) << npos << setw(20) << "Time";

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {

    c = *cc;

    out << "# [" << c->id << "]\n";

    nmax = norb;
    nbodies = c->particles.size();
    if (nskip < 1) {
      nskip = nbodies/norb;
      if (nskip < 1) {
	nmax = nbodies;
	nskip = 1;
      }
    }

    ncur = 0;
    for (int i=0; i<nmax; i++) {
      out << "# " << setw(4) << npos++ << setw(20) << " x[" << ncur << "]\n";
      out << "# " << setw(4) << npos++ << setw(20) << " y[" << ncur << "]\n";
      out << "# " << setw(4) << npos++ << setw(20) << " z[" << ncur << "]\n";
      ncur += nskip;
    }
  }
  out << endl;
  
}

void OrbTrace::initialize()
{
  string tmp;
				// Get file name
  get_value(string("filename"), filename);
  
  if (!get_value(string("norb"), tmp)) 
    norb = atoi(tmp.c_str());

  if (!get_value(string("nskip"), tmp)) 
    nskip = atoi(tmp.c_str());
}

void OrbTrace::Run(int n, bool last)
{
  if (n % nint && !last) return;

  ofstream out(filename.c_str(), ios::out | ios::app);
  if (!out) {
    ostrstream out;
    out << "OrbTrace: can't open file <" << filename << ">\n";
    return;
  }

  out << setw(15) << tnow;

  list<Component*>::iterator cc;
  Component *c;
  
  int nbodies, nskip, nmax;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {

    c = *cc;

    nmax = norb;
    nbodies = c->particles.size();
    nskip = nbodies/norb;
    if (nskip < 1) {
      nmax = nbodies;
      nskip = 1;
    }

    for (int ncur=nskip-1; ncur<nmax; ncur += nskip)
      out 
	<< setw(15) << (c->particles)[ncur].pos[0]
	<< setw(15) << (c->particles)[ncur].pos[1]
	<< setw(15) << (c->particles)[ncur].pos[2];
  }
  out << endl;
  
}
