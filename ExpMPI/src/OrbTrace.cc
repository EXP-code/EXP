#include <iostream>
#include <fstream>
#include <iomanip>
#include <strstream>

#include <expand.h>

#include <OrbTrace.H>

OrbTrace::OrbTrace(string& line) : Output(line)
{
  nint = 1;
  norb = 5;
  nbeg = 0;
  nskip = 0;
  filename = "ORBTRACE";
  tcomp = NULL;

  initialize();

  if (!tcomp) {
    if (myid==0) {
      cerr << "OrbTrace: no component to trace\n";
      MPI_Abort(MPI_COMM_WORLD, 112);
    }
  }

  if (nskip==0) nskip = (int)tcomp->particles.size()/norb;

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
  int ncur;

  out << "# " << setw(4) << npos++ << setw(20) << "Time";

  out << "# [" << tcomp->name << "]\n";

  ncur = nskip;
  for (int i=0; i<norb; i++) {
    out << "# " << setw(4) << npos++ << setw(20) << " x[" << ncur << "]\n";
    out << "# " << setw(4) << npos++ << setw(20) << " y[" << ncur << "]\n";
    out << "# " << setw(4) << npos++ << setw(20) << " z[" << ncur << "]\n";
    out << "# " << setw(4) << npos++ << setw(20) << " u[" << ncur << "]\n";
    out << "# " << setw(4) << npos++ << setw(20) << " v[" << ncur << "]\n";
    out << "# " << setw(4) << npos++ << setw(20) << " w[" << ncur << "]\n";
    ncur += nskip;
  }

  out << "# " << endl;
  
}

void OrbTrace::initialize()
{
  string tmp;
				// Get file name
  get_value(string("filename"), filename);
  
  if (get_value(string("norb"), tmp)) 
    norb = atoi(tmp.c_str());

  if (get_value(string("nbeg"), tmp))
    nbeg = atoi(tmp.c_str());

  if (get_value(string("nskip"), tmp))
    nskip = atoi(tmp.c_str());

  if (get_value(string("freq"), tmp)) 
    nint = atoi(tmp.c_str());

				// Search for desired component
  if (get_value(string("name"), tmp)) {
    list<Component*>::iterator cc;
    Component* c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if (!(c->name.compare(tmp))) tcomp  = c;
    }
  }
  
}

void OrbTrace::Run(int n, bool last)
{
  if (n % nint && !last && !tcomp) return;

  int nbodies = tcomp->particles.size();

  if (!nbodies) return;

  ofstream out(filename.c_str(), ios::out | ios::app);
  if (!out) {
    ostrstream out;
    out << "OrbTrace: can't open file <" << filename << ">\n";
    return;
  }

  out << setw(15) << tnow;

  int ncur = nbeg;
  for (int i=0; i<norb; i++) {
    if (ncur < nbodies)
      out 
	<< setw(15) << (tcomp->particles)[ncur].pos[0]
	<< setw(15) << (tcomp->particles)[ncur].pos[1]
	<< setw(15) << (tcomp->particles)[ncur].pos[2]
	<< setw(15) << (tcomp->particles)[ncur].vel[0]
	<< setw(15) << (tcomp->particles)[ncur].vel[1]
	<< setw(15) << (tcomp->particles)[ncur].vel[2];
    ncur += nskip;
  }
  out << endl;
  
}
