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

  norb = min<int>(tcomp->nbodies_tot, norb);
  if (nskip==0) nskip = (int)tcomp->nbodies_tot/norb;


				// Make orblist
  int ncur = nbeg;
  for (int i=0; i<norb; i++) {
    if (ncur<tcomp->nbodies_tot) orblist.push_back(ncur);
    ncur += nskip;
  }
  norb = (int)orblist.size();

				// Report to log
  if (myid==0) {
    cout << "OrbTrace: " << norb << " orbit(s) from component " << tcomp->name 
	 << "[" << tcomp->id << "]\n";
  }

  pbuf = vector<double>(6);

  if (myid==0 && norb) {
				// Try to open the first time . . .
    ofstream out(filename.c_str(), ios::out | ios::app);
    if (!out) {
      ostrstream outs;
      outs << "Process " << myid << ": can't open file <" << filename << ">\n";
      bomb(outs.str());
    }

    int npos = 1;

    out << "# " << setw(4) << npos++ << setw(20) << "Time\n";
    
    for (int i=0; i<norb; i++) {
      out << "# " << setw(4) << npos++ 
	  << setw(20) << " x[" << orblist[i] << "]\n";
      out << "# " << setw(4) << npos++ 
	  << setw(20) << " y[" << orblist[i] << "]\n";
      out << "# " << setw(4) << npos++ 
	  << setw(20) << " z[" << orblist[i] << "]\n";
      out << "# " << setw(4) << npos++ 
	  << setw(20) << " u[" << orblist[i] << "]\n";
      out << "# " << setw(4) << npos++ 
	  << setw(20) << " v[" << orblist[i] << "]\n";
      out << "# " << setw(4) << npos++ 
	  << setw(20) << " w[" << orblist[i] << "]\n";
    }
    out << "# " << endl;
  }
  
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

  if (get_value(string("nint"), tmp)) 
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
  if (n % nint && !last && !tcomp && norb) return;

  MPI_Status status;

				// Open output file
  ofstream out;
  if (myid==0) {
    out.open(filename.c_str(), ios::out | ios::app);
    if (!out) {
      cout << "OrbTrace: can't open file <" << filename << ">\n";
    }
    if (out) out << setw(15) << tnow;
  }

  int curproc = 0;
  for (int i=0; i<norb; i++) {
				// Identify the node that has the 
				// desired particle
    while (orblist[i] >= tcomp->nbodies_index[curproc]) { curproc++; }

    if (myid==curproc) {	// Copy particle to buffer
      for (int k=0; k<3; k++) 
	pbuf[k  ] = 
	  tcomp->particles[orblist[i]-tcomp->nbodies_index[curproc-1]].pos[k];
      for (int k=0; k<3; k++) 
	pbuf[k+3] = 
	  tcomp->particles[orblist[i]-tcomp->nbodies_index[curproc-1]].vel[k];
    } 

    if (curproc) {		// Get particle from nodes
      
      if (myid==curproc) {	// Send
	MPI_Send(&pbuf[0], 6, MPI_DOUBLE, 0, 71, MPI_COMM_WORLD);
      }

      if (myid==0) { 		// Receive
	MPI_Recv(&pbuf[0], 6, MPI_DOUBLE, curproc, 71, MPI_COMM_WORLD,
		 &status);
      }
    }
				// Print the particle buffer
    if (myid==0 && out) {
      for (int k=0; k<6; k++) out << setw(15) << pbuf[k];
    }
    
  }
  
  if (myid==0 && out) out << endl;
}
