#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <expand.h>

#include <OrbTrace.H>

OrbTrace::OrbTrace(string& line) : Output(line)
{
  nint = 1;
  norb = 5;
  nbeg = 1;
  nskip = 0;
  use_acc = false;
  use_pot = false;
  filename = outdir + "ORBTRACE." + runtag;
  orbitlist = "";
  tcomp = NULL;

  initialize();

  if (!tcomp) {
    if (myid==0) {
      bomb("OrbTrace: no component to trace\n");
    }
  }

  norb = min<int>(tcomp->nbodies_tot, norb);
  if (nskip==0) nskip = (int)tcomp->nbodies_tot/norb;

  if (orbitlist.size()>0) {
				// Read in orbit list
    ifstream iorb(orbitlist.c_str());
    if (!iorb) {
      if (myid==0) {
	bomb("OrbTrace: provided orbitlist file cannot be opened\n");
      }
    }

    int cntr;
    while (1) {
      iorb >> cntr;
      if (iorb.eof() || iorb.bad()) break;
      orblist.push_back(cntr);
    }

  } else {
				// Make orblist
    int ncur = nbeg;
    for (int i=0; i<norb; i++) {
      if (ncur<tcomp->nbodies_tot) orblist.push_back(ncur);
      ncur += nskip;
    }
  }

  norb = (int)orblist.size();

				// Report to log
  if (myid==0) {
    cout << "OrbTrace: " << norb << " orbit(s) from component " << tcomp->name 
	 << "[" << tcomp->id << "]\n";
  }

  nbuf = 6;
  if (use_acc) nbuf += 3;
  if (use_pot) nbuf += 1;
  

  pbuf = vector<double>(nbuf);

  if (myid==0 && norb) {

    if (restart) {
      
      // Backup up old file
      string backupfile = filename + ".bak";
      string command("cp ");
      command += filename + " " + backupfile;
      system(command.c_str());

      // Open new output stream for writing
      ofstream out(filename.c_str());
      if (!out) {
	ostringstream message;
	message << "OrbTrace: error opening new trace file <" 
		<< filename << "> for writing\n";
	bomb(message.str());
      }
	  
      // Open old file for reading
      ifstream in(backupfile.c_str());
      if (!in) {
	ostringstream message;
	message << "OrbTrace: error opening original trace file <" 
		<< backupfile << "> for reading\n";
	bomb(message.str());
      }

      const int cbufsiz = 16384;
      char *cbuffer = new char [cbufsiz];
      double ttim;

				// Dump header
      while (in) {
	  in.getline(cbuffer, cbufsiz);
	  if (!in) break;
	  string line(cbuffer);
	  out << cbuffer << "\n";
	  if (line.find_first_of("#") == string::npos) break;
	}
	
	while (in) {
	  string line(cbuffer);
	  
	  StringTok<string> toks(line);
	  ttim  = atof(toks(" ").c_str());
	  if (tnow < ttim) break;
	  out << cbuffer << "\n";

	  in.getline(cbuffer, cbufsiz);
	  if (!in) break;
	}
	
	in.close();


    } else {
				// Try to open the first time . . .
      ofstream out(filename.c_str(), ios::out | ios::app);
      if (!out) {
	ostringstream outs;
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
	if (use_acc) {
	  out << "# " << setw(4) << npos++ 
	      << setw(20) << " ax[" << orblist[i] << "]\n";
	  out << "# " << setw(4) << npos++ 
	      << setw(20) << " ay[" << orblist[i] << "]\n";
	  out << "# " << setw(4) << npos++ 
	      << setw(20) << " az[" << orblist[i] << "]\n";
	}
	if (use_pot) {
	  out << "# " << setw(4) << npos++ 
	      << setw(20) << " pot[" << orblist[i] << "]\n";
	}
      }
      out << "# " << endl;
    }
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

  if (get_value(string("orbitlist"), tmp))
    orbitlist = tmp;

  if (get_value(string("use_acc"), tmp)) {
    if (atoi(tmp.c_str())) use_acc = true;
    else                   use_acc = false;
  }

  if (get_value(string("use_pot"), tmp)) {
    if (atoi(tmp.c_str())) use_pot = true;
    else                   use_pot = false;
  }

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

  int curproc, nflag;
  map<unsigned long, Particle>::iterator it; 

  for (int i=0; i<norb; i++) {

    nflag = -1;
    if ( (it=tcomp->particles.find(orblist[i])) != tcomp->particles.end())
      nflag = myid;

    MPI_Allreduce(&nflag, &curproc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if (curproc == myid) {
      
      // Copy particle to buffer
      int icnt=0;
      for (int k=0; k<3; k++) pbuf[icnt++] = it->second.pos[k];
      for (int k=0; k<3; k++) pbuf[icnt++] = it->second.vel[k];
      if (use_acc) {
	for (int k=0; k<3; k++) pbuf[icnt++] = it->second.acc[k];
      }
      if (use_acc) {
	pbuf[icnt++] = it->second.pot + it->second.potext;
      }

#ifdef DEBUG
      cout << "Process " << myid << ": packing particle #" << orblist[i]
	   << "  index=" << it->second.indx;
      for (int k=0; k<3; k++) 
	cout << " " << 
	  it->second.pos[k];
      cout << endl;
#endif 
    } 

    if (curproc) {		// Get particle from nodes
      
      if (myid==curproc) {	// Send
	MPI_Send(&pbuf[0], nbuf, MPI_DOUBLE, 0, 71, MPI_COMM_WORLD);
      }

      if (myid==0) { 		// Receive
	MPI_Recv(&pbuf[0], nbuf, MPI_DOUBLE, curproc, 71, MPI_COMM_WORLD,
		 &status);
      }
    }
				// Print the particle buffer
    if (myid==0 && out) {
      for (int k=0; k<nbuf; k++) out << setw(15) << pbuf[k];
    }
    
  }
  
  if (myid==0 && out) out << endl;
}
