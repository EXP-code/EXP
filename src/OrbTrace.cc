#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <expand.H>

#include <OrbTrace.H>

const std::set<std::string>
OrbTrace::valid_keys = {
  "filename",
  "norb",
  "nbeg",
  "nskip",
  "nint",
  "nintsub",
  "orbitlist",
  "use_acc",
  "use_pot",
  "use_lev",
  "local",
  "name"
};

OrbTrace::OrbTrace(const YAML::Node& conf) : Output(conf)
{
  nint    = 1;
  nintsub = std::numeric_limits<int>::max();
  norb    = 5;
  nbeg    = 1;
  nskip   = 0;
  use_acc = false;
  use_pot = false;
  use_lev = false;
  local   = false;

  filename = outdir + "ORBTRACE." + runtag;
  orbitlist = "";
  tcomp = NULL;

  initialize();

  // Sanity check
  //
  if (nintsub <= 0) nintsub = 1;

  if (!tcomp) {
    throw GenericError("OrbTrace: no component to trace", __FILE__, __LINE__,
		       1035, false);
  }

  if (local)
    flags = Component::Local;
  else
    flags = Component::Inertial;

  norb = min<int>(tcomp->nbodies_tot, norb);
  if (nskip==0) nskip = (int)tcomp->nbodies_tot/norb;

  if (orbitlist.size()>0) {
				// Read in orbit list
    ifstream iorb(orbitlist.c_str());
    if (!iorb) {
      throw GenericError("OrbTrace: provided orbitlist file cannot be opened",
			 __FILE__, __LINE__, 1035, false);
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
      if (ncur<=tcomp->nbodies_tot) orblist.push_back(ncur);
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
  if (use_lev) nbuf += 1;
  

  pbuf = vector<double>(nbuf);

  if (myid==0 && norb) {

    if (restart) {
      
      // Backup up old file
      string backupfile = filename + ".bak";
      if (rename(filename.c_str(), backupfile.c_str())) {
	perror("OrbTrace");
	ostringstream message;
	message << "OutTrace: error creating backup file <" 
		<< backupfile << ">";
	throw GenericError(message.str(), __FILE__, __LINE__, 1035, true);
      }
	
      // Open new output stream for writing
      ofstream out(filename.c_str());
      if (!out) {
	ostringstream message;
	message << "OrbTrace: error opening new trace file <" 
		<< filename << "> for writing";
	throw GenericError(message.str(), __FILE__, __LINE__, 1035, true);
      }
	  
      // Open old file for reading
      ifstream in(backupfile.c_str());
      if (!in) {
	ostringstream message;
	message << "OrbTrace: error opening original trace file <" 
		<< backupfile << "> for reading";
	throw GenericError(message.str(), __FILE__, __LINE__, 1035, true);
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
      std::ofstream out(filename.c_str(), ios::out | ios::app);
      if (!out) {
	std::ostringstream outs;
	outs << "Process " << myid << ": can't open file <" << filename << ">\n";
	throw GenericError(outs.str(), __FILE__, __LINE__, 1035, true);
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
	if (use_lev) {
	  out << "# " << setw(4) << npos++ 
	      << setw(20) << " lev[" << orblist[i] << "]\n";
	}
      }
      out << "# " << endl;
    }
  }
  
}

void OrbTrace::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try { 
    if (conf["filename"])    filename  = conf["filename"].as<std::string>();
    if (conf["norb"])        norb      = conf["norb"].as<int>();
    if (conf["nbeg"])        nbeg      = conf["nbeg"].as<int>();
    if (conf["nskip"])       nskip     = conf["nskip"].as<int>();
    if (conf["nint"])        nint      = conf["nint"].as<int>();
    if (conf["nintsub"])     nintsub   = conf["nintsub"].as<int>();
    if (conf["orbitlist"])   orbitlist = conf["orbitlist"].as<std::string>();
    if (conf["use_acc"])     use_acc   = conf["use_acc"].as<bool>();
    if (conf["use_pot"])     use_pot   = conf["use_pot"].as<bool>();
    if (conf["use_lev"])     use_lev   = conf["use_lev"].as<bool>();
    if (conf["local"])       local     = conf["local"].as<bool>();
    
				// Sanity check
    if (nintsub <= 0) nintsub = 1;

				// Search for desired component
    if (conf["name"]) {
      std::string tmp = conf["name"].as<std::string>();
      for (auto c : comp->components) {
	if (!(c->name.compare(tmp))) tcomp  = c;
      }
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OrbTrace: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}

void OrbTrace::Run(int n, int mstep, bool last)
{
  if (n % nint && !last && !tcomp && norb) return;
  if (mstep % nintsub !=0 && !tcomp && norb) return;
  if (tnow <= prev) return;

  prev = tnow;			// Record current time

  MPI_Status status;

#ifdef HAVE_LIBCUDA
  if (use_cuda) {
    if (tcomp->force->cudaAware() and not comp->fetched[tcomp]) {
      comp->fetched[tcomp] = true;
      tcomp->CudaToParticles();
    }
  }
#endif

				// Open output file
  std::ofstream out;
  if (myid==0) {
    out.open(filename.c_str(), ios::out | ios::app);
    if (!out) {
      std::cout << "OrbTrace: can't open file <" << filename
		<< ">" << std::endl;
    }
    if (out) out << std::setw(15) << tnow;
  }

  int curproc, nflag;
  PartMapItr it; 

  for (int i=0; i<norb; i++) {

    nflag = -1;
    if ( (it=tcomp->particles.find(orblist[i])) != tcomp->particles.end())
      nflag = myid;

    MPI_Allreduce(&nflag, &curproc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if (curproc == myid) {
      
      // Copy particle to buffer
      int icnt=0;
      for (int k=0; k<3; k++) pbuf[icnt++] = tcomp->Pos(orblist[i], k, flags);
      for (int k=0; k<3; k++) pbuf[icnt++] = tcomp->Vel(orblist[i], k, flags);
      if (use_acc) {
	for (int k=0; k<3; k++) pbuf[icnt++] = tcomp->Acc(orblist[i], k, flags);
      }
      if (use_pot) {
	pbuf[icnt++] = it->second->pot + it->second->potext;
      }
      if (use_lev) {
	pbuf[icnt++] = it->second->level;
      }

#ifdef DEBUG
      std::cout << "Process " << myid << ": packing particle #" << orblist[i]
	   << "  index=" << it->second->indx;
      for (int k=0; k<3; k++) 
	std::cout << " " << 
	  it->second->pos[k];
      std::cout << std::endl;
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
  
  if (myid==0 && out) out << std::endl;
}
