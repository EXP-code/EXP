using namespace std;

#include <filesystem>
#include <sstream>
#include <cstdio>

#include "expand.H"

#include <OutLog.H>

char OutLog::lab_global[][19] = {
  "Time",
  "Mass", 
  "Bodies", 
  "R(x)", 
  "R(y)", 
  "R(z)", 
  "V(x)", 
  "V(y)", 
  "V(z)", 
  "L(x)", 
  "L(y)", 
  "L(z)", 
  "KE", 
  "PE", 
  "VC", 
  "E", 
  "2T/VC", 
  "Clock", 
  "# used"
};

char OutLog::lab_component[][20] = {
  "mass", 
  "bodies", 
  "R(x)", 
  "R(y)", 
  "R(z)", 
  "V(x)", 
  "V(y)", 
  "V(z)", 
  "L(x)", 
  "L(y)", 
  "L(z)", 
  "C(x)", 
  "C(y)", 
  "C(z)", 
  "KE", 
  "PE", 
  "VC", 
  "E", 
  "2T/VC", 
  "# used"
};



OutLog::OutLog(const YAML::Node& conf) : Output(conf)
{
  ektotxy=0.0;
  lastwtime = MPI_Wtime();
  laststep = -1;
  firstime = true;

  initialize();
}

void OutLog::initialize()
{
  try {
				// Get file name
    if (Output::conf["filename"])
      filename = Output::conf["filename"].as<std::string>();
    else {
      filename.erase();
      filename = outdir + "OUTLOG." + runtag;
    }
    
    if (Output::conf["freq"])
      nint = Output::conf["freq"].as<int>();
    else if (Output::conf["nint"])
      nint = Output::conf["nint"].as<int>();
    else
      nint = 1;

    if (Output::conf["nintsub"]) {
      nintsub = Output::conf["nintsub"].as<int>();
      if (nintsub <= 0) nintsub = 1;
    } else
      nintsub = std::numeric_limits<int>::max();

  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutLog: "
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



void OutLog::Run(int n, int mstep, bool last)
{
  // Define the stream
  std::ofstream out;

  // Field width
  const int cwid = 20;

  // Make a bigger output buffer
  //
  const int bufsize = 16384;
  char mybuffer [bufsize];
  out.rdbuf()->pubsetbuf(mybuffer, bufsize);
  
  if (myid==0) {

				// Open output stream for writing
    out.open(filename, ios::out | ios::app);

    if (out.fail()) {
      std::cout << "OutLog: failure on opening <" << filename
		<< "> for append (failbit)" << std::endl;
      out.close();
      return;
    }
  }

				// Generate header line
  if (firstime) {

    firstime = false;

    nbodies  = std::vector<int>(comp->ncomp);
    nbodies1 = std::vector<int>(comp->ncomp);

    used     = std::vector<int>(comp->ncomp);
    used1    = std::vector<int>(comp->ncomp);

    mtot     = std::vector<double>(comp->ncomp);
    mtot1    = std::vector<double>(comp->ncomp);

    com      = std::vector<dvector>(comp->ncomp);
    com1     = std::vector<dvector>(comp->ncomp);

    cov      = std::vector<dvector>(comp->ncomp);
    cov1     = std::vector<dvector>(comp->ncomp);

    angm     = std::vector<dvector>(comp->ncomp);
    angm1    = std::vector<dvector>(comp->ncomp);

    ctr = std::vector<dvector>(comp->ncomp);

    for (int i=0; i<comp->ncomp; i++) {
      com   [i] = std::vector<double>(3);
      com1  [i] = std::vector<double>(3);
      cov   [i] = std::vector<double>(3);
      cov1  [i] = std::vector<double>(3);
      angm  [i] = std::vector<double>(3);
      angm1 [i] = std::vector<double>(3);
      ctr   [i] = std::vector<double>(3);
    }

    com0      = std::vector<double>(3);
    cov0      = std::vector<double>(3);
    angmG     = std::vector<double>(3);
    angm0     = std::vector<double>(3);
    pos0      = std::vector<double>(3);
    vel0      = std::vector<double>(3);
	                         
    posL      = std::vector<double>(3);
    velL      = std::vector<double>(3);
    	                         
    comG      = std::vector<double>(3);
    covG      = std::vector<double>(3);
    
    ektot     = std::vector<double>(comp->ncomp);
    ektot1    = std::vector<double>(comp->ncomp);
    eptot     = std::vector<double>(comp->ncomp);
    eptot1    = std::vector<double>(comp->ncomp);
    eptotx    = std::vector<double>(comp->ncomp);
    eptotx1   = std::vector<double>(comp->ncomp);
    clausius  = std::vector<double>(comp->ncomp);
    clausius1 = std::vector<double>(comp->ncomp);

    if (myid==0) {

      if (restart) {

	try {
	  out.close();
	}
	catch (const ofstream::failure& e) {
	  std::cout << "OutLog:: exception closing file <" << filename
		    << " on restart: " << e.what() << std::endl;
	}

	// Backup up old file
	string backupfile = filename + ".bak";

	try {
	  std::filesystem::rename(filename, backupfile);
	} catch (const std::filesystem::filesystem_error& e) {
	  std::ostringstream message;
	  message << "OutLog::Run(): error creating backup file <" 
		  << backupfile << "> from <" << filename 
		  << ">, message: " << e.code().message();
	  bomb(message.str());
	}

	// Open new output stream for writing
	out.open(filename);

	if (!out) {
	  std::ostringstream message;
	  message << "OutLog: error opening new log file <" 
		  << filename << "> for writing";
	  bomb(message.str());
	}
	  
	// Open old file for reading
	std::ifstream in(backupfile.c_str());

	if (in.fail()) {
	  ostringstream message;
	  message << "OutLog: error opening original log file <" 
		  << backupfile << "> for reading";
	  bomb(message.str());
	}

	const int cbufsiz = 16384;
	std::shared_ptr<char> cbuffer;

	// Use this as of C++17
	// info = std::make_shared<char[]>(ninfochar+1);
	
	// C++14 workaround:
	cbuffer = std::shared_ptr<char>(new char[cbufsiz],
					std::default_delete<char[]>());

	// Null fill the buffer
	std::fill(cbuffer.get(), cbuffer.get()+cbufsiz, 0);

	double ttim;

				// Get the header
	while (in) {
	  in.getline(cbuffer.get(), cbufsiz);
	  if (!in) break;
	  std::string line(cbuffer.get());
	  out << cbuffer.get() << "\n";
	  if (line.find_first_of("Time") != string::npos) break;
	}
	
	while (in) {
	  in.getline(cbuffer.get(), cbufsiz);
	  if (!in) break;
	  string line(cbuffer.get());
	  
	  StringTok<string> toks(line);
	  ttim  = atof(toks(" ").c_str());
	  if (tnow < ttim) break;
	  out << cbuffer.get() << "\n";
	}
	
	try {
	  in.close();
	}
	catch (const ifstream::failure& e) {
	  std::cout << "OutLog:: exception closing input file <" << filename
		    << " on restart: " << e.what() << std::endl;
	}

      } else {

	string field;
				// Global stanza
	out << setfill('-') << setw(cwid) << "Global stats";
	for (int i=1; i<num_global; i++) 
	  out << "|" << setfill(' ') << setw(cwid) << " ";

      
				// Component stanzas
	for (auto c : comp->components) {
	  out << "|" << setw(cwid) << c->id.c_str();
	  for (int i=1; i<num_component; i++) 
	    out << "|" << setfill(' ') << setw(cwid) << " ";
	}
	out << endl;
    
				// Global divider
	out << setfill('-') << setw(cwid) << "-";
	for (int i=1; i<num_global; i++) 
	  out << "+" << setfill('-') << setw(cwid)  << "-";
      
				// Component dividers
	for (auto c : comp->components) {
	  for (int i=0; i<num_component; i++) 
	    out << "+" << setfill('-') << setw(cwid) << "-";
	}
	out << endl << setfill(' ');


				// Global labels
	out << setfill(' ') << setw(cwid) << lab_global[0];
	for (int i=1; i<num_global; i++) out << "|" << setw(cwid) << lab_global[i];
    
				// Component labels
	for (auto c : comp->components) {
	  for (int i=0; i<num_component; i++) {
	    string label = c->name + " " + lab_component[i];
	    if (label.size()<=cwid)
	      out << "|" << setw(cwid) << label.c_str();
	    else
	      out << "|" << label;
	  }
	}
	out << endl;

				// Global divider
	out << setfill('-') << setw(cwid) << "-";
	for (int i=1; i<num_global; i++) 
	  out << "+" << setfill('-') << setw(cwid) << "-";
	
				// Component dividers
	for (auto c : comp->components) {
	  for (int i=0; i<num_component; i++) 
	    out << "+" << setfill('-') << setw(cwid) << "-";
	}
	out << endl << setfill(' ');
	
				// Global count
	int count=0;
	{
	  ostringstream slab;
	  slab << "[" << ++count << "]";
	  out << setfill(' ') << setw(cwid) << slab.str();
	}
	for (int i=1; i<num_global; i++) {
	  ostringstream slab;
	  slab << "[" << ++count << "]";
	  out << "|" << setw(cwid) << slab.str();
	}    
				// Component count
	for (auto c : comp->components) {
	  for (int i=0; i<num_component; i++) {
	    ostringstream slab;
	    slab << "[" << ++count << "]";
	    out << "|" << setw(cwid) << slab.str();
	  }
	}
	out << endl;

				// Global divider
	out << setfill('-') << setw(cwid) << "-";
	for (int i=1; i<num_global; i++) 
	  out << "+" << setfill('-') << setw(cwid) << "-";
	
				// Component dividers
	for (auto c : comp->components) {
	  for (int i=0; i<num_component; i++) 
	    out << "+" << setfill('-') << setw(cwid) << "-";
	}
	out << endl << setfill(' ');
	
      }
    }
    
  } // END: firstime
  

  if (n % nint && !last) return;
  if (multistep>1 and mstep % nintsub !=0) return;


				// Use MPI wall clock to time step
  double wtime = 0.0;

  if (n>laststep) {
    curwtime = MPI_Wtime();
    wtime = (curwtime-lastwtime)/(n-laststep);
    lastwtime = curwtime;
    laststep = n;
  }

				// Zero out accumulators
  for (int i=0; i<comp->ncomp; i++) {

    nbodies [i] = nbodies1 [i] = 0;
    used    [i] = used1    [i] = 0;
    mtot    [i] = mtot1    [i] = 0.0;
    ektot   [i] = ektot1   [i] = 0.0;
    eptot   [i] = eptot1   [i] = 0.0;
    eptotx  [i] = eptotx1  [i] = 0.0;
    clausius[i] = clausius1[i] = 0.0;

    for (int j=0; j<3; j++) {
      com [i][j] = com1 [i][j] = 0.0;
      cov [i][j] = cov1 [i][j] = 0.0;
      angm[i][j] = angm1[i][j] = 0.0;
    }

  }

				// Global
  for (int j=0; j<3; j++) {
    comG [j] = com0 [j] = 0.0;
    covG [j] = cov0 [j] = 0.0;
    angmG[j] = angm0[j] = 0.0;
  }

				// Collect info
  unsigned ntot;
  int indx = 0;

  for (auto c : comp->components) {
  
#ifdef HAVE_LIBCUDA
    if (use_cuda) {
      if (c->force->cudaAware() and not comp->fetched[c]) {
	comp->fetched[c] = true;
	c->CudaToParticles();
      }
    }
#endif

    nbodies1[indx] = c->Number();

    PartMapItr it = c->Particles().begin();
    unsigned long i;

    for (int q=0; q<nbodies1[indx]; q++) {

      i = (it++)->first;

      if (c->freeze(i)) continue;

      Particle *p = c->Part(i);

      mtot1[indx] +=  p->mass;
      
      for (int k=0; k<3; k++) {
	pos0[k] = c->Pos(i, k, Component::Inertial);
	vel0[k] = c->Vel(i, k, Component::Inertial);
	posL[k] = c->Pos(i, k, Component::Local);
	velL[k] = c->Vel(i, k, Component::Local);
      }

      for (int k=0; k<3; k++) {
	com1[indx][k] += p->mass*posL[k];
	cov1[indx][k] += p->mass*velL[k];
	comG[k] += p->mass*pos0[k];
	covG[k] += p->mass*vel0[k];
      }

      angm1[indx][0] += p->mass*(posL[1]*velL[2] - posL[2]*velL[1]);
      angm1[indx][1] += p->mass*(posL[2]*velL[0] - posL[0]*velL[2]);
      angm1[indx][2] += p->mass*(posL[0]*velL[1] - posL[1]*velL[0]);

      angmG[0] += p->mass*(pos0[1]*vel0[2] - pos0[2]*vel0[1]);
      angmG[1] += p->mass*(pos0[2]*vel0[0] - pos0[0]*vel0[2]);
      angmG[2] += p->mass*(pos0[0]*vel0[1] - pos0[1]*vel0[0]);

      for (int k=0; k<3; k++) pos0[k] = c->Pos(i, k, Component::Centered);
      
      eptot1[indx]  += 0.5*p->mass*p->pot;
      eptotx1[indx] += p->mass*p->potext;
      for (int k=0; k<3; k++) {
	ektot1[indx]    += 0.5*p->mass*velL[k]*velL[k];
	clausius1[indx] += p->mass*posL[k]*p->acc[k];
      }
    }

    for (int k=0; k<3; k++) ctr[indx][k] = c->center[k];

    used1[indx] = c->force->Used();

    indx++;
  }
				// Send back to Process 0

  MPI_Reduce(&nbodies1[0], &nbodies[0], comp->ncomp, MPI_INT, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&mtot1[0], &mtot[0], comp->ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  for (int i=0; i<comp->ncomp; i++) {
    MPI_Reduce(&com1[i][0], &com[i][0], 3, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    MPI_Reduce(&cov1[i][0], &cov[i][0], 3, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    MPI_Reduce(&angm1[i][0], &angm[i][0], 3, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
  }

  MPI_Reduce(&comG[0], &com0[0], 3, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&covG[0], &cov0[0], 3, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&angmG[0], &angm0[0], 3, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&ektot1[0], &ektot[0], comp->ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&eptot1[0], &eptot[0], comp->ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&eptotx1[0], &eptotx[0], comp->ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&clausius1[0], &clausius[0], comp->ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&used1[0], &used[0], comp->ncomp, MPI_INT, MPI_SUM, 
	     0, MPI_COMM_WORLD);


  if (myid == 0) {

    // =============
    // Global
    // =============

				// Current time
    out << std::setw(cwid) << tnow;

    double mtot0 = 0.0;
    for (int i=0; i<comp->ncomp; i++) mtot0 += mtot[i];

				// Total mass
    out << "|" << setw(cwid) << mtot0;

				// Total number
    int nbodies0 = 0;
    for (int i=0; i<comp->ncomp; i++) nbodies0 += nbodies[i];
    out << "|" << setw(cwid) << nbodies0;

				// COM
    for (int j=0; j<3; j++)
      if (mtot0>0.0)
	out << "|" << setw(cwid) << com0[j]/mtot0;
      else
	out << "|" << setw(cwid) << 0.0;


				// COV
    for (int j=0; j<3; j++)
      if (mtot0>0.0)
	out << "|" << setw(cwid) << cov0[j]/mtot0;
      else
	out << "|" << setw(cwid) << 0.0;
	
				// Ang mom
    for (int j=0; j<3; j++)
      out << "|" << setw(cwid) << angm0[j];
    
				// KE
    double ektot0 = 0.0;
    for (int i=0; i<comp->ncomp; i++) ektot0 += ektot[i];
    out << "|" << setw(cwid) << ektot0;
      
				// PE
    double eptot0 = 0.0;
    for (int i=0; i<comp->ncomp; i++) eptot0 += eptot[i] + 0.5*eptotx[i];
    out << "|" << setw(cwid) << eptot0;
     
				// Clausius, Total, 2T/VC
    double clausius0 = 0.0;
    for (int i=0; i<comp->ncomp; i++) clausius0 += clausius[i];
    out << "|" << setw(cwid) << clausius0;
    out << "|" << setw(cwid) << ektot0 + clausius0;
    if (clausius0 != 0.0)
      out << "|" << setw(cwid) << -2.0*ektot0/clausius0;
    else
      out << "|" << setw(cwid) << 0.0;

    out << "|" << setw(cwid) << wtime;
    int usedT = 0;
    for (int i=0; i<comp->ncomp; i++) usedT += used[i];
    out << "|" << setw(cwid) << usedT;


    // =============
    // Per component
    // =============


    for (int i=0; i<comp->ncomp; i++) {

      out << "|" << setw(cwid) << mtot[i];
      out << "|" << setw(cwid) << nbodies[i];
      for (int j=0; j<3; j++)
	if (mtot[i]>0.0)
	  out << "|" << setw(cwid) << com[i][j]/mtot[i];
	else
	  out << "|" << setw(cwid) << 0.0;
      for (int j=0; j<3; j++)
	if (mtot[i]>0.0)
	  out << "|" << setw(cwid) << cov[i][j]/mtot[i];
	else
	  out << "|" << setw(cwid) << 0.0;
      for (int j=0; j<3; j++)
	out << "|" << setw(cwid) << angm[i][j];
      for (int j=0; j<3; j++)
	out << "|" << setw(cwid) << ctr[i][j];

      double vbar2=0.0;		// Kinetic energy in per component
      if (mtot[i]>0.0) {	// center of velocity frame
	for (int j=0; j<3; j++)
	  vbar2 +=  cov[i][j]*cov[i][j];
	vbar2 /=  mtot[i]*mtot[i];
      }
				// Update KE to cov frame
      if (nbodies[i]>1) ektot[i] -= 0.5*mtot[i]*vbar2;
      
      out << "|" << setw(cwid) << ektot[i];
      out << "|" << setw(cwid) << eptot[i] + eptotx[i];
      out << "|" << setw(cwid) << clausius[i];
      out << "|" << setw(cwid) << ektot[i] + clausius[i];
      if (clausius[i] != 0.0)
	out << "|" << setw(cwid) << -2.0*ektot[i]/clausius[i];
      else
	out << "|" << setw(cwid) << 0.0;
      out << "|" << setw(cwid) << used[i];
    }

    out << std::endl;

    try {
      out.close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "OutLog: exception closing file <" << filename
		<< ": " << e.what() << std::endl;
    }
  }

}

