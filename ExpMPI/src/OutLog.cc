static char rcsid[] = "$Id$";

using namespace std;

#include <cstdio>
#include <sstream>
#include "expand.h"

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



OutLog::OutLog(string& line) : Output(line)
{
  ektotxy=0.0;
  lastwtime = MPI_Wtime();
  laststep = -1;
  firstime = true;

  initialize();
}

void OutLog::initialize()
{
  string tmp;
				// Get file name
  if (!Output::get_value(string("filename"), filename)) {
    filename.erase();
    filename = outdir + "OUTLOG." + runtag;
  }

  if (Output::get_value(string("freq"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 1;
}



void OutLog::Run(int n, bool last)
{
  list<Component*>::iterator cc;
  Component *c;
  ofstream *out = 0;
  const int cwid = 20;

  if (myid==0) {
				// Open output stream for writing
    out = new ofstream(filename.c_str(), ios::out | ios::app);
    if (!out) {
      cerr << "OutLog: error opening <" << filename << "> for append\n";
      return;
    }
  }

				// Generate header line
  if (firstime) {

    firstime = false;

    nbodies = vector<int>(comp.ncomp);
    nbodies1 = vector<int>(comp.ncomp);

    used0 = vector<int>(comp.ncomp);

    mtot = vector<double>(comp.ncomp);
    mtot1 = vector<double>(comp.ncomp);

    com = vector<dvector>(comp.ncomp);
    com1 = vector<dvector>(comp.ncomp);

    cov = vector<dvector>(comp.ncomp);
    cov1 = vector<dvector>(comp.ncomp);

    angm = vector<dvector>(comp.ncomp);
    angm1 = vector<dvector>(comp.ncomp);

    ctr = vector<dvector>(comp.ncomp);

    for (int i=0; i<comp.ncomp; i++) {
      com[i] = vector<double>(3);
      com1[i] = vector<double>(3);
      cov[i] = vector<double>(3);
      cov1[i] = vector<double>(3);
      angm[i] = vector<double>(3);
      angm1[i] = vector<double>(3);
      ctr[i] = vector<double>(3);
    }

    com0 = vector<double>(3);
    cov0 = vector<double>(3);
    angmG = vector<double>(3);
    angm0 = vector<double>(3);
    pos0 = vector<double>(3);
    vel0 = vector<double>(3);

    posL = vector<double>(3);
    velL = vector<double>(3);
    
    comG = vector<double>(3);
    covG = vector<double>(3);
    
    ektot = vector<double>(comp.ncomp);
    ektot1 = vector<double>(comp.ncomp);
    eptot = vector<double>(comp.ncomp);
    eptot1 = vector<double>(comp.ncomp);
    eptotx = vector<double>(comp.ncomp);
    eptotx1 = vector<double>(comp.ncomp);
    clausius = vector<double>(comp.ncomp);
    clausius1 = vector<double>(comp.ncomp);

    if (myid==0) {

      if (restart) {

	out->close();

	// Backup up old file
	string backupfile = filename + ".bak";
	if (rename(filename.c_str(), backupfile.c_str())) {
	  perror("OutLog::Run()");
	  ostringstream message;
	  message << "OutLog: error creating backup file <" 
		  << backupfile << ">";
	  // bomb(message.str());
	}
	
	// Open new output stream for writing
	out = new ofstream(filename.c_str());
	if (!*out) {
	  ostringstream message;
	  message << "OutLog: error opening new log file <" 
		  << filename << "> for writing";
	  bomb(message.str());
	}
	  
	// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  ostringstream message;
	  message << "OutLog: error opening original log file <" 
		  << backupfile << "> for reading";
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
	  *out << cbuffer << "\n";
	  if (line.find_first_of("Time") != string::npos) break;
	}
	
	while (in) {
	  in.getline(cbuffer, cbufsiz);
	  if (!in) break;
	  string line(cbuffer);
	  
	  StringTok<string> toks(line);
	  ttim  = atof(toks(" ").c_str());
	  if (tnow < ttim) break;
	  *out << cbuffer << "\n";
	}
	
	in.close();

      } else {

	string field;
				// Global stanza
	*out << setfill('-') << setw(cwid) << "Global stats";
	for (int i=1; i<num_global; i++) 
	  *out << "|" << setfill(' ') << setw(cwid) << " ";

      
				// Component stanzas
	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  c = *cc;
	  *out << "|" << setw(cwid) << c->id.c_str();
	  for (int i=1; i<num_component; i++) 
	    *out << "|" << setfill(' ') << setw(cwid) << " ";
	}
	*out << endl;
    
				// Global divider
	*out << setfill('-') << setw(cwid) << "-";
	for (int i=1; i<num_global; i++) 
	  *out << "+" << setfill('-') << setw(cwid)  << "-";
      
				// Component dividers
	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  for (int i=0; i<num_component; i++) 
	    *out << "+" << setfill('-') << setw(cwid) << "-";
	}
	*out << endl << setfill(' ');


				// Global labels
	*out << setfill(' ') << setw(cwid) << lab_global[0];
	for (int i=1; i<num_global; i++) *out << "|" << setw(cwid) << lab_global[i];
    
				// Component labels
	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  for (int i=0; i<num_component; i++) {
	    string label = (*cc)->name + " " + lab_component[i];
	    if (label.size()<=cwid)
	      *out << "|" << setw(cwid) << label.c_str();
	    else
	      *out << "|" << label;
	  }
	}
	*out << endl;

				// Global divider
	*out << setfill('-') << setw(cwid) << "-";
	for (int i=1; i<num_global; i++) 
	  *out << "+" << setfill('-') << setw(cwid) << "-";
	
				// Component dividers
	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  for (int i=0; i<num_component; i++) 
	    *out << "+" << setfill('-') << setw(cwid) << "-";
	}
	*out << endl << setfill(' ');
	
				// Global count
	int count=0;
	{
	  ostringstream slab;
	  slab << "[" << ++count << "]";
	  *out << setfill(' ') << setw(cwid) << slab.str();
	}
	for (int i=1; i<num_global; i++) {
	  ostringstream slab;
	  slab << "[" << ++count << "]";
	  *out << "|" << setw(cwid) << slab.str();
	}    
				// Component count
	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  for (int i=0; i<num_component; i++) {
	    ostringstream slab;
	    slab << "[" << ++count << "]";
	    *out << "|" << setw(cwid) << slab.str();
	  }
	}
	*out << endl;

				// Global divider
	*out << setfill('-') << setw(cwid) << "-";
	for (int i=1; i<num_global; i++) 
	  *out << "+" << setfill('-') << setw(cwid) << "-";
	
				// Component dividers
	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  for (int i=0; i<num_component; i++) 
	    *out << "+" << setfill('-') << setw(cwid) << "-";
	}
	*out << endl << setfill(' ');
	
      }
    }
    
  }
  
  if (n % nint && !last) return;


				// Use MPI wall clock to time step
  double wtime = 0.0;

  if (n>laststep) {
    curwtime = MPI_Wtime();
    wtime = (curwtime-lastwtime)/(n-laststep);
    lastwtime = curwtime;
    laststep = n;
  }

				// Zero out accumulators
  for (int i=0; i<comp.ncomp; i++) {

    nbodies[i] = nbodies1[i] = 0;
    used0[i] = 0;
    mtot[i] = mtot1[i] = 0.0;
    ektot[i] = ektot1[i] = 0.0;
    eptot[i] =  eptot1[i] = 0.0;
    eptotx[i] =  eptotx1[i] = 0.0;
    clausius[i] = clausius1[i] = 0.0;

    for (int j=0; j<3; j++) {
      com[i][j] = com1[i][j] = 0.0;
      cov[i][j] = cov1[i][j] = 0.0;
      angm[i][j] = angm1[i][j] = 0.0;
    }

  }

				// Global
  for (int j=0; j<3; j++) {
    comG[j] = com0[j] = 0.0;
    covG[j] = cov0[j] = 0.0;
    angmG[j] = angm0[j] = 0.0;
  }


				// Collect info
  unsigned ntot;
  int indx = 0;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
  
    nbodies1[indx] = c->Number();

    PartMapItr it = c->Particles().begin();
    unsigned long i;

    for (int q=0; q<nbodies1[indx]; q++) {

      i = (it++)->first;

      if (c->freeze(i)) continue;

      mtot1[indx] +=  c->Mass(i);
      
      for (int k=0; k<3; k++) {
	pos0[k] = c->Pos(i, k, Component::Inertial);
	vel0[k] = c->Vel(i, k, Component::Inertial);
	posL[k] = c->Pos(i, k, Component::Local);
	velL[k] = c->Vel(i, k, Component::Local);
      }

      for (int k=0; k<3; k++) {
	com1[indx][k] += c->Mass(i)*posL[k];
	cov1[indx][k] += c->Mass(i)*velL[k];
	comG[k] += c->Mass(i)*pos0[k];
	covG[k] += c->Mass(i)*vel0[k];
      }

      angm1[indx][0] += c->Mass(i)*(posL[1]*velL[2] - posL[2]*velL[1]);
      angm1[indx][1] += c->Mass(i)*(posL[2]*velL[0] - posL[0]*velL[2]);
      angm1[indx][2] += c->Mass(i)*(posL[0]*velL[1] - posL[1]*velL[0]);

      angmG[0] += c->Mass(i)*(pos0[1]*vel0[2] - pos0[2]*vel0[1]);
      angmG[1] += c->Mass(i)*(pos0[2]*vel0[0] - pos0[0]*vel0[2]);
      angmG[2] += c->Mass(i)*(pos0[0]*vel0[1] - pos0[1]*vel0[0]);

      for (int k=0; k<3; k++) pos0[k] = c->Pos(i, k, Component::Centered);
      
      eptot1[indx] += 0.5*c->Mass(i)*c->Part(i)->pot;
      eptotx1[indx] += c->Mass(i)*c->Part(i)->potext;
      for (int k=0; k<3; k++) {
	ektot1[indx] += 0.5*c->Mass(i)*velL[k]*velL[k];
	clausius1[indx] += c->Mass(i)*posL[k]*c->Part(i)->acc[k];
      }
    }

    for (int k=0; k<3; k++) ctr[indx][k] = c->center[k];

    used0[indx] = c->force->Used();

    indx++;
  }
				// Send back to Process 0

  MPI_Reduce(&nbodies1[0], &nbodies[0], comp.ncomp, MPI_INT, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&mtot1[0], &mtot[0], comp.ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  for (int i=0; i<comp.ncomp; i++) {
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

  MPI_Reduce(&ektot1[0], &ektot[0], comp.ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&eptot1[0], &eptot[0], comp.ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&eptotx1[0], &eptotx[0], comp.ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&clausius1[0], &clausius[0], comp.ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);


  if (myid == 0) {

    // =============
    // Global
    // =============

				// Current time
    *out << setw(cwid) << tnow;

    double mtot0 = 0.0;
    for (int i=0; i<comp.ncomp; i++) mtot0 += mtot[i];

				// Total mass
    *out << "|" << setw(cwid) << mtot0;

				// Total number
    int nbodies0 = 0;
    for (int i=0; i<comp.ncomp; i++) nbodies0 += nbodies[i];
    *out << "|" << setw(cwid) << nbodies0;

				// COM
    for (int j=0; j<3; j++)
      if (mtot0>0.0)
	*out << "|" << setw(cwid) << com0[j]/mtot0;
      else
	*out << "|" << setw(cwid) << 0.0;


				// COV
    for (int j=0; j<3; j++)
      if (mtot0>0.0)
	*out << "|" << setw(cwid) << cov0[j];
      else
	*out << "|" << setw(cwid) << 0.0;
	
				// Ang mom
    for (int j=0; j<3; j++)
      *out << "|" << setw(cwid) << angm0[j];
    
				// KE
    double ektot0 = 0.0;
    for (int i=0; i<comp.ncomp; i++) ektot0 += ektot[i];
    *out << "|" << setw(cwid) << ektot0;
      
				// PE
    double eptot0 = 0.0;
    for (int i=0; i<comp.ncomp; i++) eptot0 += eptot[i] + 0.5*eptotx[i];
    *out << "|" << setw(cwid) << eptot0;
     
				// Clausius, Total, 2T/VC
    double clausius0 = 0.0;
    for (int i=0; i<comp.ncomp; i++) clausius0 += clausius[i];
    *out << "|" << setw(cwid) << clausius0;
    *out << "|" << setw(cwid) << ektot0 + clausius0;
    if (clausius0 != 0.0)
      *out << "|" << setw(cwid) << -2.0*ektot0/clausius0;
    else
      *out << "|" << setw(cwid) << 0.0;

    *out << "|" << setw(cwid) << wtime;
    int usedT = 0;
    for (int i=0; i<comp.ncomp; i++) usedT += used0[i];
    *out << "|" << setw(cwid) << usedT;


    // =============
    // Per component
    // =============


    for (int i=0; i<comp.ncomp; i++) {

      *out << "|" << setw(cwid) << mtot[i];
      *out << "|" << setw(cwid) << nbodies[i];
      for (int j=0; j<3; j++)
	if (mtot[i]>0.0)
	  *out << "|" << setw(cwid) << com[i][j]/mtot[i];
	else
	  *out << "|" << setw(cwid) << 0.0;
      for (int j=0; j<3; j++)
	if (mtot[i]>0.0)
	  *out << "|" << setw(cwid) << cov[i][j]/mtot[i];
	else
	  *out << "|" << setw(cwid) << 0.0;
      for (int j=0; j<3; j++)
	*out << "|" << setw(cwid) << angm[i][j];
      for (int j=0; j<3; j++)
	*out << "|" << setw(cwid) << ctr[i][j];

      double vbar2=0.0;		// Kinetic energy in per component
      if (mtot[i]>0.0) {	// center of velocity frame
	for (int j=0; j<3; j++)
	  vbar2 +=  cov[i][j]*cov[i][j];
	vbar2 /=  mtot[i]*mtot[i];
      }
      ektot[i] -= 0.5*mtot[i]*vbar2; // Update KE
      
      *out << "|" << setw(cwid) << ektot[i];
      *out << "|" << setw(cwid) << eptot[i] + eptotx[i];
      *out << "|" << setw(cwid) << clausius[i];
      *out << "|" << setw(cwid) << ektot[i] + clausius[i];
      if (clausius[i] != 0.0)
	*out << "|" << setw(cwid) << -2.0*ektot[i]/clausius[i];
      else
	*out << "|" << setw(cwid) << 0.0;
      *out << "|" << setw(cwid) << used0[i];
    }

    *out << endl;

    out->close();
    delete out;
  }

}

