#include "expand.h"

static char rcsid[] = "$Id$";

#include <OutLog.H>

char OutLog::lab_global[][19] = {
  "Current time",
  "Total mass", 
  "Number of bodies", 
  "COM (x)", 
  "COM (y)", 
  "COM (z)", 
  "COV (x)", 
  "COV (y)", 
  "COV (z)", 
  "Ang mom (x)", 
  "Ang mom (y)", 
  "Ang mom (z)", 
  "Total KE", 
  "Total PE", 
  "Virial of C", 
  "Total E", 
  "2T/VC", 
  "Time per step", 
  "Total # used"
};

char OutLog::lab_component[][20] = {
  "Total mass", 
  "Number of bodies", 
  "COM (x)", 
  "COM (y)", 
  "COM (z)", 
  "COV (x)", 
  "COV (y)", 
  "COV (z)", 
  "Ang mom (x)", 
  "Ang mom (y)", 
  "Ang mom (z)", 
  "CTR (x)", 
  "CTR (y)", 
  "CTR (z)", 
  "Total KE", 
  "Total PE", 
  "Virial of C", 
  "Total E", 
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
    filename = "OUTLOG.dat";
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
  ofstream *out;

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

    mom = vector<dvector>(comp.ncomp);
    mom1 = vector<dvector>(comp.ncomp);

    for (int i=0; i<comp.ncomp; i++) {
      com[i] = vector<double>(3);
      com1[i] = vector<double>(3);
      cov[i] = vector<double>(3);
      cov1[i] = vector<double>(3);
      angm[i] = vector<double>(3);
      angm1[i] = vector<double>(3);
      ctr[i] = vector<double>(3);
      mom[i] = vector<double>(6);
      mom1[i] = vector<double>(6);
    }

    com0 = vector<double>(3);
    cov0 = vector<double>(3);
    angm0 = vector<double>(3);
    mom0 = vector<double>(5);
    
    ektot = vector<double>(comp.ncomp);
    ektot1 = vector<double>(comp.ncomp);
    eptot = vector<double>(comp.ncomp);
    eptot1 = vector<double>(comp.ncomp);
    clausius = vector<double>(comp.ncomp);
    clausius1 = vector<double>(comp.ncomp);

    if (myid==0) {
      string field;
				// Global stanza
      *out << "--" << setw(16) << "Global stats";
      for (int i=1; i<num_global; i++) *out << "|" << setw(18) << " ";
      
				// Component stanzas
      for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	c = *cc;
	*out << "|" << setw(18) << c->id.c_str();
	for (int i=1; i<num_component; i++) *out << "|" << setw(18) << " ";
      }
      *out << endl;
    
				// Global divider
      *out << setw(18) << "------------------";
      for (int i=1; i<num_global; i++) 
	*out << "+" << setw(18) << "------------------";
      
				// Component dividers
      for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	for (int i=0; i<num_component; i++) 
	  *out << "+" << setw(18) << "------------------";
      }
      *out << endl;

				// Global labels
      *out << setw(18) << lab_global[0];
      for (int i=1; i<num_global; i++) *out << "|" << setw(18) << lab_global[i];
    
				// Component labels
      for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	for (int i=0; i<num_component; i++) 
	  *out << "|" << setw(18) << lab_component[i];
      }
      *out << endl;

				// Global divider
      *out << setw(18) << "------------------";
      for (int i=1; i<num_global; i++) 
	*out << "+" << setw(18) << "------------------";
    
				// Component dividers
      for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	for (int i=0; i<num_component; i++) 
	  *out << "+" << setw(18) << "------------------";
      }
      *out << endl;
  
    }
    
  }
  
  if (n % nint && !last) return;

  int i;
				/* Use MPI wall clock to time step */
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
    clausius[i] = clausius1[i] = 0.0;

    for (int j=0; j<3; j++) {
      com[i][j] = com1[i][j] = 0.0;
      cov[i][j] = cov1[i][j] = 0.0;
      angm[i][j] = angm1[i][j] = 0.0;
    }

    for (int j=0; j<6; j++)
      mom[i][j] = mom1[i][j] = 0.0;

  }

				// Collect info
  vector<Particle>::iterator p;
  int indx = 0;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
  
    nbodies1[indx] = c->particles.size();

    for (p=c->particles.begin(); p!=c->particles.end(); p++) {

      if (p->freeze()) continue;

      mtot1[indx] +=  p->mass;
      
      for (int j=0; j<3; j++) {
	com1[indx][j] += p->mass*p->pos[j];
	cov1[indx][j] += p->mass*p->vel[j];
	mom1[indx][j] += p->mass*p->pos[j]*p->pos[j];
      }

      // Moment matrix off-diagonal terms
      mom1[indx][3] += p->mass*p->pos[0]*p->pos[1];
      mom1[indx][4] += p->mass*p->pos[0]*p->pos[2];
      mom1[indx][5] += p->mass*p->pos[1]*p->pos[2];

      angm1[indx][0] += p->mass*(p->pos[1]*p->vel[2]-p->pos[2]*p->vel[1]);
      angm1[indx][1] += p->mass*(p->pos[2]*p->vel[0]-p->pos[0]*p->vel[2]);
      angm1[indx][2] += p->mass*(p->pos[0]*p->vel[1]-p->pos[1]*p->vel[0]);

      eptot1[indx] += 0.5*p->mass*p->pot + p->mass*p->potext;
      ektot1[indx] += 0.5*p->mass*
	(p->vel[0]*p->vel[0]+p->vel[1]*p->vel[1]+p->vel[2]*p->vel[2]);
      clausius1[indx] += p->mass*
	(p->pos[0]*p->acc[0]+p->pos[1]*p->acc[1]+p->pos[2]*p->acc[2]);
    }
    
    for (int j=0; j<3; j++) ctr[indx][j] = c->center[j];

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
    MPI_Reduce(&mom1[i][0], &mom[i][0], 6, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    MPI_Reduce(&angm1[i][0], &angm[i][0], 3, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
  }
  MPI_Reduce(&ektot1[0], &ektot[0], comp.ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&eptot1[0], &eptot[0], comp.ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&clausius1[0], &clausius[0], comp.ncomp, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);


  if (myid == 0) {

    // =============
    // Global
    // =============

				// Current time
    *out << setw(18) << tnow;

    double mtot0 = 0.0;
    for (int i=0; i<comp.ncomp; i++) mtot0 += mtot[i];

				// Total mass
    *out << "|" << setw(18) << mtot0;

				// Total number
    int nbodies0 = 0;
    for (int i=0; i<comp.ncomp; i++) nbodies0 += nbodies[i];
    *out << "|" << setw(18) << nbodies0;

				// Vectors
    vector<double> com0(3), cov0(3), angm0(3), mom0(5);
    for (int j=0; j<3; j++) com0[j] = cov0[j] = angm0[j] = 0.0;
    for (int j=0; j<5; j++) mom0[j] = 0.0;
    for (int i=0; i<comp.ncomp; i++) {
      for (int j=0; j<3; j++) {
	com0[j] += com[i][j];
	cov0[j] += cov[i][j];
	angm0[j] += angm[i][j];
      }
      for (int j=0; j<3; j++)
	mom0[j] += mom[i][j];
    }

				// COM
    for (int j=0; j<3; j++)
      *out << "|" << setw(18) << com0[j];


				// COV
    for (int j=0; j<3; j++)
      *out << "|" << setw(18) << cov0[j];
    
				// Ang mom
    for (int j=0; j<3; j++)
      *out << "|" << setw(18) << angm0[j];
    
				// KE
    double ektot0 = 0.0;
    for (int i=0; i<comp.ncomp; i++) ektot0 += ektot[i];
    *out << "|" << setw(18) << ektot0;
      
				// PE
    double eptot0 = 0.0;
    for (int i=0; i<comp.ncomp; i++) eptot0 += eptot[i];
    *out << "|" << setw(18) << eptot0;
     
				// Clausius, Total, 2T/VC
    double clausius0 = 0.0;
    for (int i=0; i<comp.ncomp; i++) clausius0 += clausius[i];
    *out << "|" << setw(18) << clausius0;
    *out << "|" << setw(18) << ektot0 + clausius0;
    *out << "|" << setw(18) << -2.0*ektot0/clausius0;

    *out << "|" << setw(18) << wtime;
    int usedT = 0;
    for (int i=0; i<comp.ncomp; i++) usedT += used0[i];
    *out << "|" << setw(18) << usedT;


    // =============
    // Per component
    // =============


    for (int i=0; i<comp.ncomp; i++) {

      *out << "|" << setw(18) << mtot[i];
      *out << "|" << setw(18) << nbodies[i];
      for (int j=0; j<3; j++)
	*out << "|" << setw(18) << com[i][j];
      for (int j=0; j<3; j++)
	*out << "|" << setw(18) << cov[i][j];
      for (int j=0; j<3; j++)
	*out << "|" << setw(18) << angm[i][j];
      for (int j=0; j<3; j++)
	*out << "|" << setw(18) << ctr[i][j];
      *out << "|" << setw(18) << ektot[i];
      *out << "|" << setw(18) << eptot[i];
      *out << "|" << setw(18) << clausius[i];
      *out << "|" << setw(18) << ektot[i] + clausius[i];
      *out << "|" << setw(18) << -2.0*ektot[i]/clausius[i];
      *out << "|" << setw(18) << used0[i];
    }

    *out << endl;

    out->close();
  }

}
