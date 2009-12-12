#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.h"

#include <UserPeriodic.H>


UserPeriodic::UserPeriodic(string &line) : ExternalForce(line)
{

  id = "PeriodicBC";		// Periodic boundary condition ID

				// Sizes in each dimension
  L = vector<double>(3, 2.0);	
				// Center offset in each dimension
  offset = vector<double>(3, 1.0);

  bc = "ppp";			// Periodic BC in all dimensions

  comp_name = "";		// Default component (must be specified)

  nbin = 0;			// Number of bins in trace (0 means no trace)
  dT = 1.0;			// Interval for trace
  tcol = -1;			// Column for temperture info (ignored if <0)
  trace = false;		// Tracing off until signalled

  initialize();

				// Look for the fiducial component
  bool found = false;
  list<Component*>::iterator cc;
  Component *c;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    if ( !comp_name.compare(c->name) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  if (!found) {
    cerr << "UserPeriodic: process " << myid 
	 << " can't find fiducial component <" << comp_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }
  
  // 
  // Initialize structures for trace data
  //
  if (nbin) {
    mbinT = vector< vector<double> >(nthrds);
    tbinT = vector< vector<double> >(nthrds);
    for (int n=0; n<nthrds; n++) {
      mbinT[n] = vector<double>(nbin, 0);
      tbinT[n] = vector<double>(nbin, 0);
    }
    mbin = vector<double>(nbin);
    tbin = vector<double>(nbin);
    Tnext = tnow;
    dX = L[0]/nbin;
  }

  userinfo();
}

UserPeriodic::~UserPeriodic()
{
}

void UserPeriodic::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine PERIODIC BOUNDARY CONDITION initialized"
       << " using component <" << comp_name << ">";
  if (nbin) cout << " with gas trace, dT=" << dT 
		 << ", nbin=" << nbin << ", tcol=" << tcol;
  cout << endl;

  cout << "   Cube sides (x , y , z) = (" 
       << L[0] << " , " 
       << L[1] << " , " 
       << L[2] << " ) " << endl; 

  cout << "Center offset (x , y , z) = (" 
       << offset[0] << " , " 
       << offset[1] << " , " 
       << offset[2] << " ) " << endl; 

  cout << "Boundary type (x , y , z) = (" 
       << bc[0] << " , " 
       << bc[1] << " , " 
       << bc[2] << " ) " << endl;

  print_divider();
}

void UserPeriodic::initialize()
{
  string val;

  if (get_value("compname", val))	comp_name = val;

  if (get_value("sx", val))	        L[0] = atof(val.c_str());
  if (get_value("sy", val))	        L[1] = atof(val.c_str());
  if (get_value("sz", val))	        L[2] = atof(val.c_str());

  if (get_value("cx", val))	        offset[0] = atof(val.c_str());
  if (get_value("cy", val))	        offset[1] = atof(val.c_str());
  if (get_value("cz", val))	        offset[2] = atof(val.c_str());

  if (get_value("dT", val))	        dT = atof(val.c_str());
  if (get_value("nbin", val))	        nbin = atoi(val.c_str());
  if (get_value("tcol", val))	        tcol = atoi(val.c_str());

  if (get_value("btype", val)) {
    if (strlen(val.c_str()) >= 3) {
      for (int k=0; k<3; k++) {
	switch (val.c_str()[k]) {
	case 'p':
	  bc[k] = 'p';		// Periodic
	  break;
	case 'r':
	  bc[k] = 'r';		// Reflection
	  break;
	default:
	  bc[k] = 'v';		// Vacuum
	  break;
	}
      }
    }
  }
  
}


void UserPeriodic::determine_acceleration_and_potential(void)
{
  if (nbin && tnow>=Tnext) trace = true;
  exp_thread_fork(false);
  if (trace) write_trace();
  print_timings("UserPeriodic: thread timings");
}

void UserPeriodic::write_trace()
{
  for (int n=1; n<nthrds; n++) {
    for (int k=0; k<nbin; k++) {
      mbinT[0][k] += mbinT[n][k];
      tbinT[0][k] += tbinT[n][k];
    }
  }

  MPI_Reduce(&mbinT[0][0], &mbin[0], nbin, MPI_DOUBLE, MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  MPI_Reduce(&tbinT[0][0], &tbin[0], nbin, MPI_DOUBLE, MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  if (myid==0) {
    ostringstream fout;
    fout << outdir << runtag << ".shocktube_trace";
    ofstream out(fout.str().c_str(), ios::app);
    for (int k=0; k<nbin; k++)
      out << setw(18) << tnow
	  << setw(18) << dX*(0.5+k) - offset[0]
	  << setw(18) << mbin[k]/(dX*L[1]*L[2])
	  << setw(18) << tbin[k]/(mbin[k]+1.0e-10)
	  << endl;
    out << endl;
  }

  //
  // Clean data structures for next call
  //
  for (int n=0; n<nthrds; n++) {
    for (int k=0; k<nbin; k++)
      mbinT[n][k] = tbinT[n][k] = 0.0;
  }

  trace = false;		// Tracing off until
  Tnext += dT;			// tnow>=Tnext
}


void * UserPeriodic::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  
  thread_timing_beg(id);

  double pos, delta;
  PartMapItr it = cC->Particles().begin();
  
  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    
				// Index for the current particle
    unsigned long i = (it++)->first;
    
    Particle *p = cC->Part(i);

				// If we are multistepping, compute BC
				// only at or above this level
    // if (!multistep || (cC->Part(i)->level >= mlevel)) {

      for (int k=0; k<3; k++) {
				// Increment so that the positions range
				// between 0 and L[k]
	pos = p->pos[k] + offset[k];

	//
	// Reflection BC
	//
	if (bc[k] == 'r') {
	  if (pos < 0.0) {
	    delta = -pos - L[k]*floor(-pos/L[k]);
	    p->pos[k] = delta;
	    p->vel[k] = -p->vel[k];
	  } 
	  if (pos >= L[k]) {
	    delta = pos - L[k]*floor(pos/L[k]);
	    p->pos[k] =  L[k] - delta;
	    p->vel[k] = -p->vel[k];
	  }
	}
	
	//
	// Periodic BC
	//
	if (bc[k] == 'p') {
	  if (pos < 0.0) {
	    p->pos[k] = p->pos[k] + L[k]*floor(1.0+fabs(pos/L[k]));
	  }
	  if (pos >= L[k]) {
	    p->pos[k] = p->pos[k] - L[k]*floor(fabs(pos/L[k]));
	  }
	}
	
				// Replace the offset
	pos = p->pos[k] - offset[k];
      }

      //
      // Sanity check for this particle
      //
      for (int k=0; k<3; k++) {
	if (bc[k] != 'v') {
	  if (p->pos[k] < -offset[k] || p->pos[k] >= L[k]-offset[k]) {
	    cout << "Process " << myid << " id=" << id 
		 << ": Error in pos[" << k << "]=" << p->pos[k] << endl;
	  }
	}
      }

      // }

    //
    // Acccumlate data for shocktube trace
    //
    if (trace) {
      int indx = static_cast<int>(floor((p->pos[0]+offset[0])/dX));
      if (indx>=0 && indx<nbin) {
	mbinT[id][indx] += p->mass;
	if (tcol>=0 && tcol<static_cast<int>(p->dattrib.size()))
	  tbinT[id][indx] += p->mass*p->dattrib[tcol];
      }
    }

  }
  
  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerPeriodic(string& line)
  {
    return new UserPeriodic(line);
  }
}

class proxypbc { 
public:
  proxypbc()
  {
    factory["userperiodic"] = makerPeriodic;
  }
};

proxypbc p;
