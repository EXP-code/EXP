/* 
  Multi-timestepping support routines
*/

#include "expand.h"
#include <sstream>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif


void sync_eval_multistep(void)
{
  comp.multistep_reset();
				// This forces interpolation to evaluate to 
				// the last set of computed coefficients
  mstep = Mstep;
  for (int M=0; M<=multistep; M++) {
    stepL[M] = 0;
    stepN[M] = Mstep;
  }
}


/// Helper class to pass info to threaded multistep update routine
struct thrd_pass_sync 
{
  //! Level
  int level;

  //! Thread counter id
  int id;

  //! Component
  Component *c;
};

/// Count offgrid particles in the threads
vector< vector<unsigned> > off1;
vector< double > mindt1;
vector< double > maxdt1;

/// Type counter
vector< vector< vector<unsigned> > > tmdt;

//
// The threaded routine
//
void * adjust_multistep_level_thread(void *ptr)
{
  // Level
  int level = static_cast<thrd_pass_sync*>(ptr)->level;

  // Thread ID
  int   id = static_cast<thrd_pass_sync*>(ptr)->id;

  // Component
  Component *c = static_cast<thrd_pass_sync*>(ptr)->c;

  // Examine all time steps at or below this level and compute timestep
  // criterion and adjust level if necessary

  double dt, dtv, dta, dtr, dsr;
  int npart = c->levlist[level].size();
  int offgrid = 0, lev;

  //
  // Compute the beginning and end points for threads
  //
  int nbeg = npart*id    /nthrds;
  int nend = npart*(id+1)/nthrds;

  //
  // The particle loop
  //
  for (int i=nbeg; i<nend; i++) {

    double rtot = 0.0;
    double vtot = 0.0;
    double atot = 0.0;
      
    int n = c->levlist[level][i];
    
    for (int k=0; k<c->dim; k++) {
      rtot += c->Part(n)->pos[k]*c->Part(n)->pos[k];
      vtot += c->Part(n)->vel[k]*c->Part(n)->vel[k];
      atot += c->Part(n)->acc[k]*c->Part(n)->acc[k];
    }
    rtot = sqrt(rtot);
    vtot = sqrt(vtot) + 1.0e-18;
    atot = sqrt(atot) + 1.0e-18;

    dsr = c->Part(n)->scale;
    if (dsr>0) rtot = min<double>(rtot, dsr);

    dtv = dynfracV*rtot/vtot;
    dta = dynfracA*sqrt(rtot/atot);
    dt  = min<double>(dtv, dta);

    dtr = c->Part(n)->dtreq;
    if (dtr>0) dt = min<double>(dt, dtr);

    mindt1[id] = min<double>(mindt1[id], dt);
    maxdt1[id] = max<double>(maxdt1[id], dt);
	
    // Tally smallest (e.g. controlling) timestep
    if (dtv<dta) {
      if (dtr<=0 || dtv<dtr)
	tmdt[id][level][0]++;  
      else
	tmdt[id][level][2]++;
    } else {
      if (dtr<=0 || dta<dtr)
	tmdt[id][level][1]++;  
      else
	tmdt[id][level][2]++;
    }
    // Counter
    tmdt[id][level][3]++;

    if (dt>dtime) lev = 0;
    else lev = (int)floor(log(dtime/dt)/log(2.0));
    
    if (lev>multistep) {
      lev = multistep;
      offgrid++;
    }
    
    unsigned plev = c->Part(n)->level;
    unsigned nlev = lev;

    // Sanity check
    if (level != plev) {
      cerr << "Process " << myid << " id=" << id 
	   << ": n=" << n << " level=" << level 
	   << " tlevel=" << plev
	   << " plevel=" << c->Part(n)->level << endl;
    }

    if (nlev != level) {
      //
      // Update coefficients
      //
      c->force->multistep_update(plev, nlev, c, n, id);
      c->Part(n)->level = lev;
    }
  }

  if (VERBOSE>0 && (this_step % 100 == 0)) {
    off1[id].push_back(offgrid);
  }
  
  return (NULL);
}


void adjust_multistep_level(bool all)
{
  if (!multistep) return;

  //
  // Begin the update
  //
  for (list<Component*>::iterator cc=comp.components.begin(); 
       cc != comp.components.end(); cc++)
    (*cc)->force->multistep_update_begin();

  //
  // Preliminary data structure and thread creation
  //
  mindt1 = vector< double > (nthrds,  1.0e20);
  maxdt1 = vector< double > (nthrds, -1.0e20);
  off1 = vector< vector<unsigned> > (nthrds);

  thrd_pass_sync* td = new thrd_pass_sync [nthrds];

  if (!td) {
    cerr << "Process " << myid
	 << ": adjust_multistep_level: error allocating thread structures\n";
    exit(18);
  }

  pthread_t* t  = new pthread_t [nthrds];

  if (!t) {
    cerr << "Process " << myid
	 << ": adjust_multistep_level: error allocating memory for thread\n";
    exit(18);
  }
  

  tmdt = vector< vector< vector<unsigned> > >(nthrds);
  for (int n=0; n<nthrds; n++) {
    tmdt[n] = vector< vector<unsigned> >(multistep+1);
    for (int k=0; k<=multistep; k++) tmdt[n][k] = vector<unsigned>(4);
  }

  for (list<Component*>::iterator cc=comp.components.begin();
       cc != comp.components.end(); cc++) {
    
    for (int n=0; n<nthrds; n++)
      for (int k=0; k<=multistep; k++) 
	for (int j=0; j<4; j++) tmdt[n][k][j] = 0;

    for (int level=0; level<=multistep; level++) {
      
      if (all || mactive[mstep-1][level]) {

	if (nthrds==1) {

	  td[0].level = level;
	  td[0].id = 0;
	  td[0].c = *cc;

	  adjust_multistep_level_thread(&td[0]);

	} else {
	  

	  //
	  // Make the <nthrds> threads
	  //
	  int errcode;
	  void *retval;
  
	  for (int i=0; i<nthrds; i++) {
	    
	    td[i].level = level;
	    td[i].id = i;
	    td[i].c = *cc;
	    
	    errcode =  pthread_create(&t[i], 0, adjust_multistep_level_thread, &td[i]);
	    
	    if (errcode) {
	      cerr << "Process " << myid
		   << " adjust_multistep_level: cannot make thread " << i
		   << ", errcode=" << errcode << endl;
	      exit(19);
	    }
#ifdef DEBUG
	    else {
	      cout << "Process " << myid << ": thread <" << i << "> created\n";
	    }
#endif
	  }
	  
	  //
	  // Collapse the threads
	  //
	  for (int i=0; i<nthrds; i++) {
	    if ((errcode=pthread_join(t[i], &retval))) {
	      cerr << "Process " << myid
		   << " adjust_multistep_level: thread join " << i
		   << " failed, errcode=" << errcode << endl;
	      exit(20);
	    }
#ifdef DEBUG    
	    cout << "Process " << myid 
		 << ": multistep thread <" << i << "> thread exited\n";
#endif
	  }
  
	}
	
      }
    }

    for (int n=0; n<nthrds; n++)
      for (int k=0; k<=multistep; k++) 
	for (int j=0; j<4; j++) 
	  (*cc)->mdt_ctr[k][j] += tmdt[n][k][j];
  }

  delete [] td;
  delete [] t;

  //
  // Finish the update
  //
  for (list<Component*>::iterator cc=comp.components.begin(); 
       cc != comp.components.end(); cc++) {
    (*cc)->reset_level_lists();
    (*cc)->force->multistep_update_finish();
  }


  //
  // Diagnostic output
  //
  if (VERBOSE>0) {

    for (int n=1; n<nthrds; n++) {
      for (unsigned j=0; j<off1[0].size(); j++) off1[0][j] += off1[n][j];
      mindt1[0] = min<double>(mindt1[0], mindt1[n]);
      maxdt1[0] = max<double>(maxdt1[0], maxdt1[n]);
    }

    vector<unsigned> off(off1[0].size(), 0);
    double mindt, maxdt;

    MPI_Reduce(&mindt1[0], &mindt, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&maxdt1[0], &maxdt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&off1[0][0], &off[0], off1.size(), MPI_UNSIGNED, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
    
    if (myid==0) {
      
      unsigned sum=0;
      for (unsigned i=0; i<off.size(); i++) sum += off[i];
      
      if (sum) {
	cout << setw(70) << setfill('-') << '-' << endl
	     << left << "--- Multistepping overrun" << endl
	     << left << "--- Min DT=" << mindt << "  Max DT=" << maxdt << endl
	     << setw(70) << setfill('-') << '-' << endl 
	     << setfill(' ') << right;

	vector<unsigned>::iterator it = off.begin();

	list<Component*>::iterator cc;
	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  //	  if (*it > 0) {
	    ostringstream sout;
	    sout << "Component <" << (*cc)->name << ">";
	    cout << setw(20) << sout.str() << ": "
		 << *it << "/" << (*cc)->nbodies_tot << endl;
	    // }
	  it++;
	}
	cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
      }
    }

  }
}


void initialize_multistep()
{
  Mstep = 1 << multistep;	// Number of substeps

				// Number of substeps at each level
  mintvl.push_back(Mstep);
  for (int n=1; n<=multistep; n++) mintvl.push_back(mintvl.back()/2);

				// Set up the step-level bool array
  for (int ms=0; ms<Mstep; ms++)
    mactive.push_back(vector<bool>(multistep+1, false));

				// Find and save the active levels at each step
  for (int ms=1; ms<=Mstep; ms++) {
    for (unsigned mlevel=0; mlevel<=multistep; mlevel++) {
      if ( (ms % (1 << (multistep-mlevel))) == 0) 
	mactive[ms-1][mlevel] = true;
    }
  }

  stepL  = vector<int>(multistep+1);
  stepN  = vector<int>(multistep+1);

  mfirst = vector<int>(Mstep);
				// Lowest active level at each step
  for (int ms=1; ms<=Mstep; ms++) {
    for (int mlevel=0; mlevel<=multistep; mlevel++) {
      if (mactive[ms-1][mlevel]) {
	mfirst[ms-1] = mlevel;
	break;
      }
    }
  }

  if (VERBOSE>10 && myid==0 && multistep) {
    cout << setw(70) << setfill('-') << '-' << endl
	 << setw(70) << left << "--- Multistep level control structure " << endl
	 << setw(70) << setfill('-') << '-' << endl << setfill(' ') << right;
    cout << "Mstep=" << Mstep << "  levels=" << multistep << endl << endl;
    cout << setw(4) << "Step" << "/" << setw(3)<< "1st" << "  ";
    for (int mlevel=0; mlevel<=multistep; mlevel++) cout << setw(3) << mlevel;
    cout << endl;
    for (int ms=1; ms<=Mstep; ms++) {
      cout << setw(4) << ms << "/";
      cout << setw(3) << mfirst[ms-1] << ": ";
      for (int mlevel=0; mlevel<=multistep; mlevel++)
	if (mactive[ms-1][mlevel]) cout << setw(3) << "+";
	else                       cout << setw(3) << "-";
      cout << endl;
    }
    cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
  }

}
