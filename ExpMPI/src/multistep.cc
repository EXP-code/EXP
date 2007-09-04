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
struct thrd_pass_sync {
  //! All levels flag
  bool all;
  //! Thread counter id
  int id;
};

// Count offgrid particles in the threads
vector< vector<unsigned> > off1;
vector< vector<int> > dlev;

//
// The threaded routine
//
void * adjust_multistep_level_thread(void *ptr)
{
  // Compute all levels?
  bool all = static_cast<thrd_pass_sync*>(ptr)->all;

  // Thread ID
  int   id = static_cast<thrd_pass_sync*>(ptr)->id;

  // Examine all time steps at or below this level and compute timestep
  // criterion and adjust level if necessary

  list<Component*>::iterator cc;
  Component *c;
  double dt;
  int npart, lev, offgrid;

  if (all) mstep = Mstep;

  //
  // Run through particles in each component
  //

  int nbeg, nend;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    
    c = *cc;

    npart = c->Number();
    offgrid = 0;

    //
    // Compute the beginning and end points for threads
    //
    nbeg = npart*id    /nthrds;
    nend = npart*(id+1)/nthrds;

    //
    // The particle loop
    //
    for (int n=nbeg; n<nend; n++) {

      if (all || mactive[mstep-1][c->Part(n)->level]) {

	double rtot = 0.0;
	double vtot = 0.0;
	double atot = 0.0;

	for (int i=0; i<c->dim; i++) {
	  rtot += c->Part(n)->pos[i]*c->Part(n)->pos[i];
	  vtot += c->Part(n)->vel[i]*c->Part(n)->vel[i];
	  atot += c->Part(n)->acc[i]*c->Part(n)->acc[i];
	}
	rtot = sqrt(rtot);
	vtot = sqrt(vtot) + 1.0e-18;
	atot = sqrt(atot) + 1.0e-18;

	dt = min<double>(dynfracV*rtot/vtot, dynfracA*sqrt(rtot/atot));
	
	if (dt>dtime) lev = 0;
	else lev = (int)floor(log(dtime/dt)/log(2.0));

	if (lev>multistep) {
	  lev = multistep;
	  offgrid++;
	}

	if ( lev != c->Part(n)->level && (all || mactive[mstep-1][lev]) ) {
	  //
	  // Adjust level counts
	  //
	  dlev[id][c->Part(n)->level]--;
	  dlev[id][lev]++;

	  //
	  // Update coefficients
	  //
	  c->force->multistep_update(c->Part(n)->level, lev, c, n, id);
	  c->Part(n)->level = lev;
	}
      }
    }

    if (VERBOSE>0 && (this_step % 100 == 0)) {
      off1[id].push_back(offgrid);
    }

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
  off1 = vector< vector<unsigned> > (nthrds);
  dlev = vector< vector<int> > (nthrds);

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

  //
  // Make the <nthrds> threads
  //
  int errcode;
  void *retval;
  
  for (int i=0; i<nthrds; i++) {

    td[i].all = all;
    td[i].id = i;
    dlev[i] = vector<int>(multistep, 0);

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
    for (int m=0; m<=multistep; m++) levpop[m] += dlev[i][m];
#ifdef DEBUG    
    cout << "Process " << myid << ": multistep thread <" << i << "> thread exited\n";
#endif
  }
  
  delete [] td;
  delete [] t;

  //
  // Finish the update
  //
  for (list<Component*>::iterator cc=comp.components.begin(); 
       cc != comp.components.end(); cc++)
    (*cc)->force->multistep_update_finish();


  //
  // Diagnostic output
  //
  if (VERBOSE>0) {

    for (int n=1; n<nthrds; n++) {
      for (unsigned j=0; j<off1[0].size(); j++) off1[0][j] += off1[n][j];
    }

    vector<unsigned> off(off1[0].size(), 0);

    MPI_Reduce(&off1[0][0], &off[0], off1.size(), MPI_UNSIGNED, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
    
    if (myid==0) {

      unsigned sum=0;
      for (unsigned i=0; i<off.size(); i++) sum += off[i];
      
      if (sum) {
	cout << setw(70) << setfill('-') << '-' << endl
	     << setw(70) << left << "--- Multistepping overrun" << endl
	     << setw(70) << setfill('-') << '-' << endl << setfill(' ') << right;

	vector<unsigned>::iterator it = off.begin();

	list<Component*>::iterator cc;
	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  if (*it > 0) {
	    ostringstream sout;
	    sout << "Component <" << (*cc)->name << ">";
	    cout << setw(20) << sout.str() << ": "
		 << *it << "/" << (*cc)->nbodies_tot << endl;
	  }
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
  levpop = vector<int>(multistep+1, 0);

  for (list<Component*>::iterator cc=comp.components.begin(); 
       cc != comp.components.end(); cc++) levpop[0] = (*cc)->Number();

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

  if (VERBOSE>0 && myid==0 && multistep) {
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
	else                          cout << setw(3) << "-";
      cout << endl;
    }
    cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
  }

}
