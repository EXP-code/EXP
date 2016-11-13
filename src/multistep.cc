/* 
  Multi-timestepping support routines
*/

#include "expand.h"
#include <sstream>
#include <map>


void sync_eval_multistep(void)
{
  comp.multistep_reset();
				// This forces interpolation to evaluate to 
				// the last set of computed coefficients
  mstep = Mstep;
}


//
// Helper class to pass info to threaded multistep update routine
//
struct thrd_pass_sync 
{
  //! Level
  int level;

  //! Thread counter id
  int id;

  //! Component
  Component *c;
};

//
// Count offgrid particles in the threads
//
static map< Component*, vector<unsigned> > offlo1, offhi1;
static vector< double > mindt1;
static vector< double > maxdt1;

// Type counter
static vector< vector< vector<unsigned> > > tmdt;

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

  double dt, dts, dtv, dta, dtA, dtd, dtr, dsr, rtot, vtot, atot, ptot;
  int npart = c->levlist[level].size();
  int offlo = 0, offhi = 0, lev;

  //
  // Compute the beginning and end points for threads
  //
  int nbeg = npart*id    /nthrds;
  int nend = npart*(id+1)/nthrds;

  //
  // Small positive constant
  //
  const double eps = 1.0e-20;

  //
  // The particle loop
  //
  for (int i=nbeg; i<nend; i++) {

    int n = c->levlist[level][i];
    Particle *p = c->Part(n);

    if (DTold) {

      // dtv = eps* r/v         -- roughly, crossing time
      // dta = eps* v/a         -- force scale
      // dtA = eps* sqrt(r/a)   -- acceleration time

      rtot = 0.0;
      vtot = 0.0;
      atot = 0.0;
      
      for (int k=0; k<c->dim; k++) {
	rtot += 
	  c->Pos(n, k, Component::Local | Component::Centered) *
	  c->Pos(n, k, Component::Local | Component::Centered) ;
	vtot += p->vel[k]*p->vel[k];
	atot += p->acc[k]*p->acc[k];
      }
      rtot = sqrt(rtot);
      vtot = sqrt(vtot) + 1.0e-18;
      atot = sqrt(atot) + 1.0e-18;

      dsr = p->scale;
      if (dsr>0) dts = dynfracS*dsr/vtot;
      else       dts = 1.0/eps;

      dtv = dynfracV*rtot/vtot;
      dta = dynfracA*vtot/atot;
      dtA = dynfracP*sqrt(rtot/atot);

    } else {

      // dtd = eps* rscale/v_i    -- char. drift time scale
      // dtv = eps* min(v_i/a_i)  -- char. force time scale
      // dta = eps* phi/(v * a)   -- char. work time scale
      // dtA = eps* sqrt(phi/a^2) -- char. "escape" time scale


      dtr  = 0.0;
      vtot = 0.0;
      atot = 0.0;

      for (int k=0; k<c->dim; k++) {
	dtr  += p->vel[k]*p->acc[k];
	vtot += p->vel[k]*p->vel[k];
	atot += p->acc[k]*p->acc[k];
      }
      ptot = fabs(p->pot + p->potext);


      dsr = p->scale;
      if (dsr>0) dts = dynfracS*dsr/fabs(sqrt(vtot)+eps);
      else       dts = 1.0/eps;
      
      dtd = dynfracD * 1.0/sqrt(vtot+eps);
      dtv = dynfracV * sqrt(vtot/(atot+eps));
      dta = dynfracA * ptot/(fabs(dtr)+eps);
      dtA = dynfracP * sqrt(ptot/(atot*atot+eps));
    }

    map<double, int> dseq;

    if (DTold) {
      dseq[dtv] = 0;
      dseq[dts] = 1;
      if ( dta > 0.0 ) dseq[dta] = 2;
      if ( dtA > 0.0 ) dseq[dtA] = 3;
      if ( (dtr=p->dtreq) > 0.0 ) dseq[dtr] = 4;
    } else {
      dseq[dtd] = 0;
      dseq[dtv] = 1;
      dseq[dts] = 2;
      if ( dta > 0.0 ) dseq[dta] = 3;
      if ( dtA > 0.0 ) dseq[dtA] = 4;
      if ( (dtr=p->dtreq) > 0.0 ) dseq[dtr] = 5;
    }



    // Smallest time step
    //
    dt = dseq.begin()->first;

    if (mstep == 0) {
      //
      // Tally smallest (e.g. controlling) timestep
      //
      tmdt[id][level][dseq.begin()->second]++;
      //
      // Counter
      //
      tmdt[id][level][mdtDim-1]++;
    }

    // Time step wants to be LARGER than the maximum
    if (dt>dtime) {
      lev = 0;
      maxdt1[id] = max<double>(dt, maxdt1[id]);
      offhi++;
    }
    else lev = (int)floor(log(dtime/dt)/log(2.0));
    
    // Time step wants to be SMALLER than the maximum
    if (lev>multistep) {
      lev = multistep;
      mindt1[id] = min<double>(dt, mindt1[id]);
      offlo++;
    }
    
    unsigned plev = p->level;
    unsigned nlev = lev;

    // Sanity check
    if (level != plev) {
      cerr << "Process " << myid << " id=" << id 
	   << ": n=" << n << " level=" << level 
	   << " tlevel=" << plev
	   << " plevel=" << p->level << endl;
    }

    if (nlev != level) {
      //
      // Update coefficients
      //
      c->force->multistep_update(plev, nlev, c, n, id);
      p->level = lev;
    }
  }

  if (VERBOSE>0) {
    offlo1[c][id] += offlo;
    offhi1[c][id] += offhi;
  }
  
  return (NULL);
}


void adjust_multistep_level(bool all)
{
  if (!multistep) return;

  // FOR DEBUGGING
  // if (mstep!=0) return;
  // END DEBUGGING

  //
  // Begin the update
  //
  for (auto c : comp.components) c->force->multistep_update_begin();

  //
  // Preliminary data structure and thread creation
  //
  mindt1 = vector< double > (nthrds,  1.0e20);
  maxdt1 = vector< double > (nthrds, -1.0e20);

  if (VERBOSE>0) {

    if (offhi1.size()==0 || mstep==0) {

      for (auto c : comp.components) {
	for (int n=0; n<nthrds; n++) {
	  offhi1[c] = vector<unsigned>(nthrds, 0);
	  offlo1[c] = vector<unsigned>(nthrds, 0);
	}
      }
    }
  }

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
  

  if (tmdt.size() == 0) {
    tmdt = vector< vector< vector<unsigned> > >(nthrds);
    for (int n=0; n<nthrds; n++) {
      tmdt[n] = vector< vector<unsigned> >(multistep+1);
      for (int k=0; k<=multistep; k++) tmdt[n][k] = vector<unsigned>(mdtDim);
    }
  }

  for (auto c : comp.components) {
    
    if (mstep == 0) {
      for (int n=0; n<nthrds; n++)
	for (int k=0; k<=multistep; k++) 
	  for (int j=0; j<mdtDim; j++) tmdt[n][k][j] = 0;
    }
    
    for (int level=0; level<=multistep; level++) {
      
      if (all || mactive[mstep][level]) {

	if (nthrds==1) {

	  td[0].level = level;
	  td[0].id = 0;
	  td[0].c = c;

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
	    td[i].c = c;
	    
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

    if (mstep == 0) {
      for (int n=0; n<nthrds; n++)
	for (int k=0; k<=multistep; k++) 
	  for (int j=0; j<mdtDim; j++) 
	    c->mdt_ctr[k][j] += tmdt[n][k][j];
    }
  }

  delete [] td;
  delete [] t;

  //
  // Finish the update
  //
  for (auto c : comp.components) {
    c->reset_level_lists();
    c->fix_positions();
    c->force->multistep_update_finish();
  }


  //
  // Diagnostic output
  //
  if (VERBOSE>0 && mstep==0) {

    //
    // Count offgrid particles in the threads
    //
    map< Component*, unsigned > offlo, offhi;

    for (int n=1; n<nthrds; n++) {
      mindt1[0] = min<double>(mindt1[0], mindt1[n]);
      maxdt1[0] = max<double>(maxdt1[0], maxdt1[n]);
    }

    double mindt, maxdt;

    MPI_Reduce(&mindt1[0], &mindt, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&maxdt1[0], &maxdt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    
    for (auto c : comp.components) {

      unsigned cofflo = 0, coffhi = 0;
      for (int n=0; n<nthrds; n++) {
	cofflo += offlo1[c][n];
	coffhi += offhi1[c][n];
      }

      MPI_Reduce(&cofflo, &offlo[c], 1,
		 MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
      
      MPI_Reduce(&coffhi, &offhi[c], 1,
		 MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    
    if (myid==0) {
      
      unsigned sumlo=0, sumhi=0;
      for (auto c : comp.components) {
	sumlo += offlo[c];
	sumhi += offhi[c];
      }
      
      if (sumlo || sumhi) {
	cout << endl
	     << setw(70) << setfill('-') << '-' << endl << setfill(' ')
	     << left << "--- Multistepping overrun" << endl;
	if (sumlo)
	  cout << left << "--- Min DT=" << setw(16) << mindt  
	       << " < " << setw(16) << dtime/(1<<multistep) 
	       << " [" << sumlo << "]" << endl;
	if (sumhi)
	  cout << left << "--- Max DT=" << setw(16) << maxdt  
	       << " > " << setw(16) << dtime 
	       << " [" << sumhi << "]" << endl;
	cout << setw(70) << setfill('-') << '-' << endl 
	     << setfill(' ') << right;

	if (sumlo) {
	  for (auto c : comp.components) {
	    ostringstream sout;
	    sout << "Component <" << c->name << ">";
	    cout << setw(30) << sout.str() << " |   low: "
		 << offlo[c] << "/" << c->nbodies_tot << endl;
	  }
	}

	if (sumhi) {
	  for (auto c : comp.components) {
	    ostringstream sout;
	    sout << "Component <" << c->name << ">";
	    cout << setw(30) << sout.str() << " |  high: "
		 << offhi[c] << "/" << c->nbodies_tot << endl;
	  }
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
  mactive.push_back(vector<bool>(multistep+1, true));
  for (int ms=1; ms<=Mstep; ms++)
    mactive.push_back(vector<bool>(multistep+1, false));

				// Find and save the active levels at each step
  for (int ms=1; ms<=Mstep; ms++) {
    for (unsigned mlevel=0; mlevel<=multistep; mlevel++) {
      if ( (ms % (1 << (multistep-mlevel))) == 0) 
	mactive[ms][mlevel] = true;
    }
  }

  mfirst = vector<int>(Mstep+1);
				// Lowest active level at each step
  for (int ms=0; ms<=Mstep; ms++) {
    for (int mlevel=0; mlevel<=multistep; mlevel++) {
      if (mactive[ms][mlevel]) {
	mfirst[ms] = mlevel;
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
    for (int ms=0; ms<=Mstep; ms++) {
      cout << setw(4) << ms << "/";
      cout << setw(3) << mfirst[ms] << ": ";
      for (int mlevel=0; mlevel<=multistep; mlevel++)
	if (mactive[ms][mlevel]) cout << setw(3) << "+";
	else                     cout << setw(3) << "-";
      cout << endl;
    }
    cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
  }

}
