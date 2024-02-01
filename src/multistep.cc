/* 
  Multi-timestepping support routines
*/

#include <expand.H>
#include <sstream>
#include <chrono>
#include <limits>
#include <map>

// #define VERBOSE_TIMING

// Cuda routines
//
#if HAVE_LIBCUDA==1
extern void cuda_initialize_multistep();
extern void cuda_compute_levels();
#endif

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
static vector< double > adjtm1;
static vector< double > adjtm2;
static vector< unsigned> numsw;
static vector< unsigned> numtt;

// Type counter
static vector< vector< vector<unsigned> > > tmdt;

//
// The threaded routine
//
void * adjust_multistep_level_thread(void *ptr)
{
  // Begin diagnostic timing
  std::chrono::high_resolution_clock::time_point start0, finish0;

  start0 = std::chrono::high_resolution_clock::now();

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
  int offlo = 0, offhi = 0;

  //
  // Compute the beginning and end points for threads
  //
  int nbeg = npart*id    /nthrds;
  int nend = npart*(id+1)/nthrds;

  //
  // Small positive constant
  //
  const double eps = 1.0e-10;

  //
  // The particle loop
  //
  for (int i=nbeg; i<nend; i++) {

    int n = c->levlist[level][i];
    Particle *p = c->Part(n);

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
    dtA = dynfracP * sqrt(ptot/(atot+eps));

    std::map<double, int> dseq;

    dseq[dtd] = 0;
    dseq[dtv] = 1;
    dseq[dts] = 2;
    if ( dta > 0.0 ) dseq[dta] = 3;
    if ( dtA > 0.0 ) dseq[dtA] = 4;

    // Smallest time step
    //
    dt = std::max<double>(eps, dseq.begin()->first);

    // Enforce minimum step per level
    //
    bool firstCall = this_step==0 and mdrft==0;

    if (c->NoSwitch()) {
      if ((c->DTreset() and mstep==0) or firstCall)
	p->dtreq = std::numeric_limits<double>::max();
      if (dt < p->dtreq)
	p->dtreq = dt;
    } else {
      p->dtreq = dt;
    }

    // Update coefficients at this substep?
    //
    bool apply = not c->NoSwitch() or mdrft==Mstep or firstCall;
    //                        ^             ^              ^
    //                        |             |              |
    // at every substep-------+             |              |
    //                                      |              |
    // otherwise: at end of full step-------+              |
    //                                                     |
    // or on the very first call to initialize levels------+

    // Only assign levels on first call; option for testing
    //
    if (not firstCall and c->FreezeLev()) apply = false;

    // Select this substep for update?
    //
    if (apply) {

      unsigned plev = p->level;
      unsigned nlev = plev;

      // Time step wants to be LARGER than the maximum
      //
      if (p->dtreq>dtime) {
	nlev = 0;
	maxdt1[id] = std::max<double>(p->dtreq, maxdt1[id]);
	offhi++;
      }
      else nlev = (int)floor(log(dtime/p->dtreq)/log(2.0));

      // Enforce n-level shifts at a time
      //
      if (shiftlevl) {
	if (nlev > plev) {
	  if (nlev - plev > shiftlevl) nlev = plev + shiftlevl;
	} else if (plev > nlev) {
	  if (plev - nlev > shiftlevl) nlev = plev - shiftlevl;
	}
      }
      
      // Time step wants to be SMALLER than the maximum
      //
      if (nlev>multistep) {
	nlev = multistep;
	mindt1[id] = std::min<double>(p->dtreq, mindt1[id]);
	offlo++;
      }
      
      // Limit new level to minimum active level
      //
      nlev = std::max<int>(nlev, mfirst[mdrft]);

      if (plev != nlev) {
	std::chrono::high_resolution_clock::time_point start1, finish1;
	start1 = std::chrono::high_resolution_clock::now();
	c->force->multistep_update(plev, nlev, c, n, id);
	finish1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::micro> duration = finish1 - start1;
	adjtm2[id] += duration.count();
	p->level = nlev;
	numsw[id]++;
      }
      numtt[id]++;
    }

    // For reporting level populations: evaluating at the final sub
    // step guarantees that every particle is active.  Also note:
    // mdrft equals Mstep for multistep=0; so multistep=0 gives the
    // desired evaluation at every step.
    //
    if (mdrft == Mstep) {
      //
      // Tally smallest (e.g. controlling) timestep
      //
      tmdt[id][p->level][dseq.begin()->second]++;
      //
      // Counter
      //
      tmdt[id][p->level][mdtDim-1]++;
    }
  }

  if (VERBOSE>0) {
    offlo1[c][id] += offlo;
    offhi1[c][id] += offhi;
  }
  
  finish0 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::micro> duration = finish0 - start0;
  adjtm1[id] += duration.count();

  return (NULL);
}


void adjust_multistep_level()
{
  if (!multistep) return;

#ifdef VERBOSE_TIMING
  std::cout << "[" << myid << "] ENTERING adjust multistep level"
	    << std::endl;
  auto dbg_start = std::chrono::high_resolution_clock::now();
#endif

#if HAVE_LIBCUDA==1
  if (use_cuda) {

    cuda_compute_levels();

#ifdef VERBOSE_TIMING
    auto dbg_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> dbg_adjust = dbg_finish - dbg_start;

    std::cout << "[" << myid << "] LEAVING adjust multistep level: "
	      << dbg_adjust.count()*1.0e-6 << std::endl;
#endif    
    return;
  }
#endif

  // Begin diagnostic timing
  std::chrono::high_resolution_clock::time_point start, finish;

  start = std::chrono::high_resolution_clock::now();

  // Always assign levels on first pass
  //
  bool firstCall = this_step==0 and mdrft==0;

  //
  // Begin the update
  //
  for (auto c : comp->components) {
    c->force->multistep_update_begin();
  }

  //
  // Preliminary data structure and thread creation
  //
  mindt1 = vector< double > (nthrds,  1.0e20);
  maxdt1 = vector< double > (nthrds, -1.0e20);
  adjtm1 = vector< double > (nthrds, 0.0);
  adjtm2 = vector< double > (nthrds, 0.0);
  numsw  = vector< unsigned > (nthrds, 0);
  numtt  = vector< unsigned > (nthrds, 0);

  if (VERBOSE>0) {

    if (offhi1.size()==0 || mdrft==Mstep) {

      for (auto c : comp->components) {
	for (int n=0; n<nthrds; n++) {
	  offhi1[c] = vector<unsigned>(nthrds, 0);
	  offlo1[c] = vector<unsigned>(nthrds, 0);
	}
      }
    }
  }

  thrd_pass_sync* td = new thrd_pass_sync [nthrds];

  if (!td) {
    std::ostringstream sout;
    sout << "Process " << myid
	 << ": adjust_multistep_level: error allocating thread structures";
    throw GenericError(sout.str(), __FILE__, __LINE__, 1024, true);
  }

  pthread_t* t  = new pthread_t [nthrds];

  if (!t) {
    std::ostringstream sout;
    sout << "Process " << myid
	 << ": adjust_multistep_level: error allocating memory for thread";
    throw GenericError(sout.str(), __FILE__, __LINE__, 1024, true);
  }
  

  if (tmdt.size() == 0) {
    tmdt = vector< vector< vector<unsigned> > >(nthrds);
    for (int n=0; n<nthrds; n++) {
      tmdt[n] = vector< vector<unsigned> >(multistep+1);
      for (int k=0; k<=multistep; k++) tmdt[n][k] = vector<unsigned>(mdtDim);
    }
  }

  for (auto c : comp->components) {
    
    // For reporting level populations: evaluating at the final sub
    // step guarantees that every particle is active.  Also note:
    // mdrft equals Mstep for multistep=0; so multistep=0 gives the
    // desired evaluation at every step.
    //
    if (mdrft == Mstep) {
      for (int n=0; n<nthrds; n++)
	for (int k=0; k<=multistep; k++) 
	  for (int j=0; j<mdtDim; j++) tmdt[n][k][j] = 0;
    }
    
    int first = mfirst[mdrft];	// First active level at drifted
				// subgrid position

    if (this_step==0 and mstep==0) first = 0; // Do all levels

    for (int level=first; level<=multistep; level++) {
      
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
	    std::ostringstream sout;
	    sout << "Process " << myid
		 << " adjust_multistep_level: cannot make thread " << i
		 << ", errcode=" << errcode;
	    throw GenericError(sout.str(), __FILE__, __LINE__, 1024, true);
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
	    std::ostringstream sout;
	    sout << "Process " << myid
		 << " adjust_multistep_level: thread join " << i
		 << " failed, errcode=" << errcode;
	    throw GenericError(sout.str(), __FILE__, __LINE__, 1024, true);
	  }
#ifdef DEBUG    
	  cout << "Process " << myid 
	       << ": multistep thread <" << i << "> thread exited\n";
#endif
	}
      }
    }

    // Accumulate counters for all threads at the master step boundary
    //
    if (mdrft == Mstep) {
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
  for (auto c : comp->components) {

    // Assign levels on every step or final drift for noswitch==True
    //
    bool apply = not c->NoSwitch() or mdrft==Mstep or firstCall;

    // Only apply on first pass (for testing)
    //
    if (not firstCall and c->FreezeLev()) apply = false;

    // Update level coefficients
    //
    if (apply) {
      c->force->multistep_update_finish();
      c->reset_level_lists();
    }
    
    c->fix_positions();
  }


  finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::micro> update_t = finish - start;
  start =  std::chrono::high_resolution_clock::now();

  //
  // Diagnostic output
  //
  if (mdrft==Mstep) {

    //
    // Count offgrid particles in the threads
    //
    map< Component*, unsigned > offlo, offhi;

    for (int n=1; n<nthrds; n++) {
      mindt1[0] = min<double>(mindt1[0], mindt1[n]);
      maxdt1[0] = max<double>(maxdt1[0], maxdt1[n]);
      adjtm1[0] += adjtm1[n];
      adjtm2[0] += adjtm2[n];
      numsw [0] += numsw [n];
      numtt [0] += numtt [n];
    }

    double mindt, maxdt, atim1, atim2;
    unsigned numtot, numadj;

    MPI_Reduce(&mindt1[0], &mindt,  1, MPI_DOUBLE,   MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&maxdt1[0], &maxdt,  1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&adjtm1[0], &atim1,  1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&adjtm2[0], &atim2,  1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&adjtm1[0], &atim1,  1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&numsw [0], &numadj, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&numtt [0], &numtot, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);

    
    for (auto c : comp->components) {

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

    
    finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> collate_t = finish - start;

    if (myid==0) {
      
      if (VERBOSE>3 and atim1>0) {
	auto pc = std::cout.precision(1);
	std::cout << endl
		  << std::setw(70) << std::setfill('-') << '-' << std::endl << std::setfill(' ')
		  << std::left << "--- Coefficient adjust stats"  << std::endl << std::fixed
		  << std::left << "--- Coef/DT = " << 100.0*atim2/atim1   << "%" << std::endl
		  << std::left << "--- Adj/Tot = " << 100.0*numadj/numtot << "%" << std::endl << std::scientific
		  << std::left << "--- Update  = " << std::setprecision(4) << update_t.count() *1.0e-6 << " sec" << std::endl
		  << std::left << "--- Collate = " << std::setprecision(4) << collate_t.count()*1.0e-6 << " sec" << std::endl
		  << std::setw(70) << std::setfill('-') << '-' << std::endl << std::setfill(' ');
	std::cout.precision(pc);
      }

      if (VERBOSE>0) {

	unsigned sumlo=0, sumhi=0;
	for (auto c : comp->components) {
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
	    for (auto c : comp->components) {
	      ostringstream sout;
	      sout << "Component <" << c->name << ">";
	      cout << setw(30) << sout.str() << " |   low: "
		   << offlo[c] << "/" << c->CurTotal() << endl;
	    }
	  }
	  
	  if (sumhi) {
	    for (auto c : comp->components) {
	      ostringstream sout;
	      sout << "Component <" << c->name << ">";
	      cout << setw(30) << sout.str() << " |  high: "
		   << offhi[c] << "/" << c->CurTotal() << endl;
	    }
	  }
	  
	  cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
	}
      }

      // Impose sanity checks and trigger stop if necessary
      //
      // Search for components that exceed max level threshold
      std::set<Component*> bad;
      for (auto c : comp->components) {
	double frac = static_cast<double>(offlo[c])/c->CurTotal();
	if (frac > max_mindt) bad.insert(c);
      }

      // If we found any bad components, set signal and print diagnostics
      //
      if (bad.size()) {
	// Set flag to stop at the end of the current step
	//
	quit_signal = 1;
	
	// Log info to stdout
	//
	std::cout << std::setw(70) << std::setfill('-') << '-'
		  << std::endl << std::setfill(' ')
		  << "---- EXP is going to stop this run for you at the end of this step" << std::endl
		  << "---- because these components have more than "
		  << floor(100.0*max_mindt)
		  << "% particles below the minimum time step" << std::endl;
	std::cout << std::setw(70) << std::setfill('-') << '-'
		  << std::endl << std::setfill(' ');
	for (auto c : bad) {
	  std::ostringstream sout;
	  sout << "Component <" << c->name << ">";
	  std::cout << std::setw(30) << sout.str() << " | "
		    << offlo[c] << "/" << c->CurTotal() << std::endl;
	}
	std::cout << std::setw(70) << std::setfill('-') << '-'
		  << std::endl << std::setfill(' ')
		  << "---- Try decreasing your 'dtime' value, increasing your 'multilevel' value, or both!" << std::endl
		  << std::setw(70) << std::setfill('-') << '-'
		  << std::endl << std::setfill(' ');
      }
    }
  }

#ifdef VERBOSE_TIMING
  auto dbg_finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::micro> dbg_adjust = dbg_finish - dbg_start;

  std::cout << "LEAVING adjust multistep level ["
	    << dbg_adjust.count()*1.0e-6 << "]" << std::endl;
#endif
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

				// Compute interpolation arrays
  std::vector<int> dstep(multistep+1);
  dstepL.resize(multistep+1);
  dstepN.resize(multistep+1);

  for (int ms=0; ms<=multistep; ms++) {
    dstep[ms] = 1<<ms;
    dstepL[ms].resize(Mstep, 0);
    dstepN[ms].resize(Mstep, 0);
  }

  for (int ms=0; ms<=multistep; ms++) {
    int rev = multistep - ms;
    for (int n=0; n<Mstep; n++) {
      dstepL[rev][n] = (n/dstep[ms])*dstep[ms];
      dstepN[rev][n] = dstepL[rev][n] + dstep[ms];
    }
  }


  if (VERBOSE>10 && myid==0 && multistep) {
    std::cout << std::setw(70) << std::setfill('-') << '-' << std::endl
	      << std::setw(70) << std::left << "--- Multistep level interpolation intervals " << std::endl
	      << std::setw(70) << setfill('-') << '-' << std::endl << std::setfill(' ') << std::right;

    for (int l=0; l<=multistep; l++) {
      for (int n=0; n<Mstep; n++) {
	std::cout << std::setw(4) << l
		  << std::setw(4) << n
		  << std::setw(8) << dstepL[l][n]
		  << std::setw(8) << dstepN[l][n]
		  << std::endl;
      }
    }
    std::cout << std::endl;

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

#if HAVE_LIBCUDA==1
  if (use_cuda) cuda_initialize_multistep();
#endif

}
